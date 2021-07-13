package de.labathome;

import java.util.ArrayList;
import java.util.List;

/**
 * Adaptive one-dimensional Gauss-Kronrod quadrature.<br>
 * <br>
 * A 7-point Gauss-Legendre quadrature is used to estimate the integral.
 * A 15-point Gauss-Kronrod quadrature using the function values from the 7-point Gauss-Legendre integral
 * is used to check for convergence.
 * The domain over which to integrate is split in half if the convergence criteria are not met. <br>
 * <br>
 * If one of the bounds is infinite, a change of variables is performed according to:
 * <pre>
 * \int_{a}^{+\infty} f(x) dx = \int_0^1 f(a+t/(1-t)) 1/(1-t)^2 dt .
 * </pre>
 * If both bounds are infinite, a different change of variables is performed according to:
 * <pre>
 * \int_{-\infty}^{+\infty} f(x) dx = \int_{-1}^1 f(t/(1-t^2)) (1 + t^2)/(1-t^2)^2 dt .
 * </pre>
 * Note that above change of variables might introduce quite some nonlinearity into the integrand
 * and it is expected that indefinite integrals require more function evaluations to converge.
 *
 * @author Jonathan Schilling (<a href="mailto:jonathan.schilling@mail.de">jonathan.schilling@mail.de</a>)
 */
public class AdaptiveQuadrature {

	/**
	 * Integrate the given integrand from {@code lowerBound} to {@code upperBound}.
	 * The numerical quadrature will be successively refined until either a given absolute tolerance
	 * or a prescribed relative tolerance is fulfilled.
	 * Note that {@code relTol} and {@code absTol} must not be both NaN.
	 *
	 * @param integrand function to integrate
	 * @param lowerBound lower boundary of integral; can be {@code Double.NEGATIVE_INFINITY}
	 * @param upperBound upper boundary of integral; can be {@code Double.POSITIVE_INFINITY}
	 * @param relTol relative tolerance; set to {@code Double.NaN} to disable
	 * @param absTol absolute tolerance; set to {@code Double.NaN} to disable
	 * @param maxEval maximum number of function evaluations; set to 0 to disable limit on number of function evaluations
	 */
	public static double[] integrate(Integrand integrand, double lowerBound, double upperBound, double relTol, double absTol, int maxEval) {

		if (Double.isNaN(lowerBound)) {
			throw new RuntimeException("lower bound must not be NaN");
		}

		if (Double.isNaN(upperBound)) {
			throw new RuntimeException("upper bound must not be NaN");
		}

		if (lowerBound == upperBound) {
			if (Double.isInfinite(lowerBound)) {
				// if bounds are equal, they must not be both negative infinity or both positive infinity
				throw new RuntimeException("integral from "+lowerBound+" to "+upperBound+" is not supported");
			}

			// otherwise quick return: equal lower and upper bound mean an integral of zero
			return new double[] {0.0, 0.0};
		}

		// invert integral value if bounds are swapped
		final boolean invert = (lowerBound > upperBound);
		if (invert) {
			final double temp = upperBound;
			upperBound = lowerBound;
			lowerBound = temp;
		}

		// apply re-scaling of integrand if one or both integration bounds are infinite
		RescaledIntegrand rescaledIntegrand = new RescaledIntegrand(integrand, lowerBound, upperBound);
		final double scaledLowerBound = rescaledIntegrand.getScaledLowerBound();
		final double scaledUpperBound = rescaledIntegrand.getScaledUpperBound();


		double center = (scaledUpperBound+scaledLowerBound)/2.0;
		double halfWidth = (scaledUpperBound-scaledLowerBound)/2.0;
		Interval rootInterval = new Interval(center, halfWidth);

		List<Interval> intervalsToEval = new ArrayList<>();
		intervalsToEval.add(rootInterval);

		GaussKronrod.evalGaussKronrod(rescaledIntegrand, intervalsToEval);

		IntervalHeap intervalHeap = new IntervalHeap();
		intervalHeap.add(rootInterval);

		int numEval = 15;
		boolean converged = true;
		double result = intervalHeap.getIntegralValue();
		double errorEstimate = intervalHeap.getErrorEstimate();

		// adaptively refine intervals until convergence
		while (maxEval == 0 || numEval < maxEval) {
			result = intervalHeap.getIntegralValue();
			errorEstimate = intervalHeap.getErrorEstimate();

			converged = false;
			if (!Double.isNaN(absTol)) {
				converged |= errorEstimate < absTol;
			}
			if (!Double.isNaN(relTol)) {
				converged |= errorEstimate < Math.abs(result)*relTol;
			}

			if (converged) {
				System.out.println("converged after "+numEval+" function evaluations");
				break;
			}

			do {
				Interval worstInterval = intervalHeap.poll();
				errorEstimate -= worstInterval.getErrorEstimate();

				Interval worstSplit = worstInterval.cutInHalf();
				intervalsToEval.add(worstInterval); // re-add halvened original interval
				intervalsToEval.add(worstSplit); // add other half of previous interval

				numEval += 2*GaussKronrod.NUM_KRONROD_POINTS;

				converged = false;
				if (!Double.isNaN(absTol)) {
					converged |= errorEstimate < absTol;
				}
				if (!Double.isNaN(relTol)) {
					converged |= errorEstimate < Math.abs(result)*relTol;
				}

				if (converged) {
					break;
				}

			} while (intervalHeap.size() > 0 && (numEval < maxEval || maxEval==0));

			GaussKronrod.evalGaussKronrod(rescaledIntegrand, intervalsToEval);

			for (Interval interval: intervalsToEval) {
				intervalHeap.add(interval);
			}
		}

		if (!converged) {
			System.out.println("Cubature did not converge after "+numEval+" function evaluations!");
		}

		result = 0.0;
		errorEstimate = 0.0;
		for (Interval interval: intervalHeap) {
			result += interval.getIntegralValue();
			errorEstimate += interval.getErrorEstimate();
		}

		if (invert) {
			return new double[] {-result, errorEstimate};
		} else {
			return new double[] {result, errorEstimate};
		}
	}


	public static void main(String[] args) {

		double a = Double.POSITIVE_INFINITY;
		double b = Double.POSITIVE_INFINITY;
		System.out.println(a == b);
		System.out.println(Double.isInfinite(a));

		a = Double.NEGATIVE_INFINITY;
		b = Double.NEGATIVE_INFINITY;
		System.out.println(a == b);
		System.out.println(Double.isInfinite(a));

		class Parabola implements Integrand {
			@Override
			public double[] eval(double... x) {
				int n = x.length;
				double[] f = new double[n];
				for (int i=0; i<n; ++i) {
					f[i] = 3.0*x[i]*x[i];
				}
				return f;
			}
		};

		Parabola integrand = new Parabola();

		double[] result = AdaptiveQuadrature.integrate(integrand, 0.0, 1.0, 1.0e-6, 1.0e-6, 0);
		System.out.printf("integral = %.3e +/- %.3e\n", result[0], result[1]);
	}


}
