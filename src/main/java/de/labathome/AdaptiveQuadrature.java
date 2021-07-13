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
	 */
	public static double[] integrate(Integrand integrand, double lowerBound, double upperBound, double relTol, double absTol) {

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

			// otherwise quick return:
			// equal lower and upper bound mean an integral of zero
			return new double[] {0.0, 0.0};
		}

		// invert integral value if bounds are swapped
		boolean invert = (lowerBound > upperBound);
		if (invert) {
			double temp = upperBound;
			upperBound = lowerBound;
			lowerBound = temp;
		}

		// apply re-scaling of integrand if one or both integration bounds are infinite
		RescaledIntegrand rescaledIntegrand = new RescaledIntegrand(integrand, lowerBound, upperBound);
		final double scaledLowerBound = rescaledIntegrand.getScaledLowerBound();
		final double scaledUpperBound = rescaledIntegrand.getScaledUpperBound();

		double result = 0.0;
		double errorEstimate = Double.POSITIVE_INFINITY;

		double center = (scaledUpperBound+scaledLowerBound)/2.0;
		double halfWidth = (scaledUpperBound-scaledLowerBound)/2.0;

		double fCenter = integrand.eval(center)[0];
		if (!Double.isFinite(fCenter)) {
			throw new RuntimeException("Evaluation of integrand at interval center failed: got "+fCenter);
		}

		Interval rootInterval = new Interval(center, halfWidth);

		List<Interval> intervalsToEval = new ArrayList<>();
		intervalsToEval.add(rootInterval);

		double[][] gkResult = GaussKronrod.evalGaussKronrod(rescaledIntegrand, intervalsToEval);


		if (invert) {
			return new double[] {-result, errorEstimate};
		} else {
			return new double[] {result, errorEstimate};
		}
	}


	public static void main(String[] args) {

		class Parabola implements Integrand {
			@Override
			public double[] eval(double... x) {
				int n = x.length;
				double[] f = new double[n];
				for (int i=0; i<n; ++i) {
					f[i] = 0.5*x[i]*x[i];
				}
				return f;
			}
		};

		Parabola integrand = new Parabola();

		double[] result = AdaptiveQuadrature.integrate(integrand, 0.0, 1.0, 1.0e-6, 1.0e-6);
	}


}
