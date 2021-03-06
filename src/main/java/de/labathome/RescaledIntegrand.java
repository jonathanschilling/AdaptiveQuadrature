package de.labathome;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.UnaryOperator;

public class RescaledIntegrand implements UnaryOperator<double[]> {

	private UnaryOperator<double[]> integrand;
	private AtomicBoolean gracefulStop;

	private double lowerBound;
	private double upperBound;

	private double scaledLowerBound;
	private double scaledUpperBound;

	public RescaledIntegrand(UnaryOperator<double[]> integrand, double lowerBound, double upperBound, AtomicBoolean gracefulStop) {
		this.integrand = integrand;
		this.gracefulStop = gracefulStop;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;

		if (Double.isFinite(lowerBound) && Double.isFinite(upperBound)) {
			scaledLowerBound = lowerBound;
			scaledUpperBound = upperBound;
		} else if (Double.isFinite(lowerBound) && Double.isInfinite(upperBound)
				|| Double.isInfinite(lowerBound) && Double.isFinite(upperBound)) {
			scaledLowerBound = 0.0;
			scaledUpperBound = 1.0;
		} else if (Double.isInfinite(lowerBound) && Double.isInfinite(upperBound)) {
			scaledLowerBound = -1.0;
			scaledUpperBound =  1.0;
		}
	}

	public double getScaledLowerBound() {
		return scaledLowerBound;
	}

	public double getScaledUpperBound() {
		return scaledUpperBound;
	}

	@Override
	public double[] apply(double[] x) {
		if (Double.isFinite(lowerBound) && Double.isFinite(upperBound)) {
			// quick return: if no rescaling is needed, evaluate directly
			return integrand.apply(x);
		}

		final int nPoints = x.length;
		final double[] rescaledX = new double[nPoints];
		final double[] jacobianFactors = new double[nPoints];

		if (Double.isFinite(lowerBound) && Double.isInfinite(upperBound)) {
			// \int_a^{+\infty} f(x) dx = \int_0^1 f(a+t/(1-t)) 1/(1-t)^2 dt
			for (int i=0; i<nPoints; ++i) {
				rescaledX[i] = lowerBound + x[i]/(1.0 - x[i]);
				jacobianFactors[i] = 1.0/((1.0 - x[i])*(1.0 - x[i]));
			}
		} else if (Double.isInfinite(lowerBound) && Double.isFinite(upperBound)) {
			// \int_{-\infty}^b f(x) dx = \int_0^1 f(b-t/(1-t)) 1/(1-t)^2 dt
			for (int i=0; i<nPoints; ++i) {
				rescaledX[i] = upperBound - x[i]/(1.0 - x[i]);
				jacobianFactors[i] = 1.0/((1.0 - x[i])*(1.0 - x[i]));
			}
		} else {
			// \int_{-\infty}^{+\infty} f(x) dx = \int_{-1}^1 f(t/(1-t^2)) (1+t^2)/(1-t^2)^2 dt
			for (int i=0; i<nPoints; ++i) {
				rescaledX[i] = x[i]/(1.0 - x[i]*x[i]);
				jacobianFactors[i] = (1.0 + x[i]*x[i])/((1.0 - x[i]*x[i])*(1.0 - x[i]*x[i]));
			}
		}

		final double[] evalIntegrand = integrand.apply(rescaledX);
		if (!gracefulStop.get()) {
			for (int i=0; i<nPoints; ++i) {
				evalIntegrand[i] *= jacobianFactors[i];
			}
			return evalIntegrand;
		} else {
			// if gracefulStop feature cancels the integration, the integrand might return invalid results
			// --> skip using the returned value in any way
			return null;
		}
	}

}
