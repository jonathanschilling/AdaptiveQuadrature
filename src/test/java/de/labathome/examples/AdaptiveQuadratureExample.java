package de.labathome.examples;

import java.util.function.UnaryOperator;

import de.labathome.AdaptiveQuadrature;

public class AdaptiveQuadratureExample {

	public static void main(String[] args) {

		// integrand: parabola 3*x^2
		UnaryOperator<double[]> parabola = x -> {
			int n = x.length;
			double[] f = new double[n];
			for (int i = 0; i < n; ++i) {
				f[i] = 3.0 * x[i] * x[i];
			}
			return f;
		};

		double a = 0.0; // lower end of integration interval
		double b = 1.0; // upper end of integration interval
		double relTol = 1.0e-6;
		double absTol = 1.0e-6;
		int maxEval = 0; // unlimited
		double[] result = AdaptiveQuadrature.integrate(parabola, a, b, relTol, absTol, maxEval);
		System.out.printf("\\int_0^1 3*x^2 dx = %.3e +/- %.3e\n", result[0], result[1]);
	}

}
