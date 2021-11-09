package de.labathome.examples;

import java.util.function.UnaryOperator;

import de.labathome.AdaptiveQuadrature;

public class AdaptiveQuadratureExample {

	public static void main(String[] args) {

		UnaryOperator<double[]> parabola = x -> {
			int n = x.length;
			double[] f = new double[n];
			for (int i=0; i<n; ++i) {
				f[i] = 3.0*x[i]*x[i];
			}
			return f;
		};

		double[] result = AdaptiveQuadrature.integrate(parabola, 0.0, 1.0, 1.0e-6, 1.0e-6, 0);
		System.out.printf("\\int_0^1 3*x^2 dx = %.3e +/- %.3e\n", result[0], result[1]);
	}

}
