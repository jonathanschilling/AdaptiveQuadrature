package de.labathome.examples;

import de.labathome.AdaptiveQuadrature;
import de.labathome.Integrand;

public class AdaptiveQuadratureExample {

	public static void main(String[] args) {

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
		System.out.printf("1/3 int_0^1 x^2 dx = %.3e +/- %.3e\n", result[0], result[1]);
	}

}
