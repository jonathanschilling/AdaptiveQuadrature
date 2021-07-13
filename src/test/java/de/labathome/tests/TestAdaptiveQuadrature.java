package de.labathome.tests;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

import de.labathome.AdaptiveQuadrature;
import de.labathome.Integrand;

public class TestAdaptiveQuadrature {

	@Test
	public void testParabola() {

		final double tolerance = 1.0e-12;

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

		double[] result = AdaptiveQuadrature.integrate(integrand, 0.0, 1.0, tolerance, tolerance, 0);

		assertAbsRelEqual(1.0, result[0], tolerance);
		assertAbsRelEqual(0.0, result[1], tolerance);
	}

	@Test
	public void testExpX2() {

		final double tolerance = 1.0e-12;

		class ExpX2 implements Integrand {
			@Override
			public double[] eval(double... x) {
				int n = x.length;
				double[] f = new double[n];
				for (int i=0; i<n; ++i) {
					f[i] = Math.exp(-x[i]*x[i]);
				}
				return f;
			}
		};

		ExpX2 integrand = new ExpX2();

		double[] result = AdaptiveQuadrature.integrate(integrand, 0.0, Double.POSITIVE_INFINITY, tolerance, tolerance, 0);

		assertAbsRelEqual(Math.sqrt(Math.PI)/2.0, result[0], tolerance);
		assertAbsRelEqual(0.0, result[1], tolerance);
	}

	/**
	 * Combined check on relative and absolute deviation between expected and actual values.
	 * From Gill, Murray & Wright, "Practial Optimization" (1984).
	 *
	 * @param expected
	 * @param actual
	 * @param relAbsTolerance
	 */
	private static void assertAbsRelEqual(double expected, double actual, double relAbsTolerance) {
		final double relAbsError = Math.abs(actual - expected) / (1.0 + Math.abs(expected));
		if (relAbsError > relAbsTolerance) {
			fail(String.format("expected %g, actual %g --> rel/abs error %g (tol %g)",
					expected, actual, relAbsError, relAbsTolerance));
		}
	}

}
