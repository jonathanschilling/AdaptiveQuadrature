package de.labathome.tests;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.UnaryOperator;

import de.labathome.AdaptiveQuadrature;

public class TestAdaptiveQuadrature {

	/**
	 * Test a definite integral:
	 * <pre>
	 * 1/3 \int_0^1 x^2 dx = 1
	 * </pre>
	 */
	@Test
	public void testParabola() {

		final double tolerance = 1.0e-12;

		UnaryOperator<double[]> parabola = x -> {
			int n = x.length;
			double[] f = new double[n];
			for (int i=0; i<n; ++i) {
				f[i] = 3.0*x[i]*x[i];
			}
			return f;
		};

		double[] result = AdaptiveQuadrature.integrate(parabola, 0.0, 1.0, tolerance, tolerance, 0);

		assertAbsRelEqual(1.0, result[0], tolerance);
		assertAbsRelEqual(0.0, result[1], tolerance);
	}

	/**
	 * Test a semi-indefinite integral:
	 * <pre>
	 * \int_0^{+\infty} exp(-x^2} dx = \sqrt(pi)/2
	 * </pre>
	 */
	@Test
	public void testExpX2() {

		final double tolerance = 1.0e-12;

		UnaryOperator<double[]> expX2 = x -> {
			int n = x.length;
			double[] f = new double[n];
			for (int i=0; i<n; ++i) {
				f[i] = Math.exp(-x[i]*x[i]);
			}
			return f;
		};

		double[] result = AdaptiveQuadrature.integrate(expX2, 0.0, Double.POSITIVE_INFINITY, tolerance, tolerance, 0);

		assertAbsRelEqual(Math.sqrt(Math.PI)/2.0, result[0], tolerance);
		assertAbsRelEqual(0.0, result[1], tolerance);
	}

	/**
	 * Test that AdaptiveQuadrature can be gracefully stopped when a certain condition is fulfilled within the integrand.
	 */
	@Test
	public void testGracefulStop() {

		final double tolerance = 1.0e-12;

		final int maxEval = 100;

		class ExpX2 implements UnaryOperator<double[]> {

			private AtomicBoolean gracefulStop;
			private int numEval;

			public ExpX2(AtomicBoolean gracefulStop) {
				this.gracefulStop = gracefulStop;
				this.numEval = 0;
			}

			@Override
			public double[] apply(double[] x) {
				int n = x.length;
				numEval += n;

				if (numEval > maxEval) {
//					System.out.printf("evaluated at %3d points --> stop\n", numEval);
					while (!gracefulStop.compareAndSet(false, true)) {
						// wait for gracefulStop to turn false again
					}
					return null;
				}
//				else {
//					System.out.printf("evaluated at %3d points --> go on\n", numEval);
//				}

				double[] f = new double[n];
				for (int i=0; i<n; ++i) {
					f[i] = Math.exp(-x[i]*x[i]);
				}
				return f;
			}

			public int getNumEval() {
				return numEval;
			}
		};

		AtomicBoolean gracefulStop = new AtomicBoolean(false);
		ExpX2 expX2 = new ExpX2(gracefulStop);

		AdaptiveQuadrature.integrate(expX2, 0.0, Double.POSITIVE_INFINITY, tolerance, tolerance, 0, gracefulStop);

		int numEval = expX2.getNumEval();

		assertEquals(true, gracefulStop.get());
		assertTrue(numEval > maxEval);
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
