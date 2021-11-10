package de.labathome;

import java.util.List;
import java.util.function.UnaryOperator;

/**
 * Poles and Weights for Gauss-Kronrod quadrature.
 *
 * The numerical values have been obtained using the method by Laurie [1]
 * using the implementation presented by Holoborodko [2].
 *
 * @see [1] https://doi.org/10.1090/S0025-5718-97-00861-2
 * @see [2] https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
 */
public class GaussKronrod {

	/** ca. 2.22e-16 for 64-bit double precision */
	private static final double DBL_EPSILON = Math.ulp(1.0);

	public static final int NUM_GAUSS_POINTS = 7;
	public static final int NUM_KRONROD_POINTS = 2*NUM_GAUSS_POINTS+1;

	// Nodes marked by (*) belong to the embedded Gauss-Legendre quadrature.

	/** poles for 7/15-point Gauss-Kronrod quadrature */
	public static double[] P_GK_7_15  = {
		// 0.0                                          0 (*)
		2.077849550078984676006894037732449e-01, // +/- 1
		4.058451513773971669066064120769615e-01, // +/- 2 (*)
		5.860872354676911302941448382587296e-01, // +/- 3
		7.415311855993944398638647732807884e-01, // +/- 4 (*)
		8.648644233597690727897127886409262e-01, // +/- 5
		9.491079123427585245261896840478513e-01, // +/- 6 (*)
		9.914553711208126392068546975263285e-01  // +/- 7
	};

	/** weights for 7-point Gauss-Legendre quadrature */
	public static double[] W_GAUSS = {
		4.179591836734693877551020408163265e-01, //     0 (*)
		3.818300505051189449503697754889751e-01, // +/- 1 (*)
		2.797053914892766679014677714237796e-01, // +/- 2 (*)
		1.294849661688696932706114326790820e-01  // +/- 3 (*)
	};

	/** weights for 15-point Kronrod quadrature */
	public static double[] W_KRNRD = {
		2.094821410847278280129991748917143e-01, //     0 (*)
		2.044329400752988924141619992346491e-01, // +/- 1
		1.903505780647854099132564024210137e-01, // +/- 2 (*)
		1.690047266392679028265834265985503e-01, // +/- 3
		1.406532597155259187451895905102379e-01, // +/- 4 (*)
		1.047900103222501838398763225415180e-01, // +/- 5
		6.309209262997855329070066318920429e-02, // +/- 6 (*)
		2.293532201052922496373200805896959e-02  // +/- 7
	};

	public static void evalGaussKronrod(UnaryOperator<double[]> integrand, List<Interval> intervals) {
		int numIntervals = intervals.size();

		double[] evaluationLocations = new double[numIntervals*NUM_KRONROD_POINTS];

		int intervalStartIdx;
		double center, halfWidth;
		for (int i=0; i<numIntervals; ++i) {

			Interval interval = intervals.get(i);
			center = interval.getCenter();
			halfWidth = interval.getHalfWidth();

			intervalStartIdx = i*NUM_KRONROD_POINTS;

			evaluationLocations[intervalStartIdx] = center;
//			System.out.printf(" center: f(   %d ) into %2d\n", 0, 0);

			// 7-point Gauss-Legendre
			for (int j=1; j < NUM_GAUSS_POINTS; j += 2) {
//				System.out.printf("  Gauss: f(-p[%d]) into %2d ; f(+p[%d]) into %2d\n", j, j, j, j+1);
				evaluationLocations[intervalStartIdx+j  ] = center-halfWidth*GaussKronrod.P_GK_7_15[j];
				evaluationLocations[intervalStartIdx+j+1] = center+halfWidth*GaussKronrod.P_GK_7_15[j];
			}

			// additional 8 points for 15-point Kronrod
			for (int j=0; j<NUM_GAUSS_POINTS; j+=2) {
//				System.out.printf("Kronrod: f(-p[%d]) into %2d ; f(+p[%d]) into %2d\n", j, NUM_GAUSS_POINTS+j, j, NUM_GAUSS_POINTS+j+1);
				evaluationLocations[intervalStartIdx+NUM_GAUSS_POINTS+j  ] = center-halfWidth*GaussKronrod.P_GK_7_15[j];
				evaluationLocations[intervalStartIdx+NUM_GAUSS_POINTS+j+1] = center+halfWidth*GaussKronrod.P_GK_7_15[j];
			}
		}

		// evaluate integrand
		final double[] functionValues = integrand.apply(evaluationLocations);

		// compute Gauss-Kronrod quadrature results

		double gVal, kVal1, kVal2;
		double resultGauss, resultKrnrd;
		double resultAbs, resultResidual;
		double errorEstimate, mean, scale;
		for (int i=0; i<numIntervals; ++i) {
			Interval interval = intervals.get(i);
			halfWidth = interval.getHalfWidth();

			intervalStartIdx = i*NUM_KRONROD_POINTS;

			resultGauss = functionValues[intervalStartIdx] * W_GAUSS[0];
			resultKrnrd = functionValues[intervalStartIdx] * W_KRNRD[0];
			resultAbs = Math.abs(resultKrnrd);

			// 7-point Gauss-Legendre
			for (int j=1; j <= NUM_GAUSS_POINTS/2; ++j) {
//				System.out.printf("  Gauss: f[%2d]*wG[%d] ; f[%2d]*wG[%d]\n", 2*j-1, j, 2*j, j);
//				System.out.printf("Kronrod: f[%2d]*wK[%d] ; f[%2d]*wK[%d]\n", 2*j-1, 2*j, 2*j, 2*j);
				gVal = (functionValues[intervalStartIdx+2*j-1] + functionValues[intervalStartIdx+2*j])*W_GAUSS[j];
				kVal1 = functionValues[intervalStartIdx+2*j-1]*W_KRNRD[2*j];
				kVal2 = functionValues[intervalStartIdx+2*j  ]*W_KRNRD[2*j];
				resultGauss += gVal;
				resultKrnrd += kVal1 + kVal2;
				resultAbs += Math.abs(kVal1) + Math.abs(kVal2);
			}

			// additional 8 points for 15-point Kronrod
			for (int j=0; j <= NUM_GAUSS_POINTS/2; ++j) {
//				System.out.printf("Kronrod: f[%2d]*wK[%d] ; f[%2d]*wK[%d]\n", NUM_GAUSS_POINTS+2*j, 2*j+1, NUM_GAUSS_POINTS+2*j+1, 2*j+1);
				kVal1 = functionValues[intervalStartIdx+NUM_GAUSS_POINTS+2*j  ]*W_KRNRD[2*j+1];
				kVal2 = functionValues[intervalStartIdx+NUM_GAUSS_POINTS+2*j+1]*W_KRNRD[2*j+1];
				resultKrnrd += kVal1 + kVal2;
				resultAbs += Math.abs(kVal1) + Math.abs(kVal2);
			}
			resultAbs *= halfWidth;
			interval.setIntegralValue(resultKrnrd*halfWidth);

			errorEstimate = Math.abs(resultKrnrd - resultGauss) * halfWidth;

			// subtract half the actual result from sum of absolute contributions
			mean = resultKrnrd * 0.5;
			resultResidual = Math.abs(functionValues[intervalStartIdx] - mean) * W_KRNRD[0];
			for (int j=1; j <= NUM_GAUSS_POINTS/2; ++j) {
				kVal1 = Math.abs(functionValues[intervalStartIdx+2*j-1] - mean) * W_KRNRD[2*j];
				kVal2 = Math.abs(functionValues[intervalStartIdx+2*j  ] - mean) * W_KRNRD[2*j];
				resultResidual += kVal1 + kVal2;
			}
			for (int j=0; j <= NUM_GAUSS_POINTS/2; ++j) {
				kVal1 = Math.abs(functionValues[intervalStartIdx+NUM_GAUSS_POINTS+2*j  ] - mean) * W_KRNRD[2*j+1];
				kVal2 = Math.abs(functionValues[intervalStartIdx+NUM_GAUSS_POINTS+2*j+1] - mean) * W_KRNRD[2*j+1];
				resultResidual += kVal1 + kVal2;
			}
			resultResidual *= halfWidth;

			if (resultResidual != 0.0 && errorEstimate != 0.0) {
				scale = Math.pow(200.0 * errorEstimate / resultResidual, 1.5);
				errorEstimate = (scale < 1.0) ? resultResidual * scale : resultResidual;
			}

			// minimum possible error due to accumulation of roundoff-errors
			if (resultAbs > Double.MIN_NORMAL / (50 * DBL_EPSILON)) {
				double min_err = 50 * DBL_EPSILON * resultAbs;
				if (min_err > errorEstimate) {
					errorEstimate = min_err;
				}
			}

			interval.setErrorEstimate(errorEstimate);
		}
	}
}
