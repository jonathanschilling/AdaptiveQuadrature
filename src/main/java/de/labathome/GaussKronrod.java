package de.labathome;

import java.util.List;

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
	public static double[] W_GAUSS_7 = {
		4.179591836734693877551020408163265e-01, //     0
		3.818300505051189449503697754889751e-01, // +/- 1
		2.797053914892766679014677714237796e-01, // +/- 2
		1.294849661688696932706114326790820e-01  // +/- 3
	};

	/** weights for 15-point Kronrod quadrature */
	public static double[] W_KRONROD_15 = {
		2.094821410847278280129991748917143e-01, //     0 (*)
		2.044329400752988924141619992346491e-01, // +/- 1
		1.903505780647854099132564024210137e-01, // +/- 2 (*)
		1.690047266392679028265834265985503e-01, // +/- 3
		1.406532597155259187451895905102379e-01, // +/- 4 (*)
		1.047900103222501838398763225415180e-01, // +/- 5
		6.309209262997855329070066318920429e-02, // +/- 6 (*)
		2.293532201052922496373200805896959e-02  // +/- 7
	};

	public static double[][] evalGaussKronrod(Integrand integrand, List<Interval> intervals) {

		int numIntervals = intervals.size();

		final int n = 7;
		double[] evaluationLocations = new double[numIntervals*(2*n+1)];

		for (int i=0; i<numIntervals; ++i) {
			final Interval interval = intervals.get(i);
			final double center = interval.getCenter();
			final double halfWidth = interval.getHalfWidth();

			evaluationLocations[0] = center;

			// 7-point Gauss-Legendre
			for (int j=1; j < n; j += 2) {
				System.out.println("Gauss: f(-p["+j+"]) into "+(j)+" ; f(+p["+j+"]) into "+(j+1));
				evaluationLocations[j  ] = center-halfWidth*GaussKronrod.P_GK_7_15[j];
				evaluationLocations[j+1] = center+halfWidth*GaussKronrod.P_GK_7_15[j];
			}

			// additional 8 points for 15-point Kronrod
			for (int j=0; j<n; j+=2) {
				System.out.println("Kronrod: f(-p["+j+"]) into "+(n+j)+" ; f(+p["+j+"]) into "+(n+j+1));
				evaluationLocations[n+j  ] = center-halfWidth*GaussKronrod.P_GK_7_15[j];
				evaluationLocations[n+j+1] = center+halfWidth*GaussKronrod.P_GK_7_15[j];
			}
		}

		// evaluate integrand
		final double[] functionValues = integrand.eval(evaluationLocations);


		// compute Gauss-Kronrod quadrature result
		final double[] gaussResults = new double[numIntervals];
		final double[] kronrodResults = new double[numIntervals];

		for (int i=0; i<numIntervals; ++i) {




		}

		return new double[][] {gaussResults, kronrodResults};
	}



}
