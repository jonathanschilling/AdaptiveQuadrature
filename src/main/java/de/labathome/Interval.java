package de.labathome;

public class Interval implements Comparable<Interval> {

	private double integralValue;
	private double errorEstimate;

	private double center;
	private double halfWidth;

	public Interval(double center, double halfWidth) {
		this.integralValue = 0.0;
		this.errorEstimate = Double.POSITIVE_INFINITY;

		this.center = center;
		this.halfWidth = halfWidth;
	}

	public Interval cutInHalf() {
		this.halfWidth /= 2.0;
		Interval otherHalf = new Interval(this.center + this.halfWidth, this.halfWidth);
		this.center -= this.halfWidth;
		return otherHalf;
	}

	public void setIntegralValue(double integralValue) {
		this.integralValue = integralValue;
	}

	public double getIntegralValue() {
		return integralValue;
	}

	public double getErrorEstimate() {
		return errorEstimate;
	}

	public double getCenter() {
		return center;
	}

	public double getHalfWidth() {
		return halfWidth;
	}


	@Override
	public int compareTo(Interval other) {
		if (this.errorEstimate == other.errorEstimate) return 0;

		/**
		 * Java implements the PriorityQueue type as a MinHeap, which means that the root
		 * is always the smallest element as judged by the compareTo method of the defining
		 * type. The integration algorithm however expects a MaxHeap (where the root is the
		 * largest element) in order to work on the worst regions first.
		 * Hence, we need to reverse the order of the Java PriorityQueue by reversing the sign
		 * of the compareTo method output.
		 */
		return (this.errorEstimate < other.errorEstimate) ? 1 : -1;
	}

}
