package de.labathome;

import java.util.PriorityQueue;

public class IntervalHeap extends PriorityQueue<Interval> {

	/** auto-generated serial version UID */
	private static final long serialVersionUID = 3211016244181536292L;

	private double integralValue;
	private double errorEstimate;

	public IntervalHeap() {
		super();
		this.integralValue = 0.0;
		this.errorEstimate = 0.0;
	}

	public double getIntegralValue() {
		return integralValue;
	}

	public double getErrorEstimate() {
		return errorEstimate;
	}

	public Interval poll() {
		Interval ret = super.poll();
		integralValue -= ret.getIntegralValue();
		errorEstimate -= ret.getErrorEstimate();
		return ret;
	}

	public boolean add(Interval e) {
		this.integralValue += e.getIntegralValue();
		this.errorEstimate += e.getErrorEstimate();
		return super.add(e);
	}

}
