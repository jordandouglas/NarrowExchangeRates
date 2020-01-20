package NERVariants;

import operators.MetaNEROperator;
import beast.util.Randomizer;

public class NEROperator_dCE extends MetaNEROperator {

	@Override
	protected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te, 
										double ra, double rb, double rc, double rd) {

		// Random proposals
		double r_alpha = 0;
		double r_beta = 0;
		double r_gamma = 0;
		double r_delta = 0;
		double r_epsilon = (Randomizer.nextFloat()-0.5) *2*rWindowSize;


		// Propose new rates + times
		this.rap = r_alpha + ra;
		this.rbp = r_beta + rb;
		this.rcp = -(te*(r_delta + rd) + rc*tc - rc*te - (r_delta + rd)*(r_epsilon + td))/(r_epsilon - tc + td);
		this.rdp = r_delta + rd;
		this.tdp = r_epsilon + td;


		// Jacobian determinant
		double JD = -(tc - te)/(r_epsilon - tc + td);
		if (JD <= 0) return Double.NEGATIVE_INFINITY;
		return Math.log(JD);

	}

}