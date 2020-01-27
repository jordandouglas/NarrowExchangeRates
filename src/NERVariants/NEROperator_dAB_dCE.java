package NERVariants;

import operators.MetaNEROperator;
import beast.util.Randomizer;

public class NEROperator_dAB_dCE extends MetaNEROperator {

	@Override
	protected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te, 
										double ra, double rb, double rc, double rd) {

		// Random proposals
		double r_alpha = 0;
		double r_beta = 0;
		double r_gamma = 0;
		double r_delta = 0;
		double r_epsilon = this.getRandomWalkStepSize(rWindowSize);


		// Propose new rates + times
		this.rap = -(tb*(r_beta + rb) - te*(r_delta + rd) - ra*ta - rb*tb + ra*td + rb*td - (r_beta + rb)*(r_epsilon + td) + (r_delta + rd)*(r_epsilon + td))/(ta - te);
		this.rbp = r_beta + rb;
		this.rcp = -(te*(r_delta + rd) + rc*tc - rc*te - (r_delta + rd)*(r_epsilon + td))/(r_epsilon - tc + td);
		this.rdp = r_delta + rd;
		this.tdp = r_epsilon + td;


		// Jacobian determinant
		double JD = -((ta - td)*(tc - te))/((ta - te)*(r_epsilon - tc + td));
		if (JD <= 0) return Double.NEGATIVE_INFINITY;
		return Math.log(JD);

	}

}
