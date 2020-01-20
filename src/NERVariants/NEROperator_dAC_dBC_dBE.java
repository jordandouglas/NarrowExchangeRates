package NERVariants;

import operators.MetaNEROperator;
import beast.util.Randomizer;

public class NEROperator_dAC_dBC_dBE extends MetaNEROperator {

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
		this.rap = (2*te*(r_delta + rd) + ra*ta - ra*td + rd*td - rd*te - 2*(r_delta + rd)*(r_epsilon + td))/(ta - te);
		this.rbp = (r_delta*r_epsilon + r_epsilon*rd + r_delta*td - rb*tb - r_delta*te + rb*td)/(r_epsilon - tb + td);
		this.rcp = (te*(r_delta + rd) - rc*tc + rc*te - (r_delta + rd)*(r_epsilon + td))/(r_epsilon - tc + td);
		this.rdp = r_delta + rd;
		this.tdp = r_epsilon + td;


		// Jacobian determinant
		double JD = ((ta - td)*(tb - td)*(tc - te))/((ta - te)*(r_epsilon - tb + td)*(r_epsilon - tc + td));
		if (JD <= 0) return Double.NEGATIVE_INFINITY;
		return Math.log(JD);

	}

}
