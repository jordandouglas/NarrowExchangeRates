package NERVariants;

import operators.MetaNEROperator;

public class NEROperator_ABCDE extends MetaNEROperator {


	
	@Override
	protected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te, 
												  double ra, double rb, double rc, double rd) {
		
		// Propose new rates + times
		rap = (ra*ta - ra*td + rd*td - rd*te)/(ta - te);
		rbp = -(te*(rd) + rb*tb - rb*td - rd*td + rd*te - (rd)*(td))/(-tb + td);
		rcp = -(te*(rd) + rc*tc - rc*te - (rd)*(td))/(-tc + td);
		rdp = rd;
		tdp = td;
		
		// Calculate logJD
		double JD = ((ta - td)*(tb - td)*(tc - te))/((ta - te)*(- tb + td)*(- tc + td));
		if (JD <= 0) return Double.NEGATIVE_INFINITY;
		return Math.log(JD);
		
	}
	

}
