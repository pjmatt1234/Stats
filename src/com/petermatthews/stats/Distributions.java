package com.petermatthews.stats;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;

public class Distributions {

	public static double normPdf(double mu, double sigma, double x) {
		double exp = -((Math.pow(x - mu, 2)) / (2 * Math.pow(sigma, 2)));
		return (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.pow(Math.E, exp);
	}

	public static double normCdf(double mu, double sigma, double x) {
		double z = (x - mu) / sigma;
		if (z < -8.0)
			return 0.0;
		if (z > 8.0)
			return 1.0;
		double sum = 0.0, term = z;
		for (int i = 3; sum + term != sum; i += 2) {
			sum = sum + term;
			term = term * z * z / i;
		}
		return 0.5 + sum * Distributions.normPdf(0.0, 1.0, z);
	}

	public static double invNorm(double area, double mu, double sigma) {
		NormalDistribution normDist = new NormalDistribution(mu, sigma);

		return normDist.inverseCumulativeProbability(area); // to be implemented
	}
	
	public static double getdf2Samp(double sigma1, double sigma2, double n1, double n2) {
		double s1 = sigma1 * sigma1;
		double s2 = sigma2 * sigma2;
		
		double num = Math.pow( ((s1)/n1) + ((s2)/n2) , 2);
		double denom = ((Math.pow((s1/n1),2) / (n1 -1))) + ((Math.pow((s2/n2),2) / (n2 -1)));
		
		
		return num/denom;
		
	}
	
	public static double tPdf(double x, double df) {
		TDistribution tDist = new TDistribution(df);
		return tDist.density(x);
	}

	public static double tCdf(double x, double df) {
		TDistribution tDist = new TDistribution(df);
		return tDist.cumulativeProbability(x);
	}

	public static double invT(double area, double df) {
		TDistribution tDist = new TDistribution(df);
		return tDist.inverseCumulativeProbability(area);
	}
	
	public static double chisqPdf(double x, double df) {
		ChiSquaredDistribution xsqDist = new ChiSquaredDistribution(df);
		return xsqDist.density(x);
	}
	
	public static double chisqCdf(double x, double df) {
		ChiSquaredDistribution xsqDist = new ChiSquaredDistribution(df);
		return xsqDist.cumulativeProbability(x);
	}
	
	public static double invChiSq(double area, double df) {
		ChiSquaredDistribution xsqDist = new ChiSquaredDistribution(df);
		return xsqDist.inverseCumulativeProbability(area);
	}
	
	public static double binomPdf(double n, double p, double x) {
		BinomialDistribution binomPdf = new BinomialDistribution((int)n, p);
		return binomPdf.probability((int)x);
	}
	
	public static double geomPdf(double x, double p){
		return Math.pow(1-p, x-1) * p;
	}

	public static double getRowExp(Double[][] arr, int index){
		double total = 0.0;

		for(int i = 0;i<arr[0].length;i++){
		    total += arr[index][i];
        }
        return total;
	}

	public static double getColExp(Double[][] arr, int index){
		double total = 0.0;
		for(int i = 0; i<arr[0].length;i++){
		    total += arr[i][index];
        }
        return total;
	}


	
}