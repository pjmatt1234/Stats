package com.petermatthews.stats;
import java.util.Scanner;

class Main {
    public static void main(String[] args) {
		while(true) {
			System.out.println("--------AP Stats Software--------");
			System.out.println("Please enter what you would like to do");
			System.out.println("1. Distributions\n" + "2. Confidence Intervals\n" + "3. Significance Tests");
			Scanner scan = new Scanner(System.in);
			int response = scan.nextInt();
			int row, col;
			double mu, sigma, x, x1, x2, area, df, n, n1, n2, p, p1, ans, lo, hi, cL, z, t, p2, sigma2, sigma3, chisq, total;
			switch (response) {
			case 1:
				// open up Distributions
				System.out.println("--------Distribution Menu--------\n" + "Please enter what you would like to do\n"
						+ "1. Normal Pdf\n" + "2. Normal Cdf\n" + "3. Inverse Normal\n" + "4. t Pdf\n" + "5. t Cdf\n"
						+ "6. Inverse t\n" + "7. X^2 Pdf\n" + "8. X^2 Cdf\n" + "9. Inverse X^2\n" + "10. Binomial Pdf\n"
						+ "11. Binomial Cdf\n" + "12. Geometric Pdf\n" + "13. Geometric Cdf\n");
				int distResponse = scan.nextInt();

				switch (distResponse) {
				case 1:
					System.out.println("Please enter mu");
					mu = scan.nextDouble();
					System.out.println("Please enter sigma as");
					sigma = scan.nextDouble();
					System.out.println("Please enter the expected X value");
					x = scan.nextDouble();
					System.out.println("Norm Pdf(mu = " + mu + ", sigma = " + sigma + ", x = " + x + ") = "
							+ Distributions.normPdf(mu, sigma, x));
					break;
				case 2:
					System.out.println("Please enter mu");
					mu = scan.nextDouble();
					System.out.println("Please enter sigma");
					sigma = scan.nextDouble();
					System.out.println("Please enter the expected x value");
					x = scan.nextDouble();

					System.out.println(
							"Norm Cdf(mu = " + mu + ", sigma = " + sigma + ", lower bound = -infinity, upper bound = "
									+ x + ") = " + Distributions.normCdf(mu, sigma, x));
					break;
				case 3:
					System.out.println("Please enter the area");
					area = scan.nextDouble();
					System.out.println("Please enter mu");
					mu = scan.nextDouble();
					System.out.println("Please enter sigma");
					sigma = scan.nextDouble();

					System.out.println("InvNorm(area = " + area + ", mu = " + mu + ", sigma = " + sigma + ") = "
							+ Distributions.invNorm(area, mu, sigma));

					break;
				case 4:
					System.out.println("Please enter your X value as a double");
					x = scan.nextDouble();
					System.out.println("Please enter your degrees of freedom");
					df = scan.nextDouble();

					System.out.println("tPdf(x= " + x + ", df= " + df + ") = " + Distributions.tPdf(x, df));
					break;
				case 5:
					System.out.println("Please enter your X value");
					x = scan.nextDouble();
					System.out.println("Please enter your degrees of freedom");
					df = scan.nextDouble();

					System.out.println("tCdf(df= " + df + ", lowerbound = -infinity, upperbound = " + x + ") = "
							+ Distributions.tCdf(x, df));
					break;
				case 6:
					System.out.println("Please enter your area");
					area = scan.nextDouble();
					System.out.println("Please enter your degrees of freedom");
					df = scan.nextDouble();

					System.out.println("invT(area = " + area + ", df = " + df + ") = " + Distributions.invT(area, df));
					break;
				case 7:
					System.out.println("Please enter your x value");
					x = scan.nextDouble();
					System.out.println("Please enter your degrees of freedom");
					df = scan.nextDouble();

					System.out.println("x^2Pdf(x= " + x + ", df= " + df + ") = " + Distributions.chisqPdf(x, df));
					break;
				case 8:
					System.out.println("Please enter your x value");
					x = scan.nextDouble();
					System.out.println("Please enter your degrees of freedom");
					df = scan.nextDouble();

					System.out.println("x^2Cdf(df= " + df + ", lowerbound = 0, upperbound = " + x + ") = "
							+ Distributions.chisqCdf(x, df));
					break;
				case 9:
					System.out.println("Please enter your area");
					area = scan.nextDouble();
					System.out.println("Please enter your degrees of freedom");
					df = scan.nextDouble();

					System.out.println(
							"invX^2(area = " + area + ", df = " + df + ") = " + Distributions.invChiSq(area, df));

					break;
				case 10:
					System.out.println("Please enter the number of trials");
					n = scan.nextDouble();
					System.out.println("Please enter the probability of success");
					p = scan.nextDouble();
					System.out.println("Please enter your X value");
					x = scan.nextDouble();

					System.out.println("binomPdf(n = " + n + ", prob = " + p + ", x = " + x + ") = "
							+ Distributions.binomPdf(n, p, x));

					break;
				case 11:
					System.out.println("Please enter the number of trials");
					n = scan.nextDouble();
					System.out.println("Please enter the probability of success");
					p = scan.nextDouble();
					System.out.println("Please enter your X value");
					x = scan.nextInt();
					ans = 0.0;
					for (int i = (int) x; i >= 0; i--) {
						ans += Distributions.binomPdf(n, p, i);
					}

					System.out.println("binomCdf(n = " + n + ", prob = " + p + ", x = " + x + ") = " + ans);

					break;
				case 12:
					System.out.println("Please enter the probability of success");
					p = scan.nextDouble();
					System.out.println("Please enter your X value");
					x = scan.nextDouble();

					System.out.println("geomPdf(p = " + p + ", x = " + x + ") = " + Distributions.geomPdf(x, p));
					break;
				case 13:
					System.out.println("Please enter the probability of success");
					p = scan.nextDouble();
					System.out.println("Please enter the low bound");
					lo = scan.nextDouble();
					System.out.println("Please enter the high bound");
					hi = scan.nextDouble();
					ans = 0.0;
					int low = (int) lo;

					for (int high = (int) hi; high >= low; high--) {
						ans += Distributions.geomPdf((double) high, p);
					}
					System.out.println(
							"geomCdf(p = " + p + ", low bound = " + lo + ", high bound = " + hi + ") = " + ans);
					break;
				}
				break;
			case 2:
				// open up Confidence Intervals
				System.out.println("--------Confidence Intervals--------\n" + "Please enter what you would like to do\n"
						+ "1. z Interval\n" + "2. t Interval\n" + "3. 1-Proportion z Interval\n"
						+ "4. 2-Proportion z Interval\n" + "5. 2 Sample t Interval");
				int ciResponse = scan.nextInt();
				switch (ciResponse) {
				case 1:

					System.out.println("Please enter sigma");
					sigma = scan.nextDouble();
					System.out.println("Please enter x-bar");
					mu = scan.nextDouble();
					System.out.println("Please enter n");
					n = scan.nextDouble();
					System.out.println("Please enter your confidence level (0-1)");
					cL = scan.nextDouble();
					z = Distributions.invNorm((1 - cL) / 2, 0, 1);
					hi = mu - z * (sigma / Math.sqrt(n));
					lo = mu + z * (sigma / Math.sqrt(n));

					System.out.println("zInterval(Lower: " + lo + ", Higher = " + hi + ", x-bar = " + mu + ", sigma = "
							+ sigma + ", n = " + n + ", ME = " + -(z * (sigma / Math.sqrt(n))) + ")");
					break;
				case 2:
					System.out.println("Please enter sigma");
					sigma = scan.nextDouble();
					System.out.println("Please enter x-bar");
					mu = scan.nextDouble();
					System.out.println("Please enter n");
					n = scan.nextDouble();
					System.out.println("Please enter your confidence level (0-1)");
					cL = scan.nextDouble();
					t = Distributions.invT((1 - cL) / 2, n - 1);

					hi = mu - t * (sigma / Math.sqrt(n));
					lo = mu + t * (sigma / Math.sqrt(n));

					System.out.println("tInterval(Lower: " + lo + ", Higher = " + hi + ", x-bar = " + mu + ", sigma = "
							+ sigma + ", n = " + n + ", ME = " + -(t * (sigma / Math.sqrt(n))) + ")");

					break;
				case 3:
					System.out.println("Please enter n");
					n = scan.nextDouble();
					System.out.println("Please enter the number of successes");
					mu = scan.nextDouble() / n;
					System.out.println("Please enter your confidence level (0-1)");
					cL = scan.nextDouble();
					z = Distributions.invNorm((1 - cL) / 2, 0, 1);

					hi = mu - z * (Math.sqrt(((mu) * (1 - mu)) / n));
					lo = mu + z * (Math.sqrt(((mu) * (1 - mu)) / n));

					System.out.println("1-Prop zInterval(Lower: " + lo + ", Higher = " + hi + ", p-hat = " + mu
							+ ", n = " + n + ", ME = " + -(z * (Math.sqrt(((mu) * (1 - mu)) / n))) + ")");

					break;
				case 4:
					System.out.println("Please enter n1");
					n = scan.nextDouble();
					System.out.println("Please enter the number successes for n1");
					p = (scan.nextDouble() / n);
					System.out.println("Please enter n2");
					n2 = scan.nextDouble();
					System.out.println("Please enter the number successes for n2");
					p2 = (scan.nextDouble() / n2);
					System.out.println("Please enter your confidence level (0-1)");
					cL = scan.nextDouble();
					z = Distributions.invNorm((1 - cL) / 2, 0, 1);

					sigma = Math.sqrt((p * (1 - p)) / n);
					sigma2 = Math.sqrt((p2 * (1 - p2)) / n2);
					sigma3 = Math.sqrt(sigma * sigma + sigma2 * sigma2);

					hi = (p - p2) - z * sigma3;
					lo = (p - p2) + z * sigma3;

					System.out.println("2-Prop zInterval(Lower: " + lo + ", Higher = " + hi + ", p-hat1 = " + p
							+ ", n1 = " + n + ", p-hat2 = " + p2 + ", n2 = " + n2 + ", ME = " + -(z * sigma3) + ")");
					break;
				case 5:
					System.out.println("Please enter x-bar one");
					x1 = scan.nextDouble();
					System.out.println("Please enter sigma one");
					sigma = scan.nextDouble();
					System.out.println("Please enter n1");
					n1 = scan.nextDouble();
					System.out.println("Please enter x-bar two");
					x2 = scan.nextDouble();
					System.out.println("Please enter sigma two");
					sigma2 = scan.nextDouble();
					System.out.println("Please enter n2");
					n2 = scan.nextDouble();
					System.out.println("Please enter the confidence level (0-1)");
					cL = scan.nextDouble();
					x = x1 - x2;

					sigma3 = Math.sqrt(((sigma * sigma) + (sigma2 * sigma2)) / (n1 + n2 - 2));


					z = Distributions.invT((1 - cL) / 2, Distributions.getdf2Samp(sigma, sigma2, n1, n2));
					hi = x - z * sigma3;
                    lo = x + z * sigma3;

					System.out.println(
							"2-Sample tInterval(Lower: " + lo + ", Higher = " + hi + ", x-bar1 = " + x1 + ", sigma1 = "
									+ sigma + ", sigma2 = " + sigma2 + ", n1 = " + n1 + ", x-bar2 = " + x2 + ", n2 = "
									+ n2 + ", x-bar = " + x + ", sigma = " + sigma3 + ", ME = " + -(z * sigma3) + ")");

					break;
				}
				break;
			case 3:
				// open up Significance Tests
				System.out.println("--------Significance Tests--------\n" + "1. z Test\n" + "2. t Test\n"
						+ "3. 1-Proportion z Test\n" + "4. 2-Proportion z Test\n" + "5. 2-Sample z Test\n"
						+ "6. 2-Sample t Test\n" + "7. X^2 Goodness of Fit\n" + "8. X^2 2-Way Test");

				int sigResponse = scan.nextInt();
				switch (sigResponse) {
				case 1:
					System.out.println("Please enter the null hypothesis' proportion");
					p = scan.nextDouble();
					System.out.println("Please enter the observed standard deviation");
					sigma = scan.nextDouble();
					System.out.println("Please enter the observed proportion");
					x = scan.nextDouble();
					System.out.println("Please enter the sample size");
					n = scan.nextDouble();
					System.out.println("Please enter a numeric response for the Alternative Hypothesis"
							+ "\n1. p-hat != p" + "\n2. p-hat > p" + "\n3. p-hat < p");
					ans = scan.nextDouble();
					switch ((int)ans) {
					case 1:
						area = 2 * Distributions.normCdf(p, sigma/Math.sqrt(n), x);
						System.out.println("zTest(alt hypo: p-hat != p, z = "+((x-p)/(sigma/Math.sqrt(n)))+", PVal = "+area+")");
						break;
					case 2:
						area = 1 - Distributions.normCdf(p, sigma/Math.sqrt(n), x);
						System.out.println("zTest(alt hypo: p-hat > p, z = "+((x-p)/(sigma/Math.sqrt(n)))+", PVal = "+area+")");
						break;
					case 3:
						area = Distributions.normCdf(p, sigma/Math.sqrt(n), x);
						System.out.println("zTest(alt hypo: p-hat < p, z = "+((x-p)/(sigma/Math.sqrt(n)))+", PVal = "+area+")");
						break;
					}
					break;
				case 2:
					System.out.println("Please enter the null hypothesis' mu");
					mu = scan.nextDouble();
					System.out.println("Please enter the x-bar");
					x = scan.nextDouble();
					System.out.println("Please enter the observed standard deviation");
					sigma = scan.nextDouble();
					System.out.println("Please enter the sample size");
					n = scan.nextDouble();
					
					System.out.println("Please enter a numeric response for the Alternative Hypothesis"
							+ "\n1. x-bar != mu" + "\n2. x-bar > mu" + "\n3. x-bar < mu");
					ans = scan.nextDouble();
					t = (x-mu) / ( sigma / Math.sqrt(n));
					switch ((int)ans) {
					case 1:
						area = 2 * Distributions.tCdf(-Math.abs(t), n-1);
						System.out.println("tTest(alt hypo: mu != x-bar, t = "+((x-mu)/(sigma/Math.sqrt(n)))+", PVal = "+area+")");		//fix the alt hypo				
						break;
					case 2:
						area = 1 - Distributions.tCdf(t, n-1);
						System.out.println("tTest(alt hypo: x-bar > mu, t = "+((x-mu)/(sigma/Math.sqrt(n)))+", PVal = "+area+")");
						break;
					case 3:
						area = Distributions.tCdf(t, n-1);
						System.out.println(area);
						break;
					}
					
					break;
				case 3:
					System.out.println("Please enter the null hypothesis' proportion");
					p = scan.nextDouble();
					System.out.println("Please enter the number of successes");
					x = scan.nextDouble();
					System.out.println("Please enter the number of trials (n)");
					n = scan.nextDouble();
					p2 = x/n;
					
					z = (p2 - p)/Math.sqrt((p*(1-p))/n);
					
					System.out.println("Please enter the Alternative Hypothesis\n1. p-hat > p0\n2. p-hat != p0\n3. p-hat < p0");
					ans = scan.nextDouble();
					
					switch((int)ans) {
					case 1:
						area = 1 - Distributions.normCdf(0, 1, z);
						System.out.println("1-prop z Test(p0 = "+p+", number of successes = "+x+", number of trials = "+n+", alt hypo: p-hat > p0 z = "+z+", P-val = "+area+")");
						break;
					case 2:
						area = 2 * Distributions.normCdf(0, 1, Math.abs(z));
						System.out.println("1-prop z Test(p0 = "+p+", number of successes = "+x+", number of trials = "+n+", alt hypo: p-hat != p0 z = "+z+", P-val = "+area+")");
						break;
					case 3:
						area = Distributions.normCdf(0, 1, z);
						System.out.println("1-prop z Test(p0 = "+p+", number of successes = "+x+", number of trials = "+n+", alt hypo: p-hat < p0 z = "+z+", P-val = "+area+")");
						break;
					}
					break;
				case 4:
					System.out.println("Please enter the first number of successes");
					x = scan.nextDouble();
					System.out.println("Please enter the first number of trials");
					n = scan.nextDouble();
					System.out.println("Please enter the second number of successes");
					x1 = scan.nextDouble();
					System.out.println("Please enter the second number of trials");
					n2 = scan.nextDouble();
					
					p1 = x/n;
					p2 = x1/n2;
					p = (x+x1)/(n+n2);
					z = (p1 - p2)/( Math.sqrt((p*(1-p))*(Math.pow(n2, -1) + Math.pow(n, -1))));
					
					System.out.println("Please enter the Alternative Hypothesis\n1. p1 > p2\n2. p1 != p2\n3. p1 < p2");			
					ans = scan.nextDouble();
					
					switch((int)ans) {
					case 1:
						area = 1 - Distributions.normCdf(0, 1, z);
						System.out.println("2-prop zTest(p1 = "+p1+", p2 = "+p2+", p = "+p+", z = "+z+", alt hypo: p1 > p2, P-val = "+area+")");
						break;
					case 2:
						area = 2 * Distributions.normCdf(0, 1, z);
						System.out.println("2-prop zTest(p1 = "+p1+", p2 = "+p2+", p = "+p+", z = "+z+", alt hypo: p1 != p2, P-val = "+area+")");
						break;
					case 3:
						area = Distributions.normCdf(0, 1, z);
						System.out.println("2-prop zTest(p1 = "+p1+", p2 = "+p2+", p = "+p+", z = "+z+", alt hypo: p1 < p2, P-val = "+area+")");
						break;
					}
					break;
				case 5:
					System.out.println("Please enter sigma1");
					sigma = scan.nextDouble();
					System.out.println("Please enter sigma2");
					sigma2 = scan.nextDouble();
					System.out.println("Please enter the first x-bar");
					x = scan.nextDouble();
					System.out.println("Please enter the first sample size");
					n = scan.nextDouble();
					System.out.println("Please enter the second x-bar");
					x1 = scan.nextDouble();
					System.out.println("Please enter the second sample size");
					n1 = scan.nextDouble();
					
					z = (x-x1)/(Math.sqrt((Math.pow(sigma, 2)/n) + (Math.pow(sigma2, 2)/n1)));
					System.out.println(z);
					System.out.println("Please enter the Alternative Hypothesis\n1. mu1 > mu2\n2. mu1 != mu2\n3. mu1 < mu2");
					ans = scan.nextDouble();
					
					switch((int)ans) {
					case 1:
						area = 1 - Distributions.normCdf(0, 1, z);
						System.out.println("2-samp zTest(x1 = "+x+", sigma1 = "+sigma+", n1 = "+n+", x2 = "+x1+", sigma2 = "+sigma2+", n2 = "+n1+"m alt hypo mu1 > mu2, P-val = "+area+")");
						break;
					case 2:
						area = 2 * Distributions.normCdf(0, 1, z);
						System.out.println("2-samp zTest(x1 = "+x+", sigma1 = "+sigma+", n1 = "+n+", x2 = "+x1+", sigma2 = "+sigma2+", n2 = "+n1+"m alt hypo mu1 != mu2, P-val = "+area+")");
						break;
					case 3:
						area = Distributions.normCdf(0, 1, z);
						System.out.println("2-samp zTest(x1 = "+x+", sigma1 = "+sigma+", n1 = "+n+", x2 = "+x1+", sigma2 = "+sigma2+", n2 = "+n1+"m alt hypo mu1 < mu2, P-val = "+area+")");
						break;
					}
					
					break;
				case 6:
					System.out.println("Please enter the first x-bar");
					x = scan.nextDouble();
					System.out.println("Please enter the first sigma");
					sigma = scan.nextDouble();
					System.out.println("Please enter the first sample size");
					n1 = scan.nextDouble();
					System.out.println("Please enter the second x-bar");
					x1 = scan.nextDouble();
					System.out.println("Please enter the second sigma");
					sigma2 = scan.nextDouble();
					System.out.println("Please enter the second sample size");
					n2 = scan.nextDouble();
                    t = (x-x1) / Math.sqrt( ((sigma * sigma)/n1) + ((sigma2 * sigma2)/n2));
                    df = Distributions.getdf2Samp(sigma, sigma2, n1, n2);
                    System.out.println("Please enter the Alternative Hypothesis\n1. mu1 > mu2\n2. mu1 != mu2\n3. mu1 < mu2");
                    ans = scan.nextDouble();

                    switch((int)ans) {
                        case 1:
                            area  = 1 - Distributions.tCdf(t, df);
                            System.out.println("2-samp tTest(x1 = "+x+", sigma1 = "+sigma+", n1 = "+n1+", x2 = "+x1+", sigma2 = "+sigma2+", n2 = "+n2+"m alt hypo mu1 > mu2, P-val = "+area+")");
                            break;
                        case 2:
                            area = 2 * Distributions.tCdf(t, df);
                            System.out.println("2-samp tTest(x1 = "+x+", sigma1 = "+sigma+", n1 = "+n1+", x2 = "+x1+", sigma2 = "+sigma2+", n2 = "+n2+"m alt hypo mu1 != mu2, P-val = "+area+")");
                            break;
                        case 3:
                            area = Distributions.tCdf(t, df);
                            System.out.println("2-samp tTest(x1 = "+x+", sigma1 = "+sigma+", n1 = "+n1+", x2 = "+x1+", sigma2 = "+sigma2+", n2 = "+n2+"m alt hypo mu1 < mu2, P-val = "+area+")");
                            break;
                    }
					break;
				case 7:
					System.out.println("Please enter the length of the observed list");
					Double[] obs = new Double[scan.nextInt()];
					System.out.println("Enter every entry of the observed list one at a time");
					for(int i = 0; i<obs.length;i++)
						obs[i] = scan.nextDouble();

					Double [] exp = new Double[obs.length];
					System.out.println("Enter every entry of the expected list one at a time");
					for(int i = 0; i<exp.length;i++){
						exp[i] = scan.nextDouble();
					}

					chisq = 0.0;
					for(int i = 0; i<exp.length;i++)
						chisq += Math.pow(obs[i] - exp[i],2)/exp[i];
					System.out.println("Please enter the degrees of freedom");
					df = scan.nextDouble();
					area = 1 - Distributions.chisqCdf(chisq, df);
					System.out.print("X^2 GOF Test(obs = {");
					for(int i = 0; i<obs.length-1;i++)
						System.out.print(obs[i]+", ");
					System.out.print(obs[obs.length-1]+"}, exp = {");
					for(int i = 0; i<exp.length-1;i++)
						System.out.print(exp[i]+", ");
					System.out.print("}, df = "+df+", Chi Sq = "+chisq+", P-val = "+area+")");
					System.out.println();
					break;
				case 8:
					System.out.println("Please enter the number of rows");
					row = scan.nextInt();
					System.out.println("Please enter the number of columns");
					col = scan.nextInt();
					df = (row-1)*(col-1);
					Double[][] observed = new Double[row][col];
					for(int i = 0; i<row; i++){
						System.out.println("Please enter the values in row "+(i+1)+" one at a time");
						for(int y = 0;y<col ;y++)
							observed[i][y] = scan.nextDouble();
					}
					total = 0.0;
					for(int i=0;i<row;i++){
					    for(int y = 0;y<col;y++){
					    	total+=observed[i][y];
					    }
                    }

					Double[][] expected = new Double[row][col];

                    for(int xco=0;xco<expected.length;xco++){
                        for(int yco = 0;yco<expected[0].length;yco++){
                            expected[xco][yco] = (Distributions.getRowExp(observed, xco)*Distributions.getColExp(observed, yco))/total;
                        }
                    }
                    chisq = 0.0;
                    for(int i=0;i<row;i++){
                        for(int y = 0;y<col;y++){
                            chisq += Math.pow(observed[i][y] - expected[i][y], 2) / expected[i][y];
                        }
                    }
                    area = 1 - Distributions.chisqCdf(chisq, df);
                    System.out.println("X^2 Two Way Test(X^2 = "+chisq+", df = "+df+", P-val = "+area+")");
					break;
				}
				break;
			}
		}
	}
}