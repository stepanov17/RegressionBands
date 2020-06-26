package org.vniim.regbands;

import java.math.RoundingMode;
import java.text.DecimalFormat;
//import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

/**
 * @author astepanov
 */
public class CalculateK {
    
    private final static boolean COUNT_TRIALS = true;
    private final static boolean PRINT_DEBUG = false;
    //private final static int N_AVG = 20;
    
    private double zeroInBorders(double forecast[], double u[], double K) {
        
        int nx = forecast.length;
        if (nx != u.length) { throw new RuntimeException("incompatible dimensions"); }
        
        for (int i = 0; i < nx; ++i) {
            double L = forecast[i] - K * u[i];
            double H = forecast[i] + K * u[i];
            if ((L > 0.) || (H < 0.)) { return 0.; }
        }
        
        return 1.;
    }
    
    private double[] getProbs(int n, int n0, double alpha, double ACF[], double kNorm, double K[], int nSim) {
        
        ForecastCalculator fc = new ForecastCalculator(n, n0, ACF);
        double u[] = fc.getU();
        
        NoiseGenerator ng = new NoiseGenerator(n, alpha, kNorm);
        
        int nK = K.length;
        for (int i = 1; i < nK; ++i) {
            if (K[i - 1] >= K[i]) { throw new RuntimeException("K must be increasing"); }
        }
        
        double P[] = new double[nK];
        
        double chk[] = new double[nSim];
        
        for (int sim = 0; sim < nSim; ++sim) {
            
            double noise[] = ng.ColoredNoise();
            double forecast[] = fc.getForecast(noise);
            
            double dP[] = new double[nK];
            dP[0] = zeroInBorders(forecast, u, K[0]);
            P[0] += dP[0];
            
            for (int i = 1; i < nK; ++i) {
                if (dP[i - 1] > 1.e-3) {
                    dP[i] = 1.; // may speed up here as the K is increasing
                } else {
                    dP[i] = zeroInBorders(forecast, u, K[i]);
                }
                P[i] += dP[i];
            }
            
            chk[sim] = noise[n / 2];
            
            if (COUNT_TRIALS && ((sim + 1) % 10_000 == 0)) {
                System.out.println(">> " + (sim + 1) + " trials");
            }
        }
        
        double stdev = (new StandardDeviation()).evaluate(chk);
        if (Math.abs(1. - stdev) > 2.e-2) {
            //throw new RuntimeException("stdev check failed, stdev = " + stdev);
            System.err.println("stdev check failed, stdev = " + stdev);
        }
        
        for (int i = 0; i < P.length; ++i) { P[i] /= nSim; }
        return P;
    }

    private double getArg95(double x1, double y1, double x2, double y2) {
        
        double k = (y1 - y2) / (x1 - x2);
        double b = y1 - k * x1;
        return (0.95 - b) / k;
    }
    
    private double getK(int n, int n0, double alpha, double ACF[], double kNorm, double k1, double k2, int nSim) {
        
        double dk = 0.02 * (k2 - k1);
        
        int nK = 51;
        double K[] = new double[nK];
        for (int i = 0; i < nK; ++i) {
            K[i] = k1 + i * dk;
        }
        
        double P[] = getProbs(n, n0, alpha, ACF, kNorm, K, nSim);
        
        int i0 = 0;
        for (; i0 < nK; ++i0) {
            if (P[i0] > 0.95) { break; }
        }
        
        if (i0 == 0 ) { ++i0; }
        if (i0 == nK) { --i0; }
        
        double kk1 = K[i0 - 1], kk2 = K[i0];
        double p1  = P[i0 - 1], p2  = P[i0];
        
        double k95 = getArg95(kk1, p1, kk2, p2);
        
        if (PRINT_DEBUG) {
            System.out.println(">> " + kk1 + " -> " + p1);
            System.out.println(">> " + kk2 + " -> " + p2);
            System.out.println(">> k95 = " + k95);
            System.out.println("");
        }
                
        return k95;
    }

    private double getKIter(int n, int n0, double alpha, double ACF[], double kNorm, double k1, double k2, int nSim) {
        
        double k = getK(n, n0, alpha, ACF, kNorm, k1, k2, nSim);
        
        double dk = 0.05 * (k2 - k1);
        k = getK(n, n0, alpha, ACF, kNorm, k - dk, k + dk, nSim);

        return k;
    }


//*    
    public static void main(String[] args) {
        
        CalculateK calc = new CalculateK();
        
        int n0 = 650, nSim = 1_000_000;
        
        double k1 = 1.5, k2 = 2.;
        
        double data[] = NoiseData.DATA;
        int n = data.length - 2;
        double alpha = data[0];
        double kNorm = data[1];
        double ACF[] = Arrays.copyOfRange(data, 2, n + 2);
        
        System.out.println(">> alpha = " + alpha + ", n = " + n + ", n0 = " + n0);
        
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.HALF_UP);
        
        double k = calc.getKIter(n, n0, alpha, ACF, kNorm, k1, k2, nSim);
        System.out.println("k = " + df.format(k));
        
        // check
        double p[] = calc.getProbs(n, n0, alpha, ACF, kNorm, new double[]{k}, nSim);
        System.out.println("check: p = " + df.format(p[0]));
    }
//*/
    
    /*
    public static void main(String[] args) {
        
        CalculateK calc = new CalculateK();
        
        int nSim = 1_000_000;
        
        double k1 = 1.5, k2 = 2.5;
        
        double data[] = NoiseData.DATA;
        
        int n = data.length - 2;
        
        double alpha = data[0];
        double kNorm = data[1];
        
        double ACF[] = Arrays.copyOfRange(data, 2, n + 2);
        
        ArrayList<Double>  K  = new ArrayList<>();
        ArrayList<Double>  P  = new ArrayList<>();
        ArrayList<Integer> N0 = new ArrayList<>();
        
        for (int n0 = 250; n0 <= 1000; n0 += 50) {
            
            System.out.println(">>> " + n0 + " >>>");
            
            N0.add(n0);
            double k = calc.getKIter(n, n0, alpha, ACF, kNorm, k1, k2, nSim);
            K.add(k);
            double p[] = calc.getProbs(n, n0, alpha, ACF, kNorm, new double[]{k}, nSim);
            P.add(p[0]);
        }
        
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.HALF_UP);
        
        System.out.println(">>> alpha = " + alpha + ", n = " + n);
        
        System.out.println(">>> K");
        for (int i = 0; i < K.size(); ++i) {
            System.out.println(N0.get(i) + "  " + df.format(K.get(i)) + "  " + df.format(P.get(i)));
        }
    }
    */
}
