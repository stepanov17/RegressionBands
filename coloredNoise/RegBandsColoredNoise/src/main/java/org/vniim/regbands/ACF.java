
package org.vniim.regbands;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;


public class ACF {
    
    private final NoiseGenerator NG;
    private final static StandardDeviation STDEV = new StandardDeviation();
    
    private final static int LEN = 1024 * 1024;
    
    public ACF(double alpha) {
        
        NG = new NoiseGenerator(LEN, alpha);
    }

    
    private double[] getACF(int len) {

        int k = 1_000_000;
        
        if (len > LEN - k) { throw new RuntimeException("len is too large"); }
        
        double kk = 1. / k;
        
        double acf[] = new double[len];
        
        double noise[] = NG.ColoredNoise();
        
        double d1 = STDEV.evaluate(noise, 0, k);
        
        for (int i = 0; i < len; ++i) {
            
            double p = 0.;
            for (int j = 0; j < k; ++j) {
                p += noise[j] * noise[i + j];
            }
            acf[i] = kk * p / (d1 * STDEV.evaluate(noise, i, k + i));
        }
        
        return acf;
    }
    
    private double[] getACF(int len, int nAvg) {
        
        double r[] = new double[len];
        for (int l = 0; l < len; ++l) { r[l] = 0.; }; // don't need, but just in case...
        
        for (int i = 0; i < nAvg; ++i) {
            
            double f[] = getACF(len);
            for (int l = 0; l < len; ++l) {
                r[l] += f[l];
            }
            System.out.println(">> averaging: " + (i + 1) + " / " + nAvg);
        }
        
        for (int l = 0; l < len; ++l) { r[l] /= nAvg; }
        
        for (int l = 1; l < len; ++l) { r[l] = Math.min(r[l], r[l - 1]); }
        
        return r;
    }
    
    
    public static void main(String[] args) {
        
        double alpha = 0.5;
        int len = 1500, nAvg = 500;
        
        double r[] = (new ACF(alpha)).getACF(len, nAvg);

        System.out.println(">> alpha = " + alpha);
        for (double v: r) {
            System.out.println(String.format("%.4f", v));
        }
    }
}
