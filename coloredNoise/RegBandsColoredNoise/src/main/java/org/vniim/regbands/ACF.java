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

        long t0 = System.nanoTime();

        for (int i = 0; i < nAvg; ++i) {

            double f[] = getACF(len);
            for (int l = 0; l < len; ++l) {
                r[l] += f[l];
            }

            long t = System.nanoTime();
            double dt = 1.e-9 * (t - t0);

            System.out.println(">> averaging: " + (i + 1) + " / " + nAvg + ", " + Math.round(dt) + " s");
            //if ((i > 0) && (i % 50 == 1)) {
            //    for (double v: r) {
            //        System.out.print(String.format("%.4f", v / (i + 1.)) + " ");
            //    }
            //    System.out.println();
            //}
        }

        for (int l = 0; l < len; ++l) { r[l] /= nAvg; }

        for (int l = 1; l < len; ++l) { r[l] = Math.min(r[l], r[l - 1]); }

        return r;
    }


    public static void main(String[] args) {

        double alpha = 1.2;
        int len = 2048, nAvg = 2000;

        double r[] = (new ACF(alpha)).getACF(len, nAvg);

        System.out.println();
        System.out.println(">> alpha = " + alpha);
        for (double v: r) {
            System.out.println(String.format("%.4f", v));
        }
    }
}
