package org.vniim.regbands;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class NormFactor {

    private static double getFact(int n, double alpha) {

        int nSim = 2_000_000;

        NoiseGenerator NG = new NoiseGenerator(n, alpha);
        StandardDeviation STDEV = new StandardDeviation();

        double s = 0.;
        for (int i = 0; i < nSim; ++i) {
            s += STDEV.evaluate(NG.ColoredNoise());
            if (i % 100_000 == 0) { System.out.println(">> " + i); }
        }
        s /= nSim;

        return 1. / s;
    }

    private static double checkFact(int n, double alpha, double fact) {

        int nSim = 500_000;

        NoiseGenerator NG = new NoiseGenerator(n, alpha, fact);
        StandardDeviation STDEV = new StandardDeviation();

        double s = 0.;
        for (int i = 0; i < nSim; ++i) {
            s += STDEV.evaluate(NG.ColoredNoise());
            if (i % 100_000 == 0) { System.out.println(">> " + i); }
        }
        s /= nSim;

        return 1. / s;

    }

    public static void main(String[] args) {

        int n = 2048;
        double alpha = 1.;

        double k = getFact(n, alpha);

        double check = checkFact(n, alpha, k);

        System.out.println(alpha + ":  " + k + "  " + check);
    }
}
