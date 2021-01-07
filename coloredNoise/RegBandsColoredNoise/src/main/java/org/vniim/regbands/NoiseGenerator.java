package org.vniim.regbands;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class NoiseGenerator {

    private int nearestPowerOf2(int n) {

        int p = 1;
        while (p < n) { p *= 2; }
        return p;
    }

    private final int n;
    private final int p2;
    private final double alpha;
    private final double kNorm;

    public NoiseGenerator(int n, double alpha) {

        this(n, alpha, 1.);
    }

    public NoiseGenerator(int n, double alpha, double kNorm) {

        // todo: check alpha

        this.n = n;

        p2 = nearestPowerOf2(n);
        if (p2 > n) {
            System.out.println(">> NoiseGenerator: using p2 = " + p2 + " (power of 2)");
        }

        this.alpha = alpha;
        this.kNorm = kNorm;
    }

    private final static Random RND = new Random(/*123456L*/);

    private final static FastFourierTransformer FFTR = new FastFourierTransformer(DftNormalization.STANDARD);

    // 1 / f^alpha noise
    public double[] ColoredNoise() {

        double x[] = new double[p2];
        for (int i = 0; i < p2; ++i) { x[i] = RND.nextGaussian(); }

        int n2 = p2 / 2 + 1;

        Complex t[] = FFTR.transform(x, TransformType.FORWARD);

        t[0] = Complex.valueOf(0.);

        for (int i = 1; i < n2; ++i) {
            t[i] = t[i].divide(Math.pow(i + 1., 0.5 * alpha));
            if (i < n2 - 1) { t[p2 - i] = t[i].conjugate(); }
        }

        Complex y[] = FFTR.transform(t, TransformType.INVERSE);

        double noise[] = new double[p2];
        for (int i = 0; i < p2; ++i) { noise[i] = kNorm * y[i].getReal(); }

        return Arrays.copyOfRange(noise, 0, n);
    }
}
