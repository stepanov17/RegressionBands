
package org.vniim.regbands;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;


public class PSD {

    // check: the noise is indeed 1 / f^alpha noise, i.e. its PSD ~ 1 / f^alpha
    public static void main(String[] args) {

        int nAvg = 1_000_000;

        int n = 128;
        double alpha = 1.;

        NoiseGenerator NG = new NoiseGenerator(n, alpha);

        FastFourierTransformer FFTR = new FastFourierTransformer(DftNormalization.STANDARD);

        double psd[] = new double[n];

        for (int i = 0; i < nAvg; ++i) {

            double noise[] = NG.ColoredNoise();
            Complex t[] = FFTR.transform(noise, TransformType.FORWARD);
            for (int j = 0; j < n; ++j) {
                psd[j] += Math.pow(t[j].abs(), 2.); // PSD ~ square of FFT magnitude
            }
        }

        for (int j = 0; j < n; ++j) {
            psd[j] /= (n * nAvg);
        }

        System.out.printf("psd\tref\n");
        for (int i = 0; i < 20; ++i) {

            double psdRef = Math.pow(1. + i, -alpha);
            System.out.printf("%.3f\t%.3f\n", psd[i], psdRef);
        }
    }
}
