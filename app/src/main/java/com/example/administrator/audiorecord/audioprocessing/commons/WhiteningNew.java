package com.example.administrator.audiorecord.audioprocessing.commons;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.util.Arrays;
import java.util.stream.IntStream;

import static java.lang.Math.pow;

public class WhiteningNew {

    private int nChannels, nFrames, nFreq;
    private Complex inFrames;

    private Complex[][][] STFT;

    private Array2DRowFieldMatrix<Complex>[] STFTwhite;

    private final double EPSILON = 1.0e-15;

    /*
    IMPORTANT:

    STFTin is nChannels by nFrames by nFreq to match STFT's output

    However,

    STFTout is nFreq by nChannels by nFrames to simplify further calculations
     */

    public WhiteningNew(Complex[][][] STFTin) {
        nChannels = STFTin.length;
        nFrames = STFTin[0].length;
        nFreq = STFTin[0][0].length;

        STFT = new Complex[nFreq][nChannels][nFrames];

        /* swapping axes */
        for (int c = 0; c < nChannels; c++) {
            for (int t = 0; t < nFrames; t++) {
                for (int f = 0; f < nFreq; f++) {
                    STFT[f][c][t] = STFTin[c][t][f];
                }
            }
        }

        inFrames = Complex.ONE.divide(nFrames);

        STFTwhite = new Array2DRowFieldMatrix[nFreq];

        Arrays.fill(STFTwhite, new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nFrames));
    }

    public void run() {

        Log.i("DEBUG", "Whitening");

        IntStream.range(0, nFreq).parallel().forEach(f -> { //temporarily removing  to debug
            //Log.i("DEBUG", "freq = " + f);
            STFTwhite[f] = whitenEachBin(STFT[f]);
            //checkWhiteness(STFTwhite[f]);
        });
    }

    private void checkWhiteness(Array2DRowFieldMatrix<Complex> A) {

        Array2DRowFieldMatrix<Complex> checkCov;

        checkCov = (Array2DRowFieldMatrix<Complex>) A.multiply((Array2DRowFieldMatrix<Complex>) conjugate(A).transpose()).scalarMultiply(inFrames);

        Log.i("DEBUG", "checkCov: " + checkCov.toString());
    }

    public Array2DRowFieldMatrix<Complex>[] getWhitenedMatrixArray() {
        return STFTwhite;
    }

    private Array2DRowFieldMatrix<Complex> whitenEachBin(Complex[][] Xf) {
        /* NORMALIZING */

        Complex[][] Xf_copy = new Complex[nChannels][nFrames];

        for (int c = 0; c < nChannels; c++) {
            System.arraycopy(Xf[c], 0, Xf_copy[c], 0, nFrames);
        }


        Complex[][] Xf_conj = new Complex[nChannels][nFrames];
        Complex sum;

        for (int c = 0; c < nChannels; c++) {
            sum = Complex.ZERO;
            for (int t = 0; t < nFrames; t++) {
                sum = sum.add(Xf_copy[c][t]);
            }
            Complex average = sum.divide(nFrames);

            for (int t = 0; t < nFrames; t++) {
                Xf_copy[c][t] = Xf_copy[c][t].subtract(average);
                Xf_conj[c][t] = Xf_copy[c][t].conjugate();
            }
        }

        /* Find covariance matrix*/

        Array2DRowFieldMatrix<Complex> X = new Array2DRowFieldMatrix<>(Xf_copy);
        Array2DRowFieldMatrix<Complex> X_H = (Array2DRowFieldMatrix<Complex>) new Array2DRowFieldMatrix<>(Xf_conj).transpose();
        Array2DRowFieldMatrix<Complex> cov = (Array2DRowFieldMatrix<Complex>) X.multiply(X_H).scalarMultiply(inFrames);

        ComplexSingularValueDecomposition csvd = new ComplexSingularValueDecomposition(cov, true);

        Array2DRowFieldMatrix<Complex> S_isqrt = csvd.getPoweredS(-0.5);
        Array2DRowFieldMatrix<Complex> U_H = csvd.getUH();

        Log.i("DEBUG", "S^(-1/2): " + S_isqrt.toString());
        Log.i("DEBUG", "U^H: " + U_H.toString());

        /* Whitening */

        Array2DRowFieldMatrix<Complex> out = S_isqrt.multiply(U_H).multiply(X);
/*

        Array2DRowFieldMatrix<Complex> checkCov = (Array2DRowFieldMatrix<Complex>) out.multiply(conjugate(out).transpose()).scalarMultiply(inFrames);
        Log.i("DEBUG", "checkCov: " + checkCov.toString());
*/

        return out;
    }

    private Array2DRowFieldMatrix<Complex> conjugate(Array2DRowFieldMatrix<Complex> z) {

        Array2DRowFieldMatrix<Complex> zconj = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), z.getRowDimension(), z.getColumnDimension());

        for (int r = 0; r < z.getRowDimension(); r++) {
            for (int c = 0; c < z.getColumnDimension(); c++) {
                zconj.setEntry(r, c, z.getEntry(r, c).conjugate());
            }
        }

        return zconj;
    }
}
