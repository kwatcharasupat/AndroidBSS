package com.example.administrator.audiorecord.audioprocessing.commons;

import android.util.Log;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;

import java.util.Arrays;
import java.util.stream.IntStream;

import static com.example.administrator.audiorecord.audioprocessing.commons.Conjugate.transjugate;


public class Whitening {

    private int nChannels, nFrames, nFreq;
    private Complex inFrames;

    private Complex[][][] STFT;

    private Array2DRowFieldMatrix<Complex>[] STFTwhite;

    private Array2DRowFieldMatrix<Complex>[] Q_PCA;

    /*
    IMPORTANT:
    STFTin is nChannels by nFrames by nFreq to match STFT's output
    However, STFTout is nFreq by nChannels by nFrames to simplify further calculations
     */

    public Whitening(Complex[][][] STFTin) {
        nChannels = STFTin.length;
        nFrames = STFTin[0].length;
        nFreq = STFTin[0][0].length;

        inFrames = Complex.ONE.divide((double)nFrames);

        STFT = new Complex[nFreq][nChannels][nFrames];
        /* swapping axes */
        for (int c = 0; c < nChannels; c++) {
            for (int t = 0; t < nFrames; t++) {
                for (int f = 0; f < nFreq; f++) {
                    STFT[f][c][t] = STFTin[c][t][f];
                }
            }
        }

        STFTwhite = new Array2DRowFieldMatrix[nFreq];
        Arrays.fill(STFTwhite, new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nFrames));

        Q_PCA = new Array2DRowFieldMatrix[nFreq];
    }

    public void run() {
        //Log.i("DEBUG", "Whitening");
        IntStream.range(0, nFreq).parallel().forEach(f -> {
            //Log.i("DEBUG", "freq = " + f);
            STFTwhite[f] = whitenEachBin(STFT[f], f);
            //checkWhiteness(STFTwhite[f], f);
        });
    }

    private void checkWhiteness(Array2DRowFieldMatrix<Complex> A, int f) {
        Array2DRowFieldMatrix<Complex> checkCov;
        checkCov = (Array2DRowFieldMatrix<Complex>) A.multiply(transjugate(A)).scalarMultiply(inFrames);
        Log.i("DEBUG", "checkCov at " + f + ": " + checkCov.toString());
    }

    public Array2DRowFieldMatrix<Complex>[] getWhitenedMatrixArray() {
        return STFTwhite;
    }

    private Array2DRowFieldMatrix<Complex> whitenEachBin(Complex[][] Xf, int f) {
        /* normalizing */
        Complex[][] Xf_copy = SerializationUtils.clone(Xf);
        Complex[][] Xf_conj = new Complex[nChannels][nFrames];
        Complex sum;

        for (int c = 0; c < nChannels; c++) {
            sum = Complex.ZERO;
            for (int t = 0; t < nFrames; t++) {
                sum = sum.add(Xf_copy[c][t]);
            }
            Complex average = sum.divide((double)nFrames);

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
        Array2DRowFieldMatrix<Complex> U = csvd.getU();
        Array2DRowFieldMatrix<Complex> U_H = csvd.getUH();

        /* Whitening */
        Q_PCA[f] = U.multiply(S_isqrt).multiply(U_H);

        Array2DRowFieldMatrix<Complex> out = Q_PCA[f].multiply(X);
        return out;
    }

}
