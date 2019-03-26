package com.example.administrator.audiorecord.audioprocessing.bss;

import android.util.Log;

import com.example.administrator.audiorecord.audioprocessing.commons.ComplexSingularValueDecomposition;
import com.example.administrator.audiorecord.audioprocessing.commons.Whitening;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import java.util.Arrays;
import java.util.stream.IntStream;

import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static java.lang.Math.random;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;

public class AbsCosSimSCA {

    private Complex[][][] STFTin, STFTout;

    private int nSrc, nItr;
    private int nChannels, nFrames, nFreqs;
    private boolean isProjBack, derivCheck, isRowDecoupling;

    private final double EPSILON = 1e-8;
    private final double PROGRESS_TOLERANCE = 1e-8;

    Array2DRowFieldMatrix<Complex>[] STFTwhite;
    Array2DRowFieldMatrix<Complex>[] STFTper;

    Array2DRowFieldMatrix<Complex>[] Q_PCA;
    Array2DRowFieldMatrix<Complex>[] demix;

    Complex eta;

    final double beta = 12.5;

    double[][][] permutations;

    //Complex[][] randomInit;

    public AbsCosSimSCA(Complex[][][] STFTin, int nItr, boolean isProjBack, int nSrc, double eta, boolean derivCheck, boolean isRowDecoupling) {

        this.STFTin = STFTin;

        this.nItr = nItr;

        this.eta = new Complex(eta);

        this.derivCheck = derivCheck;
        this.isRowDecoupling = isRowDecoupling;

        this.nSrc = nSrc;

        this.nChannels = STFTin.length;
        this.nFrames = STFTin[0].length;
        this.nFreqs = STFTin[0][0].length;

        STFTout = new Complex[this.nSrc][nFrames][nFreqs];

        STFTwhite = new Array2DRowFieldMatrix[nFreqs];
        STFTper = new Array2DRowFieldMatrix[nFreqs];
        demix = new Array2DRowFieldMatrix[nFreqs];

        permutations = new double[nFreqs][nFrames][nSrc];

        /*
        randomInit = new Complex[nChannels][this.nSrc];
        */

        /*
        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                if (c == s) {
                    randomInit[c][s] = Complex.ONE;
                } else {
                    randomInit[c][s] = Complex.ZERO;
                }
            }
        }
        */
    }

    public void run() {
        whiten();

        IntStream.range(0, nFreqs).parallel().forEach(this::runThisFreq); //

        //Log.i("DEBUG", "COMPLETED ALL FREQUENCIES");

        alignPermutation();
    }

    private void whiten() {

        /*
        IMPORTANT:

        STFTin is nChannels by nFrames by nFreq to match STFT's output

        However,

        STFTout is nFreq by nChannels by nFrames to simplify further calculations
        */

        /*if (derivCheck) {

            Complex[][][] STFT = new Complex[nFreqs][nChannels][nFrames];

            for (int c = 0; c < nChannels; c++) {
                for (int t = 0; t < nFrames; t++) {
                    for (int f = 0; f < nFreqs; f++) {
                        STFT[f][c][t] = STFTin[c][t][f];
                    }
                }
            }

            for (int f = 0; f < nFreqs; f++) {
                STFTwhite[f] = new Array2DRowFieldMatrix<>(STFT[f]);
            }

        } else {*/
        WhiteningNew whitening = new WhiteningNew(STFTin);
        whitening.run();
        STFTwhite = whitening.getWhitenedMatrixArray(); // nFreqs x (nChannels x nFrames)

        Q_PCA = whitening.getQ_PCA();

    }

    private void runThisFreq(int f) {

        //Log.i("DEBUG", "bin = " + f);

        /* normalizing the columns of X */

        Array2DRowFieldMatrix<Complex> X_bar = STFTwhite[f];

        //Log.i("DEBUG", "X = " + X_bar);

        double colNorm;

        for (int t = 0; t < nFrames; t++) {

            colNorm = 0.0;

            for (int c = 0; c < nChannels; c++) {
                colNorm += pow(X_bar.getEntry(c, t).abs(), 2.0);
            }

            colNorm = sqrt(colNorm + EPSILON);

            for (int c = 0; c < nChannels; c++) {
                X_bar.multiplyEntry(c, t, Complex.ONE.divide(colNorm));
            }
        }

        //Log.i("DEBUG", "X_bar = " + X_bar);
        /* checking */
        /*
            double checkColNorm;

            for (int t = 0; t < nFrames; t++) {

                checkColNorm = 0.0;

                for (int c = 0; c < nChannels; c++) {
                    checkColNorm += pow(X_bar.getEntry(c, t).abs(), 2.0);
                }

                checkColNorm = sqrt(checkColNorm);

                //Log.i("DEBUG", "frame " + t + " norm = " + checkColNorm);
            }
        */

        /* Initializing mixing matrix */

        Complex[][] randomInit = new Complex[nChannels][nSrc];

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                randomInit[c][s] = new Complex(random(), random());
            }
        }


        /* decouple the rows of A */

        Array2DRowFieldMatrix<Complex> A = new Array2DRowFieldMatrix<>(randomInit);
