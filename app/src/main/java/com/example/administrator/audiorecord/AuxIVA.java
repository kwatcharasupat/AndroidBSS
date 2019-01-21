package com.example.administrator.audiorecord;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;

import static java.lang.Math.cosh;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.Math.tanh;

public class AuxIVA {

    Complex[][][] STFTin, STFTout;
    Complex[][][] W;
    double[][] r, G_r;
    Complex[][][][] V;

    int nSrc, nItr;
    int nChannels, nFrames, nFreq;
    boolean isProjBack;
    String contrastFunc;
    double C, m;

    AuxIVA(Complex[][][] STFTin,
           int nItr, boolean isProjBack,
           Complex[][][] W0,
           String contrastFunc, double C, double m) {
        /*
        Implementation of AuxIVA algorithm for BSS presented in
        N. Ono, *Stable and fast update rules for independent vector analysis based
        on auxiliary function technique*, Proc. IEEE, WASPAA, 2011.

        Reference
        ---------
        Robin Scheibler (2018), Blind Source Separation using Independent Vector Analysis with Auxiliary Function

        Parameters
        ----------
        reSTFTin: double[][][] (nChannels by nFrames by nFreq)
            real part of the STFT representation of observed signal

        imSTFTin: double[][][] (nChannels by nFrames by nFreq)
            imaginary part of the STFT representation of observed signal

        nItr: int
            number of iterations

        isProjBack: boolean
            scaling on the first microphone by back projection

        demix_init: double[][][] (nSrc by nChannels by nFreq)

        ContrastFunc:  String
            contrast function

        Output
        ------
        reSTFTout: double[][][] (nChannels by nFrames by nFreq)
            real part of the STFT representation of separated signal

        imSTFTout: double[][][] (nChannels by nFrames by nFreq)
            imaginary part of the STFT representation of separated signal
         */

        Log.i("DEBUG", "AuxIVA now preparing");

        this.nChannels = STFTin.length;
        this.nFrames = STFTin[0].length;
        this.nFreq = STFTin[0][0].length;

        this.nSrc = nChannels;  // force equal number of sources and sensors
        // using nSrc as a separate variable to allow future addition of the overdetermined case

        this.STFTin = new Complex[this.nChannels][this.nFrames][this.nFreq];

        // preemptively prevents the original STFT data from being modified
        Log.i("DEBUG", "Creating a duplicate of STFT data");
        for (int c = 0; c < this.nChannels; c++) {
            for (int i = 0; i < this.nFrames; i++) {
                System.arraycopy(STFTin[c][i], 0, this.STFTin[c][i], 0, this.nFreq);
            }
        }

        // initializing the demixing matrix
        Log.i("DEBUG", "Initializing demixing matrix");
        if (W0 == null) {
            Log.i("DEBUG", "Initializing with identity");
            this.W = new Complex[this.nFreq][this.nSrc][this.nChannels];

            for (int f = 0; f < this.nFreq; f++) {
                for (int s = 0; s < this.nSrc; s++) {
                    for (int c = 0; c < this.nChannels; c++) {
                        if (s == c) {
                            this.W[f][s][c] = Complex.ONE;
                        }
                    }
                }
            }
        } else {
            Log.i("DEBUG", "Initializing with specified value");
            this.W = W0;
        }

        // initializing the separated source matrix
        this.STFTout = new Complex[this.nSrc][this.nFrames][this.nFreq];

        // initializing contrast function and related variables
        this.contrastFunc = contrastFunc;

        this.r = new double[this.nSrc][this.nFrames];
        this.G_r = new double[this.nSrc][this.nFrames];
        this.V = new Complex[this.nSrc][this.nFreq][this.nChannels][this.nChannels];
    }

    public void run() {
        for (int runCount = 0; runCount < nItr; runCount++) {

            // demixing
            demix(STFTout, STFTin, W);

            // calculate r and G/r
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    double sum = 0;
                    for (int f = 0; f < nFreq; f++) {
                        sum += pow(STFTout[s][t][f].abs(), 2.0);
                    }
                    r[s][t] = sqrt(sum);
                    G_r[s][t] = calcContrastFunc(contrastFunc, "df", r[s][t]) / r[s][t];
                }
            }

            // calculate weighted covariance matrix V
            for (int s = 0; s < nSrc; s++) {
                for (int f = 0; f < nFreq; f++) {
                    for (int i = 0; i < nChannels; i++) {
                        for (int j = 0; j < nChannels; j++) {
                            for (int t = 0; t < nFrames; t++) {
                                V[s][f][i][j].add(
                                        STFTin[i][t][f].multiply(
                                                STFTin[i][t][f].conjugate()
                                        ).multiply(G_r[s][t]));
                            }
                            V[s][f][i][j].divide((double) nFrames);
                        }
                    }
                }
            }

            //updating the demixing matrix
            for (int s = 0; s < nSrc; s++) {
                Complex[][][] WV = new Complex[nFreq][nSrc][nChannels];
                for (int f = 0; f < nFreq; f++) {
                    for (int i = 0; i < nSrc; i++) {
                        for (int j = 0; j < nChannels; j++) {
                            WV[f][i][j] = Complex.ZERO;
                            for (int k = 0; k < nChannels; k++) {
                                WV[f][i][j].add(W[f][i][k].conjugate().multiply(V[s][f][k][j]));
                            }
                        }
                    }

                    FieldMatrix<Complex> WVmat = new Array2DRowFieldMatrix<>(WV[f]);
                    FieldMatrix<Complex> WVinv = new FieldLUDecomposition<Complex>(WVmat).getSolver().getInverse();
                    W[f] = WVinv.getData();

                    Complex[] Vw = new Complex[nChannels];

                    for (int i = 0; i < nChannels; i++) {
                        Vw[i] = Complex.ZERO;
                        for (int k = 0; k < nChannels; k++) {
                            Vw[i].add(V[s][f][i][k].multiply(W[f][s][k]));
                        }
                    }

                    Complex norm = Complex.ZERO;
                    for (int k = 0; k < nChannels; k++) {
                        norm.add(W[f][s][k].conjugate().multiply(Vw[k]));
                    }

                    norm = norm.sqrt();

                    for (int i = 0; i < nChannels; i++) {
                        W[f][s][i] = W[f][s][i].divide(norm);
                    }
                }
            }
        }

        demix(STFTout, STFTin, W);
        if (isProjBack) {
            Complex[][] scale = new Complex[nFreq][nSrc];
            scale = projectBack(STFTout, STFTin[0]);

            for (int s = 0; s < nSrc; s++){
                for (int f = 0; f < nFreq; f++){
                    for (int t = 0; t < nFrames; t++){
                        STFTout[s][t][f].multiply(scale[s][f].conjugate());
                    }
                }
            }

        }
    }

    Complex[][][] getSourceEstimatesSTFT(){
        return STFTout;
    }


    private void demix(Complex[][][] Y, Complex[][][] X, Complex[][][] W) {

        // Y(f) = X(f) * W*(f)
        // Y: (nSrc by nFrames by nFreq)
        // X: (nChannels by nFrames by nFreq)
        // W: (nSrc by nChannels by nFreq)
        //(frame by channel) * (channel by source)

        for (int f = 0; f < nFreq; f++) {
            for (int i = 0; i < this.nSrc; i++) {
                for (int j = 0; j < this.nFrames; j++) {
                    Y[j][i][f] = Complex.ZERO;
                    for (int k = 0; k < this.nChannels; k++) {
                        Y[j][i][f] = W[f][i][k].conjugate().multiply(X[k][j][f]);
                    }
                }
            }
        }
    }

    private double calcContrastFunc(String func, String deriv, double r) {

        double result = 0.0;

        switch (func) {
            case "norm":
                switch (deriv) {
                    case "f":
                        result = (C * r);
                        break;
                    case "df":
                        result = C;
                        break;
                }
                break;
            case "cosh":
                switch (deriv) {
                    case "f":
                        result = (m * log(cosh(C * r)));
                        break;
                    case "df":
                        result = (C * m * tanh(C * r));
                        break;
                }
                break;
        }

        return result;
    }

    private Complex[][] projectBack(Complex[][][] out, Complex[][] ref) {
        Complex[][] scale = new Complex[nFreq][nSrc];
        Complex[][] num = new Complex[nFreq][nSrc];
        double[][] denom = new double[nFreq][nSrc];

        for (int f = 0; f < nFreq; f++) {
            for (int s = 0; s < nSrc; s++) {
                num[f][s] = Complex.ZERO;
                scale[f][s] = Complex.ZERO;
                for (int t = 0; t < nFrames; t++) {
                    num[f][s].add(ref[f][s].conjugate().multiply(out[s][t][f]));
                    denom[f][s] += (pow(out[s][t][f].abs(), 2.0));
                }
                if (denom[f][s] > 0) {
                    scale[f][s] = num[f][s].divide(denom[f][s]);
                }
            }
        }

        return scale;
    }
}
