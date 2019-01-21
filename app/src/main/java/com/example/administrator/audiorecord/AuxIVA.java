package com.example.administrator.audiorecord;

import android.util.Log;

import static java.lang.Math.cosh;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.Math.tanh;

public class AuxIVA {

    double[][][] reSTFTin, imSTFTin, reSTFTout, imSTFTout;
    double[][][] reW, imW;
    double[][] r, G_r;
    double[][][][] reV, imV;

    int nSrc, nItr;
    int nChannels, nFrames, nFreq;
    boolean isProjBack;
    String contrastFunc;
    double C, m;

    AuxIVA(double[][][] reSTFTin, double[][][] imSTFTin,
           int nItr, boolean isProjBack,
           double[][][] reW0, double[][][] imW0,
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

        this.nChannels = reSTFTin.length;
        this.nFrames = reSTFTin[0].length;
        this.nFreq = reSTFTin[0][0].length;

        this.nSrc = nChannels;  // force equal number of sources and sensors
        // using nSrc as a separate variable to allow future addition of the overdetermined case

        this.reSTFTin = new double[this.nChannels][this.nFrames][this.nFreq];
        this.imSTFTin = new double[this.nChannels][this.nFrames][this.nFreq];


        // preemptively prevents the original STFT data from being modified
        Log.i("DEBUG", "Creating a duplicate of STFT data");
        for (int c = 0; c < this.nChannels; c++) {
            for (int i = 0; i < this.nFrames; i++) {
                System.arraycopy(reSTFTin[c][i], 0, this.reSTFTin[c][i], 0, this.nFreq);
                System.arraycopy(imSTFTin[c][i], 0, this.imSTFTin[c][i], 0, this.nFreq);
            }
        }

        // initializing the demixing matrix
        Log.i("DEBUG", "Initializing demixing matrix");
        if ((reW0 == null) || (imW0 == null)) {
            Log.i("DEBUG", "Initializing with identity");
            this.reW = new double[this.nSrc][this.nChannels][this.nFreq];
            this.imW = new double[this.nSrc][this.nChannels][this.nFreq];

            for (int f = 0; f < this.nFreq; f++) {
                for (int s = 0; s < this.nSrc; s++) {
                    for (int c = 0; c < this.nChannels; c++) {
                        if (s == c) {
                            this.reW[s][c][f] = 1.0;
                        }
                    }
                }
            }
        } else {
            Log.i("DEBUG", "Initializing with specified value");
            this.reW = reW0;
            this.imW = imW0;
        }

        // initializing the separated source matrix
        this.reSTFTout = new double[this.nSrc][this.nFrames][this.nFreq];
        this.imSTFTout = new double[this.nSrc][this.nFrames][this.nFreq];

        // initializing contrast function and related variables
        this.contrastFunc = contrastFunc;

        this.r = new double[this.nSrc][this.nFrames];
        this.G_r = new double[this.nSrc][this.nFrames];
        this.reV = new double[this.nSrc][this.nFreq][this.nChannels][this.nChannels];
        this.imV = new double[this.nSrc][this.nFreq][this.nChannels][this.nChannels];
    }

    public void run() {
        for (int runCount = 0; runCount < nItr; runCount++) {

            // demixing
            demix(this.reSTFTout, this.imSTFTout,
                    this.reSTFTin, this.imSTFTin,
                    this.reW, this.imW);

            // calculate r and G/r
            for (int s = 0; s < nSrc; s++) {
                for (int i = 0; i < nFrames; i++) {
                    double sum = 0;
                    for (int f = 0; f < nFreq; f++) {
                        sum += pow(reSTFTout[s][i][f], 2) + pow(imSTFTout[s][i][f], 2);
                    }
                    r[s][i] = sqrt(sum);

                    G_r[s][i] = calcContrastFunc(contrastFunc, "df", r[s][i]) / r[s][i];
                }
            }


            // calculate weighted covariance matrix V
            for (int s = 0; s < nSrc; s++) {
                for (int f = 0; f < nFreq; f++) {
                    for (int i = 0; i < nChannels; i++) {
                        for (int j = 0; j < nChannels; j++) {
                            for (int t = 0; t < nFrames; t++) {
                                reV[s][f][i][j] += complexMultiply(
                                        reSTFTin[i][t][f], imSTFTin[i][t][f],
                                        reSTFTin[j][t][f], -imSTFTin[j][t][f])[0]
                                        * G_r[s][t];
                                imV[s][f][i][j] += complexMultiply(
                                        reSTFTin[i][t][f], imSTFTin[i][t][f],
                                        reSTFTin[j][t][f], -imSTFTin[j][t][f])[1]
                                        * G_r[s][t];
                                ;
                            }
                            reV[s][f][i][j] = reV[s][f][i][j] / (double) nFrames;
                            imV[s][f][i][j] = imV[s][f][i][j] / (double) nFrames;
                        }
                    }
                }
            }

            //updating the demixing matrix
            for (int s = 0; s < nSrc; s++) {
                double[][][] reWV = new double[nFreq][nSrc][nChannels];
                double[][][] imWV = new double[nFreq][nSrc][nChannels];
                for (int f = 0; f < nFreq; f++) {
                    for (int i = 0; i < nSrc; i++) {
                        for (int j = 0; j < nChannels; j++) {
                            for (int k = 0; k < nChannels; k++) {
                                reWV[f][i][j] += complexMultiply(
                                        reW[i][k][f], -imW[i][k][f],
                                        reV[s][f][k][j], imV[s][f][k][j])[0];
                                imWV[f][i][j] += complexMultiply(
                                        reW[i][k][f], -imW[i][k][f],
                                        reV[s][f][k][j], imV[s][f][k][j])[1];
                            }
                        }
                    }



                    if (nSrc == 2) {
                        double[][][] invWV = new double[nChannels][nSrc][2];
                        invWV = complex2Dinv(reWV[f], imWV[f]);
                    }
                }

                //inverting WV

            }

            if (isProjBack) {
                //project back
            }

        }
    }

    private double[][][] complex2Dinv(double[][] real, double[][] imag) {
        double reDet, imDet;

        double[][][] out = new double[real.length][real.length][2];

        reDet = real[0][0] * real[0][0] - real[0][1] * real[1][0];
        imDet = imag[0][0] * imag[0][0] - imag[0][1] * imag[1][0];

        double[] invDet;

        invDet = complexReciprocal(reDet, imDet);

        out[0][0][0] = real[1][1];
        out[0][1][0] = -real[1][0];
        out[1][0][0] = -real[0][1];
        out[1][1][0] = real[0][0];

        out[0][0][1] = imag[1][1];
        out[0][1][1] = -imag[1][0];
        out[1][0][1] = -imag[0][1];
        out[1][1][1] = imag[0][0];

        for (int i = 0; i < real.length; i++) {
            for (int j = 0; j < real.length; j++) {
                out[i][j] = complexMultiply(out[i][j][0], out[i][j][1],
                        invDet[0], invDet[1]);
            }
        }

        return out;
    }

    private double[] complexReciprocal(double real, double imag) {
        double[] out = new double[2];

        double mod = pow(real, 2.0) + pow(imag, 2.0);

        out[0] = real / mod;
        out[1] = -imag / mod;

        return out;
    }

        /*
    private double[][] matrix2DMultiply(double[][] A, double[][] B) {
        int inner;
        if (A[0].length != B.length) {
            throw new RuntimeException("Invalid dimensions");
        } else {
            inner = A[0].length
        }

        int outRow = A.length;
        int outCol = B[0].length;
        double[][] out = new double[outRow][outCol];

        for (int i = 0; i < outRow; i++) {
            for (int j = 0; j < outCol; j++) {
                for (int k = 0; k < inner; k++){
                    out[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return out;
    }
    */

    private double[] complexMultiply(double reA, double imA, double reB, double imB) {
        double[] out = new double[2];
        // (a + bi) * (c + di) = (ac - bd) + i (ad + bc)
        out[0] = reA * reB - imA * imB;
        out[1] = reA * imB + imA * reB;
        return out;
    }

    private void demix(double[][][] reY, double[][][] imY, double[][][] reX, double[][][] imX, double[][][] reW, double[][][] imW) {

        // Y(f) = X(f) * W*(f)
        // X: (nChannels by nFrames by nFreq)
        // W: (nSrc by nChannels by nFreq)
        //(frame by channel) * (channel by source)

        for (int f = 0; f < nFreq; f++) {
            for (int i = 0; i < this.nSrc; i++) {
                for (int j = 0; j < this.nFrames; j++) {
                    reY[j][i][f] = 0.0;
                    imY[j][i][f] = 0.0;
                    for (int k = 0; k < this.nChannels; k++) {
                        double[] out = new double[2];
                        out = complexMultiply(reW[i][k][f], -imW[i][k][f], reX[k][j][f], imX[k][j][f]);
                        reY[j][i][f] += out[0];
                        imY[j][i][f] += out[1];
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


}
