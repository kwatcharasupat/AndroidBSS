package com.example.administrator.audiorecord;

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

        this.nChannels = reSTFTin.length;
        this.nFrames = reSTFTin[0].length;
        this.nFreq = reSTFTin[0][0].length;

        this.nSrc = nChannels;  // force equal number of sources and sensors
        // using nSrc as a separate variable to allow future addition of the overdetermined case

        this.reSTFTin = new double[this.nChannels][this.nFrames][this.nFreq];
        this.imSTFTin = new double[this.nChannels][this.nFrames][this.nFreq];

        // preemptively prevents the original STFT data from being modified
        for (int c = 0; c < this.nChannels; c++) {
            for (int i = 0; i < this.nFrames; i++) {
                System.arraycopy(reSTFTin[c][i], 0, this.reSTFTin[c][i], 0, this.nFreq);
                System.arraycopy(imSTFTin[c][i], 0, this.imSTFTin[c][i], 0, this.nFreq);
            }
        }

        // initializing the demixing matrix
        if ((reW0 == null) || (imW0 == null)) {
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
        this.reV = new double[this.nChannels][this.nChannels][this.nSrc][this.nFreq];
        this.imV = new double[this.nChannels][this.nChannels][this.nSrc][this.nFreq];
    }

    public void run() {
        for (int runCount = 0; runCount < nItr; runCount++) {

            demix(this.reSTFTout, this.imSTFTout,
                    this.reSTFTin, this.imSTFTin,
                    this.reW, this.imW);

            if (isProjBack) {
                //project back
            }

            for (int s = 0; s < nSrc; s++) {
                for (int i = 0; i < nFrames; i++) {
                    double sum = 0;
                    for (int f = 0; f < nFreq; f++) {
                        sum += pow(reSTFTout[s][i][f], 2) + pow(imSTFTout[s][i][f], 2);
                    }
                    r[s][i] = sqrt(sum);

                    G_r[s][i] = calcContrastFunc(contrastFunc,"df",r[s][i]) / r[s][i];
                }
            }


        }
    }

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
}
