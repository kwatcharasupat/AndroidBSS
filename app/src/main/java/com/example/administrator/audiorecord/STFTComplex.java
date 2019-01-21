package com.example.administrator.audiorecord;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.FastFourierTransformer;

import static java.lang.Math.PI;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static org.apache.commons.math3.transform.DftNormalization.STANDARD;
import static org.apache.commons.math3.transform.TransformType.FORWARD;
import static org.apache.commons.math3.transform.TransformType.INVERSE;

class STFTComplex {

    /*
    Multichannel Short Time Fourier Transform and Inverse Short Time Fourier Transform for Java
    2018 (c) Karn Watcharasupat, GNU Public License v3.0
     */

    private int winLen, hopSize, nOverlap, nFrames;
    private int nSamples, nChannels, nFreq;
    private int paddedLength;

    private double[] window;
    private double[] iswin;

    private double[][] input;
    private Complex[][][] output;

    private double[][] invOutput;

    STFTComplex() {
    }

    void stftm(double[][] input, int winLen, int nOverlap, String winFunc) {

        /*
        Multichannel Short Time Fourier Transform

        Parameters
        ----------
        input: double[][] (nChannels by nSamples)
            The raw multichannel signal
        winLen: int
            The window length
        nOverlap: int
            The number of overlapping samples
        winFunc: String
            The window function
        Output
        ------
        output: Complex[][][] (nChannels by nFrames by nFreq)
            The STFT representation of input signal
            call getSTFT() to obtain the matrix
         */

        Log.i("DEBUG", "STFT running");

        this.input = input;

        this.nChannels = input.length;
        this.nSamples = input[0].length;  // time is in the second dimension

        this.winLen = winLen;
        this.nOverlap = nOverlap;
        this.hopSize = winLen - nOverlap;

        this.nFrames = (int) ceil(nSamples / nOverlap) + (int) ceil(nOverlap / hopSize);
        this.paddedLength = nFrames * hopSize + 2 * nOverlap;

        double[][] paddedInput = new double[nChannels][paddedLength];

        this.nFreq = winLen / 2 + 1;

        this.output = new Complex[nChannels][nFrames][this.nFreq];

        // padding
        for (int c = 0; c < nChannels; c++) {
            System.arraycopy(this.input[c], 0, paddedInput[c], nOverlap, nSamples);
        }

        this.window = new double[winLen];
        getSTFTwindow(winFunc);

        setScalingVec(true);    // true for forward fft
        for (int c = 0; c < nChannels; c++) {
            for (int t = 0; t < nFrames; t++) {
                double[] frame = new double[winLen];
                double[] wFrame = new double[winLen];
                Complex[] inputFrame = new Complex[winLen];

                System.arraycopy(paddedInput[c], t * hopSize, frame, 0, winLen);

                for (int j = 0; j < winLen; j++) {
                    wFrame[j] = frame[j] * window[j] * iswin[t * hopSize + j];
                }

                FastFourierTransformer FFT = new FastFourierTransformer(STANDARD);

                output[c][t] = FFT.transform(wFrame, FORWARD);
            }
        }
        Log.i("DEBUG", "STFT completed");
    }

    void istftm(Complex[][][] STFT, int winLen, int nOverlap, String winFunc, int nSamples) {
        /*
        Multichannel Inverse Short Time Fourier Transform

        Parameters
        ----------
        input: Complex[][][] (nChannels by nSamples by nFreq)
            The processed complex STFT data
        winLen: int
            The window length
        nOverlap: int
            The number of overlapping samples
        winFunc: String
            The window function
        nSamples: int
            Number of samples in the original data
        Output
        ------
        output: double[][] (nChannels by nSamples)
            Reconstructed real signal
         */

        Log.i("DEBUG", "Inverse STFT running");

        this.nChannels = STFT.length;
        int inv_nFrames = STFT[0].length;

        this.winLen = winLen;
        this.nOverlap = nOverlap;
        this.hopSize = winLen - nOverlap;

        this.nSamples = nSamples;  // time is in the second dimension

        this.paddedLength = inv_nFrames * hopSize + 2 * nOverlap;

        this.window = new double[winLen];
        double[][] invPaddedOutput = new double[nChannels][paddedLength];

        int nFreq = winLen / 2 + 1;

        this.input = new double[nChannels][nSamples];

        getSTFTwindow(winFunc);

        setScalingVec(false);   // false for inverse fft

        this.invOutput = new double[nChannels][nSamples];

        for (int c = 0; c < nChannels; c++) {
            for (int i = 0; i < inv_nFrames; i++) {
                Complex[] inputFrame = new Complex[winLen];
                System.arraycopy(STFT[c][i], 0, inputFrame, 0, nFreq);

                for (int j = 0; j < winLen / 2 - 1; j++) {
                    inputFrame[nFreq + j] = STFT[c][i][nFreq - 2 - j].conjugate();
                }

                FastFourierTransformer FFT = new FastFourierTransformer(STANDARD);

                Complex[] tempFrame;
                tempFrame = FFT.transform(inputFrame, INVERSE);

                double[] fframe = new double[winLen];
                for (int t = 0; t < winLen; t++) {
                    fframe[t] = tempFrame[t].getReal();
                }

                double[] wFrame = new double[winLen];

                for (int j = 0; j < winLen; j++) {
                    wFrame[j] = window[j] * fframe[j] * iswin[i * hopSize + j];
                    invPaddedOutput[c][i * hopSize + j] = invPaddedOutput[c][i * hopSize + j] + wFrame[j];
                }
            }

            System.arraycopy(invPaddedOutput[c], nOverlap, this.invOutput[c], 0, nSamples);

            Log.i("DEBUG", "Inverse STFT completed");
        }
    }

    private void setScalingVec(boolean isFwd) {

        final double EPSILON = 1e-8;
        double[] swin = new double[paddedLength];
        iswin = new double[paddedLength];

        for (int i = 0; i < nFrames; i++) {
            for (int j = 0; j < winLen; j++) {
                swin[i * hopSize + j] = swin[i * hopSize + j] + window[j] * window[j];
            }
        }

        for (int i = 0; i < paddedLength; i++) {
            if (swin[i] == 0) {
                swin[i] = EPSILON;
            } else {
                if (isFwd) {
                    swin[i] = sqrt(swin[i] * winLen);
                } else {
                    swin[i] = sqrt(swin[i] / winLen);
                }
            }
            iswin[i] = 1.0 / swin[i];
        }
    }

    private void getSTFTwindow(String winFunc) {
        switch (winFunc) {
            case "sine":
                // sine window
                for (int i = 0; i < winLen; i++) {
                    window[i] = sin(PI / winLen * ((double) i + 0.5));
                }
                break;
            case "hann":
                // hann window
                for (int i = 0; i < winLen; i++) {
                    window[i] = 0.5 * (1 - cos(2 * PI * i / (winLen - 1)));
                }

        }
    }

    Complex[][][] getSTFT() {
        return output;
    }

    double[][] getRealSigFromInvSTFT() {
        return this.invOutput;
    }

    int get_nFrames() {
        return nFrames;
    }

    int get_nFreq() {
        return nFreq;
    }
}
