package com.example.administrator.audiorecord;

import android.util.Log;

import java.util.Arrays;

import static java.lang.Math.PI;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

public class STFT {


    private int winLen;
    private int hopSize;
    private int nOverlap;
    private int nFrames;
    private int nSamples;
    private int nChannels;
    private int paddedLength;
    private int nFreq;
    private int inv_nFreq;

    private double[] window;
    private double[] iswin;
    private double[][] input;
    private double[][] paddedInput;

    private double[][][] reOutput;
    private double[][][] imOutput;

    private double[][] invPaddedOutput;
    private double[][] invOutput;

    STFT() {
    }

    double[][] getInvPaddedOutput() {
        return this.invPaddedOutput;
    }

    int getPaddedLength() {
        return paddedLength;
    }

    double[][][] getReSTFT() {
        return this.reOutput;
    }

    double[][][] getImSTFT() {
        return this.imOutput;
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

    void stftm(double[][] input, int winLen, int nOverlap, String winFunc) {
        Log.i("DEBUG", "STFT running");

        this.input = input;

        this.nChannels = input.length;
        this.nSamples = input[0].length;  // time is in the second dimension

        this.winLen = winLen;
        this.nOverlap = nOverlap;
        this.hopSize = winLen - nOverlap;

        this.nFrames = (int) ceil(nSamples / nOverlap) + (int) ceil(nOverlap / hopSize);
        this.paddedLength = nFrames * hopSize + 2 * nOverlap;

        this.paddedInput = new double[nChannels][paddedLength];

        this.nFreq = winLen / 2 + 1;

        this.reOutput = new double[nChannels][nFrames][this.nFreq];
        this.imOutput = new double[nChannels][nFrames][this.nFreq];

        // padding
        for (int c = 0; c < nChannels; c++) {
            //Log.i("DEBUG", "input: " + Arrays.toString(this.input[c]));
            System.arraycopy(this.input[c], 0, paddedInput[c], nOverlap, nSamples);
            //Log.i("DEBUG", "paddedInput: " + Arrays.toString(paddedInput[c]));
        }


        this.window = new double[winLen];
        getSTFTwindow(winFunc);

        //Log.i("DEBUG", "window: " + Arrays.toString(window));

        setScalingVec(true);    // true for forward fft
        for (int c = 0; c < nChannels; c++) {
            for (int i = 0; i < nFrames; i++) {
                double[] frame = new double[winLen];
                double[] wFrame = new double[winLen];

                System.arraycopy(paddedInput[c], i * hopSize, frame, 0, winLen);

                for (int j = 0; j < winLen; j++) {
                    wFrame[j] = frame[j] * window[j] * iswin[i * hopSize + j];
                    //Log.i("DEBUG","frame[j]" + frame[j] + ", window[j]" + window[j] + ", iswin" + iswin[i * hopSize + j] + ", iswin * window[j]" + window[j]*iswin[i * hopSize + j]);
                }

                double[] reInputFrame = new double[winLen];
                double[] imInputFrame = new double[winLen];

                System.arraycopy(wFrame, 0, reInputFrame, 0, winLen);

                int result = fft(reInputFrame, imInputFrame);

                if (result < 0) {
                    throw new RuntimeException("STFT failed. Real input and imaginary input are of different lengths.");
                }

                System.arraycopy(reInputFrame, 0, this.reOutput[c][i], 0, this.nFreq);
                System.arraycopy(imInputFrame, 0, this.imOutput[c][i], 0, this.nFreq);
            }
        }
    }

    void istftm(double[][][] reSTFT, double[][][] imSTFT, int winLen, int nOverlap, String winFunc, int nSamples) {

        this.reOutput = reSTFT;
        this.imOutput = imSTFT;

        this.nChannels = reSTFT.length;
        int inv_nFrames = reSTFT[0].length;
        //Log.i("DEBUG", "nFrames = " + inv_nFrames);

        this.winLen = winLen;
        this.nOverlap = nOverlap;
        this.hopSize = winLen - nOverlap;

        this.nSamples = nSamples;  // time is in the second dimension

        //Log.i("DEBUG", "reSTFT[0].length = " + reSTFT[0].length);

        this.paddedLength = inv_nFrames * hopSize + 2 * nOverlap;

        //Log.i("DEBUG", "Initializing window");

        this.window = new double[winLen];
        this.invPaddedOutput = new double[nChannels][paddedLength];

        int nFreq = winLen / 2 + 1;

        //Log.i("DEBUG", "Preparing inverse output");
        this.input = new double[nChannels][nSamples];

        //Log.i("DEBUG", "Preparing window");
        getSTFTwindow(winFunc);

        //Log.i("DEBUG", "Preparing scaling vector");
        setScalingVec(false);   // false for inverse fft

        this.invOutput = new double[nChannels][nSamples];

        Log.i("DEBUG", "Running InvSTFT");
        for (int c = 0; c < nChannels; c++) {
            Log.i("DEBUG", "Channel " + c);
            for (int i = 0; i < inv_nFrames; i++) {
                Log.i("DEBUG", "Frame " + i + " out of " + inv_nFrames);
                double[] reInputFrame = new double[winLen];
                double[] imInputFrame = new double[winLen];

                //Log.i("DEBUG", "Copying real STFT to frame");
                System.arraycopy(reSTFT[c][i], 0, reInputFrame, 0, nFreq);

                //Log.i("DEBUG", "Reconstructing real");

                for (int j = 0; j < winLen / 2 - 1; j++) {
                    reInputFrame[nFreq + j] = reSTFT[c][i][nFreq - 2 - j];
                }

                //Log.i("DEBUG", "Copying imag STFT to frame");
                System.arraycopy(imSTFT[c][i], 0, imInputFrame, 0, nFreq);

                //Log.i("DEBUG", "Reconstructing imag");
                for (int j = 0; j < winLen / 2 - 1; j++) {
                    imInputFrame[nFreq + j] = -imSTFT[c][i][nFreq - 2 - j];
                }

                //Log.i("DEBUG", "Displaying real");
                //Log.i("DEBUG", Arrays.toString(reInputFrame));

                //Log.i("DEBUG", "Displaying imag");
                //Log.i("DEBUG", Arrays.toString(imInputFrame));

                //Log.i("DEBUG", "Running iFFT");

                int result = ifft(reInputFrame, imInputFrame);

                //Log.i("DEBUG", "reInputFrame size = " + reInputFrame.length);
                //Log.i("DEBUG", "reInputFrame: " + Arrays.toString(reInputFrame));
                //Log.i("DEBUG", "imInputFrame: " + Arrays.toString(imInputFrame));

                if (result < 0) {
                    throw new RuntimeException("Inverse STFT failed.");
                }

                //Log.i("DEBUG", "Copying inv result to frame");
                double[] fframe = new double[winLen];
                System.arraycopy(reInputFrame, 0, fframe, 0, winLen);

                double[] wFrame = new double[winLen];
                //Log.i("DEBUG", "Reconstructing signal");

                for (int j = 0; j < winLen; j++) {
                    wFrame[j] = window[j] * fframe[j] * iswin[i * hopSize + j];
                    invPaddedOutput[c][i * hopSize + j] = invPaddedOutput[c][i * hopSize + j] + wFrame[j];
                }
                //Log.i("DEBUG", "invPaddedOutput: " + Arrays.toString(invPaddedOutput[c]));
            }
            //Log.i("DEBUG", "Trimming padding");

            System.arraycopy(invPaddedOutput[c], nOverlap, this.invOutput[c], 0, nSamples);

            //Log.i("DEBUG", "invOutput: " + Arrays.toString(invOutput[c]));
            Log.i("DEBUG", "istft done");
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
        Log.i("DEBUG", Arrays.toString(iswin));
    }

    private void getSTFTwindow(String winFunc) {
        if (winFunc.equals("sine")) {
            getSineWindow();
        } else if (winFunc.equals("hann")) {
            getHannWindow();
        }
    }

    private void getHannWindow() {
        for (int i = 0; i < winLen; i++) {
            window[i] = 0.5 * (1 - cos(2 * PI * i / (winLen - 1)));
        }
    }

    private void getSineWindow() {
        for (int i = 0; i < winLen; i++) {
            window[i] = sin(PI / winLen * ((double) i + 0.5));
        }
    }

    int fft(double x[], double y[]) {

        /*--------------------------------------------------------*/
        /* Modified by Karn Watcharasupat                         */
        /*--------------------------------------------------------*/
        /* fft.c                                                  */
        /* (c) Douglas L. Jones                                   */
        /* University of Illinois at Urbana-Champaign             */
        /* January 19, 1992                                       */
        /*                                                        */
        /*   fft: in-place radix-2 DIT DFT of a complex input     */
        /*                                                        */
        /*   input:                                               */
        /* n: length of FFT: must be a power of two               */
        /* m: n = 2**m                                            */
        /*   input/output                                         */
        /* x: double array of length n with real part of data     */
        /* y: double array of length n with imag part of data     */
        /*                                                        */
        /*   Permission to copy and use this program is granted   */
        /*   under a Creative Commons "Attribution" license       */
        /*   http://creativecommons.org/licenses/by/1.0/          */
        /*--------------------------------------------------------*/

        Log.i("DEBUG", "FFT running");

        int nFrames;

        if (x.length != y.length) {
            return -1;
        } else {
            nFrames = x.length;
        }

        int m = (int) (log(nFrames) / log(2));

        int i, j, k, n1, n2;
        double c, s, e, a, t1, t2;

        j = 0; /* bit-reverse */
        n2 = nFrames / 2;
        for (i = 1; i < nFrames - 1; i++) {
            n1 = n2;
            while (j >= n1) {
                j = j - n1;
                n1 = n1 / 2;
            }
            j = j + n1;

            if (i < j) {
                t1 = x[i];
                x[i] = x[j];
                x[j] = t1;
                t1 = y[i];
                y[i] = y[j];
                y[j] = t1;
            }
        }

        /* FFT */
        n1 = 0;
        n2 = 1;

        for (i = 0; i < m; i++) {
            n1 = n2;
            n2 = n2 + n2;
            e = -2 * PI / n2;
            a = 0.0;

            for (j = 0; j < n1; j++) {
                c = cos(a);
                s = sin(a);
                a = a + e;

                for (k = j; k < nFrames; k = k + n2) {
                    t1 = c * x[k + n1] - s * y[k + n1];
                    t2 = s * x[k + n1] + c * y[k + n1];
                    x[k + n1] = x[k] - t1;
                    y[k + n1] = y[k] - t2;
                    x[k] = x[k] + t1;
                    y[k] = y[k] + t2;
                }
            }
        }
        return 1;
    }

    int ifft(double[] reInputFrame, double[] imInputFrame) {
        //Inverse FFT using conjugation trick
        //Log.i("DEBUG", "IFFT running");

        if (reInputFrame.length != imInputFrame.length) {
            //Log.i("DEBUG", "Real and imag are of unequal lengths!");
            return -1;
        } else {
            inv_nFreq = reInputFrame.length;
        }

        // conjugating the input
        //Log.i("DEBUG", "Conjugating");
        for (int i = 0; i < inv_nFreq; i++) {
            imInputFrame[i] = -imInputFrame[i];
        }

        // fft
        //Log.i("DEBUG", "FFT on conjugate");
        int fftResult = fft(reInputFrame, imInputFrame);

        if (fftResult < 0) {
            return -2;
        }

        // conjugating the output
        //Log.i("DEBUG", "Conjugate again");
        for (int i = 0; i < inv_nFreq; i++) {
            imInputFrame[i] = -imInputFrame[i] / inv_nFreq;
            reInputFrame[i] = reInputFrame[i] / inv_nFreq;
        }

        return 1;

    }
}
