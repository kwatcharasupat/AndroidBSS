package com.example.administrator.audiorecord.audioprocessing.commons;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.FastFourierTransformer;

import java.util.stream.IntStream;

import static java.lang.Math.PI;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static org.apache.commons.math3.transform.DftNormalization.STANDARD;
import static org.apache.commons.math3.transform.TransformType.FORWARD;
import static org.apache.commons.math3.transform.TransformType.INVERSE;

public class STFT {

    /*
    Multichannel Short Time Fourier Transform and Inverse Short Time Fourier Transform for Java
    2018 (c) Karn Watcharasupat, GNU Public License v3.0
     */

    private int nFrames;
    private int nFreqs;

    private int inv_nFrames;
    private int inv_nFreqs;

    private double[] iswin;

    private Complex[][][] output;

    private double[][] invOutput;

    public STFT() {
    }

    public void stftm(double[][] input, int winLen, int nOverlap, String winFunc) {

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
        output: Complex[][][] (nChannels by nFrames by nFreqs)
            The STFT representation of input signal
            call getSTFT() to obtain the matrix
         */

        Log.i("DEBUG", "STFT running");

        int nChannels = input.length;
        int nSamples = input[0].length;

        int hopSize = winLen - nOverlap;

        this.nFrames = (int) ceil((double) nSamples / nOverlap) + (int) ceil((double) nOverlap / hopSize);
        int paddedLength = nFrames * hopSize + 2 * nOverlap;

        double[][] paddedInput = new double[nChannels][paddedLength];

        this.nFreqs = winLen / 2 + 1;

        this.output = new Complex[nChannels][nFrames][this.nFreqs];

        // padding
        for (int c = 0; c < nChannels; c++) {
            System.arraycopy(input[c], 0, paddedInput[c], nOverlap, nSamples);
        }

        double[] window = new double[winLen];
        getSTFTwindow(winFunc, window, winLen);

        setScalingVec(true, paddedLength, winLen, hopSize, nFrames, window);    // true for forward fft
        for (int c = 0; c < nChannels; c++) {
            int finalC = c;
            IntStream.range(0, nFrames).parallel().forEach(t -> {
                double[] frame = new double[winLen];
                double[] wFrame = new double[winLen];

                System.arraycopy(paddedInput[finalC], t * hopSize, frame, 0, winLen);

                for (int j = 0; j < winLen; j++) {
                    wFrame[j] = frame[j] * window[j] * iswin[t * hopSize + j];
                }

                FastFourierTransformer FFT = new FastFourierTransformer(STANDARD);

                System.arraycopy(FFT.transform(wFrame, FORWARD), 0, output[finalC][t], 0, nFreqs);
            });
        }
        Log.i("DEBUG", "STFT completed");
    }

    public void istftm(Complex[][][] STFT, int winLen, int nOverlap, String winFunc, int nSamples) {
        /*
        Multichannel Inverse Short Time Fourier Transform

        Parameters
        ----------
        input: Complex[][][] (nChannels by nSamples by nFreqs)
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

        int inv_nChannels = STFT.length;
        int inv_nFrames = STFT[0].length;

        int inv_hopSize = winLen - nOverlap;

        int inv_paddedLength = inv_nFrames * inv_hopSize + 2 * nOverlap;

        double[] inv_window = new double[winLen];
        double[][] invPaddedOutput = new double[inv_nChannels][inv_paddedLength];

        int nFreq = winLen / 2 + 1;

        getSTFTwindow(winFunc, inv_window, winLen);

        setScalingVec(false, inv_paddedLength, winLen, inv_hopSize, inv_nFrames, inv_window);   // false for inverse fft

        this.invOutput = new double[inv_nChannels][nSamples];
        for (int c = 0; c < inv_nChannels; c++) {
            int finalC = c;
            IntStream.range(0, inv_nFrames).parallel().forEach(i -> {
                Complex[] inputFrame = new Complex[winLen];
                System.arraycopy(STFT[finalC][i], 0, inputFrame, 0, nFreq);

                for (int j = 0; j < winLen / 2 - 1; j++) {
                    inputFrame[nFreq + j] = STFT[finalC][i][nFreq - 2 - j].conjugate();
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
                    wFrame[j] = inv_window[j] * fframe[j] * iswin[i * inv_hopSize + j];
                    invPaddedOutput[finalC][i * inv_hopSize + j] = invPaddedOutput[finalC][i * inv_hopSize + j] + wFrame[j];
                }
            });
        }

        for (int c = 0; c < inv_nChannels; c++) {
            System.arraycopy(invPaddedOutput[c], nOverlap, this.invOutput[c], 0, nSamples);
        }

        Log.i("DEBUG", "Inverse STFT completed");
    }

    private void setScalingVec(boolean isFwd, int paddedLength, int winLen, int hopSize,
                               int nFrames, double[] window) {

        final double EPSILON = 1e-8;
        double[] swin = new double[paddedLength];
        iswin = new double[paddedLength];

        IntStream.range(0, nFrames).parallel().forEach(i -> {
            for (int j = 0; j < winLen; j++) {
                swin[i * hopSize + j] = swin[i * hopSize + j] + window[j] * window[j];
            }
        });

        IntStream.range(0, paddedLength).parallel().forEach(i -> {
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
        });
    }

    private void getSTFTwindow(String winFunc, double[] window, int winLen) {
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

    public Complex[][][] getSTFT() {
        return output;
    }

    public double[][] getRealSigFromInvSTFT() {
        return this.invOutput;
    }

    public int get_nFrames() {
        return nFrames;
    }

    public int get_nFreqs() {
        return nFreqs;
    }
}
