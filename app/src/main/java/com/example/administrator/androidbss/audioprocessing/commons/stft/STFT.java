/**
 * Android Blind Source Separation Project
 * Copyright (C) 2019 Karn Watcharasupat
 * <p>
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package com.example.administrator.androidbss.audioprocessing.commons.stft;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.FastFourierTransformer;

import java.util.stream.IntStream;

import static com.example.administrator.androidbss.audioprocessing.commons.Epsilon.EPSILON;
import static org.apache.commons.math3.transform.DftNormalization.STANDARD;
import static org.apache.commons.math3.transform.TransformType.FORWARD;
import static org.apache.commons.math3.transform.TransformType.INVERSE;
import static org.apache.commons.math3.util.FastMath.ceil;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Implementation of short-time Fourier Transform and its inverse for multichannel signal
 *
 * @author Karn Watcharasupat
 * School of Electrical and Electronic Engineering,
 * Nanyang Technological University, Singapore
 * @version April 2019
 */
public class STFT {
    public static final boolean FORWARD_SCALE = true;
    public static final boolean INVERSE_SCALE = false;

    /**
     * Class constructor
     */
    public STFT() {
    }

    /**
     * Perform forward STFT
     *
     * @param input        multichannel time-domain input, nChannels by nSamples
     * @param stftSettings settings for STFT via {@link STFTsettings}
     * @return STFT of input
     */
    public static Complex[][][] forwardSTFT(double[][] input, STFTsettings stftSettings) {
        int winLen = stftSettings.winLen;
        int nOverlap = stftSettings.nOverlap;
        int hopSize = winLen - nOverlap;

        int nFreqs = winLen / 2 + 1;

        int winFunc = stftSettings.winFunc;

        int nChannels = input.length;
        int nSamples = input[0].length;
        int nFrames = (int) ceil((double) (nSamples + nOverlap) / hopSize);
        int paddedLength = nFrames * hopSize + 2 * nOverlap;

        double[][] paddedInput = new double[nChannels][paddedLength];
        for (int c = 0; c < nChannels; c++) {
            System.arraycopy(input[c], 0, paddedInput[c], nOverlap, nSamples);
        }

        double[] window = WindowFunctions.calculate(winFunc, winLen);
        double[] swin = calculateScalingVector(FORWARD_SCALE, paddedLength, winLen, hopSize, nFrames, window);    // true for forward fft

        Complex[][][] output = new Complex[nChannels][nFrames][nFreqs];
        for (int c = 0; c < nChannels; c++) {
            int finalC = c;
            IntStream.range(0, nFrames).parallel().forEach(t -> {
                double[] frame = new double[winLen];
                double[] wFrame = new double[winLen];

                System.arraycopy(paddedInput[finalC], t * hopSize, frame, 0, winLen);

                for (int j = 0; j < winLen; j++) {
                    wFrame[j] = frame[j] * window[j] * swin[t * hopSize + j];
                }

                FastFourierTransformer FFT = new FastFourierTransformer(STANDARD);
                System.arraycopy(FFT.transform(wFrame, FORWARD), 0, output[finalC][t], 0, nFreqs);
            });
        }

        return output;
    }

    /**
     * Perform inverse STFT
     *
     * @param input        multichannel STFT input, nChannels by nFrames by nFreqs
     * @param nSamples     number of samples in the original time domain data,
     *                     if zero, nSamples is set to {@code paddedLength}
     * @param stftSettings settings for STFT via {@link STFTsettings}
     * @return inverse STFT of input
     */
    public static double[][] inverseSTFT(Complex[][][] input, int nSamples, STFTsettings stftSettings) {
        int winLen = stftSettings.winLen;
        int nOverlap = stftSettings.nOverlap;
        int hopSize = winLen - nOverlap;

        int nFreq = winLen / 2 + 1;

        int winFunc = stftSettings.winFunc;

        int nChannels = input.length;
        int nFrames = input[0].length;
        int paddedLength = nFrames * hopSize + 2 * nOverlap;

        if (nSamples == 0) {
            nSamples = paddedLength;
        }

        double[] window = WindowFunctions.calculate(winFunc, winLen);
        double[] swin = calculateScalingVector(INVERSE_SCALE, paddedLength, winLen, hopSize, nFrames, window);   // false for inverse fft

        double[][] paddedOutput = new double[nChannels][paddedLength];
        for (int c = 0; c < nChannels; c++) {
            int finalC = c;
            IntStream.range(0, nFrames).parallel().forEach(i -> {
                Complex[] inputFrame = new Complex[winLen];
                System.arraycopy(input[finalC][i], 0, inputFrame, 0, nFreq);

                for (int j = 0; j < winLen / 2 - 1; j++) {
                    inputFrame[nFreq + j] = input[finalC][i][nFreq - 2 - j].conjugate();
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
                    wFrame[j] = window[j] * fframe[j] * swin[i * hopSize + j];
                    paddedOutput[finalC][i * hopSize + j] = paddedOutput[finalC][i * hopSize + j] + wFrame[j];
                }
            });
        }

        double[][] output = new double[nChannels][nSamples];
        for (int c = 0; c < nChannels; c++) {
            System.arraycopy(paddedOutput[c], nOverlap, output[c], 0, nSamples);
        }

        return output;
    }

    /**
     * Calculate scaling vector for perfect reconstruction
     *
     * @param isForward    whether the scaling vector is for forward or inverse transform
     * @param paddedLength number of samples inclusive of zero padding
     * @param winLen       window length
     * @param hopSize      hop size
     * @param nFrames      number of frames
     * @param window       array containing the values of the window function
     * @return scaling vector
     */
    private static double[] calculateScalingVector(boolean isForward, int paddedLength, int winLen, int hopSize,
                                                   int nFrames, double[] window) {
        double[] swin = new double[paddedLength];

        for (int i = 0; i < nFrames; i++) {
            for (int j = 0; j < winLen; j++) {
                swin[i * hopSize + j] = swin[i * hopSize + j] + window[j] * window[j];
            }
        }

        if (isForward) {
            for (int i = 0; i < paddedLength; i++) {
                swin[i] = 1.0 / (sqrt(swin[i] * winLen) + EPSILON);
            }
        } else {
            for (int i = 0; i < paddedLength; i++) {
                swin[i] = 1.0 / (sqrt(swin[i] / winLen) + EPSILON);
            }
        }

        return swin;
    }
}
