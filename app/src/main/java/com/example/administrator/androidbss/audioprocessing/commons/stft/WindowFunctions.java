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

import static java.lang.Math.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;

/**
 * A class for window function calculation and default overlap
 *
 * @author Karn Watcharasupat
 * School of Electrical and Electronic Engineering,
 * Nanyang Technological University, Singapore
 * @version April 2019
 */
public class WindowFunctions {
    public static final int SINE_WINDOW = 1000;
    public static final int PERIODIC_HAMMING_WINDOW = 1010;

    /**
     * Calculate a window of specified length
     *
     * @param winFunc the window function
     * @param winLen  window length
     * @return array containing values of the window function
     */
    static double[] calculate(int winFunc, int winLen) {

        double[] window = new double[winLen];

        switch (winFunc) {
            case SINE_WINDOW:
                for (int i = 0; i < winLen; i++) {
                    window[i] = sin(PI / winLen * (double) i);
                }
                break;
            case PERIODIC_HAMMING_WINDOW:
                for (int i = 0; i < winLen; i++) {
                    window[i] = 0.54 - 0.46 * cos(2 * PI * i / winLen);
                }
                break;
        }

        return window;
    }

    /**
     * Get constant overlap-add overlap as a fraction of the window length
     *
     * @param winFunc the window function
     * @return constant overlap-add overlap as a fraction of the window length
     */
    static double getDefaultOverlap(int winFunc) {
        switch (winFunc) {
            case SINE_WINDOW:
                return 0.50;
            case PERIODIC_HAMMING_WINDOW:
                return 0.25;
            default:
                return 0;
        }
    }
}
