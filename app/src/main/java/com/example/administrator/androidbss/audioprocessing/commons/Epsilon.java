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

package com.example.administrator.androidbss.audioprocessing.commons;

import org.apache.commons.math3.complex.Complex;

/**
 * A small class for containing the value of 1e-8 in double and Complex data types
 * to prevent division by zero
 *
 * @author Karn Watcharasupat
 * School of Electrical and Electronic Engineering,
 * Nanyang Technological University, Singapore
 * @version April 2019
 */
public class Epsilon {
    public static final double EPSILON = 1e-8;
    public static final Complex COMPLEX_EPSILON = new Complex(EPSILON);
}
