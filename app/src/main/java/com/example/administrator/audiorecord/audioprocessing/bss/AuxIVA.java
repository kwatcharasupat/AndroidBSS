package com.example.administrator.audiorecord.audioprocessing.bss;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.FieldDecompositionSolver;
import org.apache.commons.math3.linear.FieldLUDecomposition;

import java.util.Arrays;
import java.util.stream.IntStream;

import static java.lang.Math.cosh;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.Math.tanh;

public class AuxIVA {

    private double EPSILON = 1e-15;

    private Complex[][][] STFTin, STFTout;
    //private Complex[][][] W;
    private double[][] r;
    //private Complex[][][][] V;
    //private Complex[][] Identity;

    private int nSrc, nItr;
    private int nChannels, nFrames, nFreqs;
    private boolean isProjBack;
    private String cFunc;
    private double[] cFuncParam;

    private Array2DRowFieldMatrix<Complex>[][] V; // freq by src by (chan x chan)
    private Array2DRowFieldMatrix<Complex>[] W; // freq by (chan by src)
    private Array2DRowFieldMatrix<Complex>[] Wconj; // freq by (chan by src)
    private Array2DRowFieldMatrix<Complex>[] X; // frame by freq by (frm x chan)
    private Array2DRowFieldMatrix<Complex>[] Xcopy; // frame by freq by (frm x chan)

    private Array2DRowFieldMatrix<Complex>[] Xconj;
    private Array2DRowFieldMatrix<Complex>[] Y; // frame by freq by (frm x chan)

    private Complex[][] G_r; // src by (frm by frm)


    private Array2DRowFieldMatrix<Complex> w, WV;
    private Complex norm, inorm, inFrames;
    double rSum;
    private FieldDecompositionSolver<Complex> solver;

    /*
    NOTE

    W = (w1* w2* ... wk*)^T
     */

    private Array2DRowFieldMatrix Identity;

    public AuxIVA(Complex[][][] STFTin,
                  int nItr, boolean isProjBack,
                  Complex[][][] W0,
                  String cFunc, double[] cFuncParam) {
        /*
        Implementation of AuxIVA algorithm for BSS presented in
        N. Ono, *Stable and fast update rules for independent vector analysis based
        on auxiliary function technique*, Proc. IEEE, WASPAA, 2011.

        Reference
        ---------
        Robin Scheibler (2018), Blind Source Separation using Independent Vector Analysis with Auxiliary Function

        Parameters
        ----------
        STFTin: Complex][][] (nChannels by nFrames by nFreqs)

        nItr: int
            number of iterations

        isProjBack: boolean
            scaling on the first microphone by back projection

        W0: Complex[][][] (nFreqs by nChannels by nSrc )
            initial mixing matrix, optional

        cFunc:  String
            contrast function

        cFunc: double[]
            contrast function parameters = {C, m}

        Output
        ------
        STFTout: Complex[][][] (nSrc by nFrames by nFreqs)
            STFT representation of separated signal
        */

        Log.i("DEBUG", "AuxIVA now preparing");

        //Log.i("DEBUG", "Getting STFT dimension - channels");
        this.nChannels = STFTin.length;
        //Log.i("DEBUG", "Getting STFT dimension - frames");
        this.nFrames = STFTin[0].length;
        this.inFrames = Complex.ONE.divide(nFrames);
        //Log.i("DEBUG", "Getting STFT dimension - freqs");
        this.nFreqs = STFTin[0][0].length;
        //Log.i("DEBUG", "Getting STFT dimension - done");

        this.nItr = nItr;
        this.isProjBack = isProjBack;
        this.nSrc = nChannels;  // force equal number of sources and sensors
        // using nSrc as a separate variable to allow future addition of the overdetermined case

        //Log.i("DEBUG", "Initializing STFTin");
        this.STFTin = new Complex[this.nChannels][this.nFrames][this.nFreqs];

        // preemptively prevents the original STFT data from being modified
        //Log.i("DEBUG", "Creating a duplicate of STFT data");
        for (int c = 0; c < this.nChannels; c++) {
            for (int i = 0; i < this.nFrames; i++) {
                System.arraycopy(STFTin[c][i], 0, this.STFTin[c][i], 0, this.nFreqs);
            }
        }

        Complex[][] temp = new Complex[nFrames][nChannels];
        Complex[][] tempConj = new Complex[nFrames][nChannels];

        X = new Array2DRowFieldMatrix[nFreqs];
        Xcopy = new Array2DRowFieldMatrix[nFreqs];
        Xconj = new Array2DRowFieldMatrix[nFreqs];
        Y = new Array2DRowFieldMatrix[nFreqs];
        for (int f = 0; f < this.nFreqs; f++) {
            for (int t = 0; t < this.nFrames; t++) {
                for (int c = 0; c < this.nChannels; c++) {
                    temp[t][c] = STFTin[c][t][f];
                    tempConj[t][c] = STFTin[c][t][f].conjugate();
                    //Log.i("DEBUG", "STFTin[c][t][f] = " + STFTin[c][t][f] + ", temp[c] = " + temp[c]);
                }
            }
            X[f] = new Array2DRowFieldMatrix<>(temp);
            Xconj[f] = new Array2DRowFieldMatrix<>(tempConj);
            Xcopy[f] = (Array2DRowFieldMatrix<Complex>) X[f].copy();
            Y[f] = (Array2DRowFieldMatrix<Complex>) X[f].copy();
        }

        Identity = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nSrc);
        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                if (s == c) {
                    Identity.setEntry(c, s, Complex.ONE);
                } else {
                    Identity.setEntry(c, s, Complex.ZERO);
                }
            }
        }

        // initializing the demixing matrix
        //Log.i("DEBUG", "Initializing demixing matrix");
        W = new Array2DRowFieldMatrix[nFreqs];
        Wconj = new Array2DRowFieldMatrix[nFreqs];
        if (W0 == null) {
            //Log.i("DEBUG", "Initializing with identity");

            for (int f = 0; f < this.nFreqs; f++) {
                W[f] = (Array2DRowFieldMatrix<Complex>) Identity.copy();
                Wconj[f] = (Array2DRowFieldMatrix<Complex>) Identity.copy();
            }
            //Log.i("DEBUG", "Initializing with identity - done");
        } else {
            //Log.i("DEBUG", "Initializing with specified value");
            for (int f = 0; f < nFreqs; f++) {
                W[f] = new Array2DRowFieldMatrix<>(W0[f]);
            }
        }

        // initializing the separated source matrix
        this.STFTout = new Complex[this.nSrc][this.nFrames][this.nFreqs];

        // initializing contrast function and related variables
        this.cFunc = cFunc;
        this.cFuncParam = new double[cFuncParam.length];
        System.arraycopy(cFuncParam, 0, this.cFuncParam, 0, cFuncParam.length);

        this.r = new double[this.nFrames][this.nSrc];
        this.G_r = new Complex[this.nFrames][this.nSrc];

        //PRE-ALLOCATING MEMORY---------------------------------------------------------------------

        V = new Array2DRowFieldMatrix[nFreqs][nSrc];
        for (int f = 0; f < nFreqs; f++) {
            for (int s = 0; s < nSrc; s++) {
                V[f][s] = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nChannels);
            }
        }

        for (int t = 0; t < nFrames; t++) {
            for (int s = 0; s < nSrc; s++) {
                G_r[t][s] = Complex.ZERO;
            }
        }
        //------------------------------------------------------------------------------------------
    }

    public void run() {
        Log.i("DEBUG", "AuxIVA now running");

        for (int runCount = 0; runCount < nItr; runCount++) {
            Log.i("DEBUG", "Iteration: " + runCount);
            //demixing
            demix(Y, X, W);
            calculateG();
            calculateV();
            updateDemix();
        }
        demix(Y, X, W);
        OutMatrixToOutArray();

        if (isProjBack) {
            projectBack();
        }
        Log.i("DEBUG", "AuxIVA - done");
    }

    private void projectBack() {
        //Log.i("DEBUG", "Projecting back");
        Complex[][] scale = projectBack(STFTout, STFTin[0]);
        Complex thisSrcFrameScaleConj;

        for (int s = 0; s < nSrc; s++) {
            for (int f = 0; f < nFreqs; f++) {
                thisSrcFrameScaleConj = scale[f][s].conjugate();
                for (int t = 0; t < nFrames; t++) {
                    STFTout[s][t][f] = STFTout[s][t][f].multiply(thisSrcFrameScaleConj);
                }
            }
        }
    }

    private void updateDemix() {
        //Log.i("DEBUG", "Updating demixing matrix");
        //updating the demixing matrix
        for (int f = 0; f < nFreqs; f++) {
            //Log.i("DEBUG", "bin = " + f);
            for (int s = 0; s < nSrc; s++) {
                //Log.i("DEBUG", "src = " + s);
                WV = (Array2DRowFieldMatrix<Complex>) Wconj[f].transpose().multiply(V[f][s]);
                //Log.i("DEBUG", "WV = " + WV);

                solver = new FieldLUDecomposition<>(WV).getSolver();

                boolean isNonSingular = solver.isNonSingular();
                if (isNonSingular) {
                    w = (Array2DRowFieldMatrix<Complex>) solver.solve((Array2DRowFieldMatrix<Complex>) (Identity.getColumnMatrix(s)));
                    //Log.i("DEBUG", "WV = " + WV);
                    //Log.i("DEBUG", "I("+s+") = " + Identity.getColumnMatrix(s));
                    //Log.i("DEBUG", "w("+s+") = " + w);
                    //Log.i("DEBUG", "V[" + f + "][" + s + "] = " + V[f][s]);
                    norm = (conjugate(w).transpose().multiply(V[f][s].multiply(w))).getEntry(0, 0).sqrt();
                    //Log.i("DEBUG", "norm = " + norm);
                    inorm = Complex.ONE.divide(norm);
                    //Log.i("DEBUG", "inorm = " + inorm);
                    w = (Array2DRowFieldMatrix<Complex>) w.scalarMultiply(inorm);
                    //Log.i("DEBUG", "w("+s+") = " + w);
                    W[f].setColumnMatrix(s, w);
                    //Log.i("DEBUG", "W["+f+"]("+s+") = " + W[f]);
                } else {
                    Log.i("DEBUG", "WV is singular.");
                    break;
                }
            }
        }
        //Log.i("DEBUG", "Updating demixing matrix - done");
    }

    private void calculateV() {
        //Log.i("DEBUG", "Calculating V");
        // calculate weighted covariance matrix V

        for (int f = 0; f < nFreqs; f++) {
            //Log.i("DEBUG", "bin = " + f);
            for (int s = 0; s < nSrc; s++) {
                // Log.i("DEBUG", "src = " + s);
                for (int t = 0; t < nFrames; t++) {
                    Xcopy[f].multiplyEntry(t, s, G_r[t][s]);
                }
                V[f][s] = (Array2DRowFieldMatrix<Complex>) (Xcopy[f].transpose().multiply(Xconj[f])).scalarMultiply(inFrames);
                //Log.i("DEBUG", "V["+f+"]["+s+"] = " + V[f][s] + ", nFrames = " + nFrames);
            }
            Xcopy[f] = (Array2DRowFieldMatrix<Complex>) X[f].copy();    //resetting Xcopy values to original
        }
        //Log.i("DEBUG", "Calculating V - done");
    }

    private void calculateG() {
        //Log.i("DEBUG", "Calculating r and G/r");
        // calculate r and G/r
        // r(s) = sqrt(sum(|y(f)|^2))

        for (int t = 0; t < nFrames; t++) {
            for (int s = 0; s < nSrc; s++) {
                rSum = 0;
                for (int f = 0; f < nFreqs; f++) {
                    rSum += pow(Y[f].getEntry(0, s).abs(), 2.0);
                }
                r[t][s] = sqrt(rSum);    // nSrc by nFrames
                //Log.i("DEBUG", "r["+t+"]["+s+"] = " + r[t][s]);
                if (r[t][s] == 0) {
                    Log.i("DEBUG", "r is zero!");
                    r[t][s] = EPSILON;
                }
                G_r[t][s] = (Complex.ZERO).add(cFuncCalc(cFunc, "df", r[t][s]) / r[t][s]); // nSrc by nFrames
                //Log.i("DEBUG", "G_r["+t+"]["+s+"] = " + G_r[t][s]);
            }
        }
        //Log.i("DEBUG", "Calculating r and G/r - done");
    }

    public void runDemixTest() {

        //2.0 * (random() - 0.5)

        for (int f = 0; f < nFreqs; f++) {
            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    if (s == c) {
                        W[f].setEntry(c, s, new Complex(2.0, 0.0));
                    } else {
                        W[f].setEntry(c, s, new Complex(0.0, 0.0));
                    }
                }
            }
            Log.i("DEBUG", "W[" + f + "] = " + W[f]);
        }


        demix(Y, X, W);

        OutMatrixToOutArray();
        for (int s = 0; s < nSrc; s++) {
            for (int t = 0; t < nFrames; t++) {
                Log.i("DEBUG", "STFTin[" + s + "][" + t + "] : " + Arrays.deepToString(STFTin[s][t]));
                Log.i("DEBUG", "STFTout[" + s + "][" + t + "] : " + Arrays.deepToString(STFTout[s][t]));
            }
        }
    }

    public Complex[][][] getSourceEstimatesSTFT() {
        return STFTout;
    }


    private void demix(Array2DRowFieldMatrix<Complex>[] Y, Array2DRowFieldMatrix<Complex>[] X, Array2DRowFieldMatrix<Complex>[] W) {
        //Log.i("DEBUG", "Demixing");


        for (int f = 0; f < nFreqs; f++) {
            Wconj[f] = conjugate(W[f]);
        }
        for (int f = 0; f < nFreqs; f++) {
            Y[f] = X[f].multiply(Wconj[f]);
        }

        //Log.i("DEBUG", "Y = " + Arrays.deepToString(Y));
        //Log.i("DEBUG", "Demixing - done");
    }

    private void OutMatrixToOutArray() {
        for (int t = 0; t < this.nFrames; t++) {
            for (int f = 0; f < this.nFreqs; f++) {
                for (int c = 0; c < this.nSrc; c++) {
                    STFTout[c][t][f] = Y[f].getEntry(t, c);
                }
            }
        }
    }

    private double cFuncCalc(String func, String deriv, double r) {

        double result = 0.0;

        double C, m = 0.0;

        C = cFuncParam[0];
        if (cFuncParam.length > 1) {
            m = cFuncParam[1];
        }

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
        //Log.i("DEBUG", "Calculating scale matrix");
        Complex[][] scale = new Complex[nFreqs][nSrc];
        Complex[][] num = new Complex[nFreqs][nSrc];
        double[][] denom = new double[nFreqs][nSrc];

        //Log.i("DEBUG", "out.length (src)= " + out.length);
        //Log.i("DEBUG", "out.length[0] (frm)= " + out[0].length);
        //Log.i("DEBUG", "out.length[0][0] (bin)= " + out[0][0].length);

        //Log.i("DEBUG", "ref.length (frm)= " + ref.length);
        //Log.i("DEBUG", "ref.length[0] (bin)= " + ref[0].length);

        for (int f = 0; f < nFreqs; f++) {
            //Log.i("DEBUG", "bin = " + f);
            for (int s = 0; s < nSrc; s++) {
                //Log.i("DEBUG", "src = " + s);
                num[f][s] = Complex.ZERO;
                scale[f][s] = Complex.ZERO;
                for (int t = 0; t < nFrames; t++) {
                    //Log.i("DEBUG", "frm = " + t);
                    //Log.i("DEBUG", "ref[t][f] = " + ref[t][f]);
                    //Log.i("DEBUG", "out[s][t][f] = " + out[s][t][f]);

                    num[f][s] = num[f][s].add(ref[t][f].conjugate().multiply(out[s][t][f]));

                    //Log.i("DEBUG", "num = " + num[f][s]);
                    denom[f][s] += (pow(out[s][t][f].abs(), 2.0));
                    //Log.i("DEBUG", "denom = " + denom[f][s]);
                }
                if (denom[f][s] > 0) {
                    //Log.i("DEBUG", "denom > 0");
                    scale[f][s] = num[f][s].divide(denom[f][s]);
                    //Log.i("DEBUG", "scale = " + scale[f][s]);
                }
            }
        }

        //Log.i("DEBUG", "Calculating scale matrix - done");

        return scale;
    }

    private Array2DRowFieldMatrix<Complex> conjugate(Array2DRowFieldMatrix<Complex> z) {

        Array2DRowFieldMatrix<Complex> zconj = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), z.getRowDimension(), z.getColumnDimension());

        for (int r = 0; r < z.getRowDimension(); r++) {
            for (int c = 0; c < z.getColumnDimension(); c++) {
                zconj.setEntry(r, c, z.getEntry(r, c).conjugate());
            }
        }

        return zconj;
    }
}
