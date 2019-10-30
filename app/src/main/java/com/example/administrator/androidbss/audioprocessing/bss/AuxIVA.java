package com.example.administrator.androidbss.audioprocessing.bss;

import android.util.Log;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.FieldDecompositionSolver;
import org.apache.commons.math3.linear.FieldLUDecomposition;

import java.util.Arrays;
import java.util.stream.IntStream;

import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.conjugate;
import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.transjugate;
import static com.example.administrator.androidbss.audioprocessing.commons.Epsilon.EPSILON;
import static org.apache.commons.math3.util.FastMath.cosh;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.tanh;

public class AuxIVA {

    private Complex[][][] STFTin, STFTout;
    private double[][] r;

    private int nSrc, nItr;
    private int nChannels, nFrames, nFreqs;
    private boolean isProjBack;
    private String cFunc;
    private double[] cFuncParam;

    private Array2DRowFieldMatrix<Complex>[][] V; // freq by src by (chan x chan)
    private Array2DRowFieldMatrix<Complex>[] W; // freq by (chan by src)
    private Array2DRowFieldMatrix<Complex>[] X; // freq by (frm x chan)

    private Array2DRowFieldMatrix<Complex>[] Xconj;
    private Array2DRowFieldMatrix<Complex>[] Y; // freq by (frm x chan)

    private Complex[][] G_r; // src by (frm by frm)

    private Complex[] norm, inorm;
    private Array2DRowFieldMatrix<Complex>[] w, WV;
    private FieldDecompositionSolver<Complex>[] solver;

    private Complex inFrames;

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

        //Log.i("DEBUG", "AuxIVA now preparing");

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
        //this.STFTin = new Complex[this.nChannels][this.nFrames][this.nFreqs];

        // preemptively prevents the original STFT data from being modified
        //Log.i("DEBUG", "Creating a duplicate of STFT data");

        this.STFTin = SerializationUtils.clone(STFTin);

        X = new Array2DRowFieldMatrix[nFreqs];
        Xconj = new Array2DRowFieldMatrix[nFreqs];
        Y = new Array2DRowFieldMatrix[nFreqs];
        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            Complex[][] temp = new Complex[nFrames][nChannels];
            for (int t = 0; t < this.nFrames; t++) {
                for (int c = 0; c < this.nChannels; c++) {
                    temp[t][c] = STFTin[c][t][f];
                }
            }
            X[f] = new Array2DRowFieldMatrix<>(temp, false);
            Y[f] = (Array2DRowFieldMatrix<Complex>) X[f].copy();
            Xconj[f] = transjugate(X[f]);
        });


        /*
        System.arraycopy(X, 0, Xcopy, 0, nFreqs);
        System.arraycopy(X, 0, STFTper, 0, nFreqs);
        */

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
        if (W0 == null) {
            //Log.i("DEBUG", "Initializing with identity");
            IntStream.range(0, nFreqs).parallel().forEach(f -> {
                W[f] = (Array2DRowFieldMatrix<Complex>) Identity.copy();
            });
            //Log.i("DEBUG", "Initializing with identity - done");
        } else {
            //Log.i("DEBUG", "Initializing with specified value");
            IntStream.range(0, nFreqs).parallel().forEach(f -> {
                W[f] = new Array2DRowFieldMatrix<>(W0[f], false);
            });
        }

        // initializing contrast function and related variables
        this.cFunc = cFunc;
        this.cFuncParam = new double[cFuncParam.length];
        System.arraycopy(cFuncParam, 0, this.cFuncParam, 0, cFuncParam.length);

        this.r = new double[this.nFrames][this.nSrc];
        this.G_r = new Complex[this.nFrames][this.nSrc];

        //PRE-ALLOCATING MEMORY---------------------------------------------------------------------

        V = new Array2DRowFieldMatrix[nFreqs][nSrc];

        for (int f = 0; f < nFreqs; f++) {
            Arrays.fill(V[f], new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nChannels));
        }

        for (int t = 0; t < nFrames; t++) {
            Arrays.fill(G_r[t], Complex.ZERO);
        }

        norm = new Complex[nFreqs];
        inorm = new Complex[nFreqs];
        w = new Array2DRowFieldMatrix[nFreqs];
        WV = new Array2DRowFieldMatrix[nFreqs];
        solver = new FieldDecompositionSolver[nFreqs];


        //------------------------------------------------------------------------------------------
    }

    public void run() {
        //Log.i("DEBUG", "AuxIVA now running");

        for (int epoch = 0; epoch < nItr; epoch++) {
            Log.i("DEBUG", "iter = " + epoch);
            demix(Y, X, W);
            Log.i("DEBUG", "G");
            calculateG();
            Log.i("DEBUG", "V");
            calculateV();
            Log.i("DEBUG", "update");
            updateDemix();
        }

        demix(Y, X, W);
        outMatrixToOutArray();
        if (isProjBack) {
            projectBack();
        }
        //Log.i("DEBUG", "AuxIVA - done");
    }

    private void calculateG() {
        Complex[][][] Yarray = new Complex[nFreqs][][];
        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            Yarray[f] = Y[f].getDataRef();
        });

        IntStream.range(0, nFrames).parallel().forEach(t -> {
            for (int s = 0; s < nSrc; s++) {
                r[t][s] = 0.0;
                for (int f = 0; f < nFreqs; f++) {
                    r[t][s] += Yarray[f][t][s].multiply(Yarray[f][t][s].conjugate()).getReal();
                }
                r[t][s] = sqrt(r[t][s]); // nFrames by nSrc

                if (r[t][s] == 0) {
                    r[t][s] = EPSILON;
                }

                G_r[t][s] = new Complex(cFuncCalc(cFunc, "df", r[t][s]) / r[t][s]); // nFrames by nSrc
            }
        });
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

    private void calculateV() {

        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            Complex[][] tempX = X[f].getData();
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    for (int c = 0; c < nChannels; c++) {
                        tempX[t][c] = tempX[t][c].multiply(G_r[t][s]);
                    }
                }

                V[f][s] = (Array2DRowFieldMatrix<Complex>) (Xconj[f].multiply(new Array2DRowFieldMatrix<>(tempX, false))).transpose();//.scalarMultiply(inFrames);
            }
        });
    }

    private void updateDemix() {
        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            for (int s = 0; s < nSrc; s++) {
                WV[f] = transjugate(W[f]).multiply(V[f][s]);
                solver[f] = new FieldLUDecomposition<>(WV[f]).getSolver();

                if (solver[f].isNonSingular()) {
                    w[f] = (Array2DRowFieldMatrix<Complex>) solver[f].solve((Array2DRowFieldMatrix<Complex>) (Identity.getColumnMatrix(s)));
                    norm[f] = (transjugate(w[f]).multiply(V[f][s].multiply(w[f]))).getEntry(0, 0).sqrt();
                    inorm[f] = Complex.ONE.divide(norm[f]);
                    w[f] = (Array2DRowFieldMatrix<Complex>) w[f].scalarMultiply(inorm[f]);
                    W[f].setColumnMatrix(s, w[f]);
                } else {
                    Log.i("DEBUG", "WV is singular.");
                }
            }
        });
    }

    private void demix(Array2DRowFieldMatrix<Complex>[] Y, Array2DRowFieldMatrix<Complex>[] X, Array2DRowFieldMatrix<Complex>[] W) {
        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            Y[f] = X[f].multiply(conjugate(W[f]));
        });
    }

    private void projectBack() {
        Complex[][] scale = projBackScale(STFTout, STFTin[0]);

        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            for (int s = 0; s < nSrc; s++) {
                Complex thisSrcFrameScaleConj;
                thisSrcFrameScaleConj = scale[s][f].conjugate();
                for (int t = 0; t < nFrames; t++) {
                    STFTout[s][t][f] = STFTout[s][t][f].multiply(thisSrcFrameScaleConj);
                }
            }
        });
    }

    private Complex[][] projBackScale(Complex[][][] out, Complex[][] ref) {
        Complex[][] scale = new Complex[nSrc][nFreqs];
        Complex[][] num = new Complex[nSrc][nFreqs];
        double[][] denom = new double[nSrc][nFreqs];

        for (int s = 0; s < nSrc; s++) {
            int finalS = s;
            IntStream.range(0, nFreqs).parallel().forEach(f -> {
                num[finalS][f] = Complex.ZERO;
                scale[finalS][f] = Complex.ZERO;
                for (int t = 0; t < nFrames; t++) {
                    num[finalS][f] = num[finalS][f].add(ref[t][f].conjugate().multiply(out[finalS][t][f]));
                    denom[finalS][f] += out[finalS][t][f].multiply(out[finalS][t][f].conjugate()).getReal();
                }
                if (denom[finalS][f] > 0) {
                    scale[finalS][f] = num[finalS][f].divide(denom[finalS][f]);
                }
            });
        }
        return scale;
    }

    private void outMatrixToOutArray() {

        this.STFTout = new Complex[this.nSrc][this.nFrames][this.nFreqs];

        IntStream.range(0, nFreqs).parallel().forEach(f -> {

            Complex[][] Yarray = Y[f].getDataRef();

            for (int t = 0; t < nFrames; t++) {
                for (int c = 0; c < this.nSrc; c++) {
                    STFTout[c][t][f] = Yarray[t][c];
                }
            }
        });
    }

    public Complex[][][] getSourceEstimatesSTFT() {
        return STFTout;
    }

    public void printDemix() {
        for (int f = 0; f < nFreqs; f++) {
            Log.i("DEBUG", "W[" + f + "]= " + W[f].toString());
        }
    }
}

