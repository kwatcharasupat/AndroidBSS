package com.example.administrator.androidbss.audioprocessing.commons;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;

import VisualNumerics.math.ComplexSVD;

import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.transjugate;

public class ComplexSingularValueDecomposition {

    ComplexSVD complexSVD;

    VisualNumerics.math.Complex[][] Uarray, Varray;
    VisualNumerics.math.Complex[] Sarray;

    Array2DRowFieldMatrix<Complex> U, S, VH;

    double[] SentryDouble;
    Complex[] SentryComplex;

    boolean econ;

    int nRow, nCol, nColS, nColV;

    public ComplexSingularValueDecomposition(Array2DRowFieldMatrix<Complex> input, boolean econ) {
        this.econ = econ;

        nRow = input.getRowDimension();
        nCol = input.getColumnDimension();

        if (econ) {
            nColS = nRow;
            nColV = nRow;
        } else {
            nColS = nCol;
            nColV = nCol;
        }

        VisualNumerics.math.Complex complexArray[][] = new VisualNumerics.math.Complex[nRow][nCol];
        Complex entry;
        for (int i = 0; i < nRow; i++) {
            for (int j = 0; j < nCol; j++) {
                entry = input.getEntry(i, j);
                if (entry != null) {
                    complexArray[i][j] = new VisualNumerics.math.Complex(entry.getReal(), entry.getImaginary());
                } else {
                    Log.i("DEBUG", "NULL ENTRY!");
                    complexArray[i][j] = new VisualNumerics.math.Complex(0.0);
                }
            }
        }

        complexSVD = new ComplexSVD(complexArray);

        Uarray = complexSVD.U();
        Sarray = complexSVD.S();
        Varray = complexSVD.V();

        Complex[][] tempU = new Complex[nRow][nRow];
        Complex[][] tempS = new Complex[nRow][nColS];
        Complex[][] tempVH = new Complex[nColV][nCol];

        VisualNumerics.math.Complex tempEntry;

        for (int i = 0; i < nRow; i++) {
            for (int j = 0; j < nRow; j++) {
                tempEntry = Uarray[i][j];
                tempU[i][j] = new Complex(tempEntry.re, tempEntry.im);
            }
        }

        U = new Array2DRowFieldMatrix<>(tempU,false);

        SentryDouble = new double[nRow];
        SentryComplex = new Complex[nRow];

        for (int i = 0; i < nRow; i++) {
            SentryDouble[i] = Sarray[i].re;

            SentryComplex[i] = new Complex(SentryDouble[i]);
            tempS[i][i] = SentryComplex[i];

            for (int j = 0; j < nColS; j++) {
                if (i != j) {
                    tempS[i][j] = Complex.ZERO;
                }
            }
        }

        S = new Array2DRowFieldMatrix<>(tempS,false);

        for (int i = 0; i < nCol; i++) {
            for (int j = 0; j < nColV; j++) {

                tempEntry = Varray[i][j];

                tempVH[j][i] = new Complex(tempEntry.re, -tempEntry.im);
            }
        }

        VH = new Array2DRowFieldMatrix<>(tempVH, false);
    }

    public Array2DRowFieldMatrix<Complex> getU() {
        return U;
    }

    public Array2DRowFieldMatrix<Complex> getUH() {
        return transjugate(U);
    }

    public Array2DRowFieldMatrix<Complex> getV() {
        return transjugate(VH);
    }

    public Array2DRowFieldMatrix<Complex> getVH() {
        return VH;
    }

    public Array2DRowFieldMatrix<Complex> getS() {
        return S;
    }

    public Array2DRowFieldMatrix<Complex> getColumnVectorS() {
        return new Array2DRowFieldMatrix<>(SentryComplex);
    }

    public Complex[] getSentryComplex() {
        return SentryComplex;
    }

    public double[] getSentryDouble() {
        return SentryDouble;
    }

    public Array2DRowFieldMatrix<Complex> getPoweredS(double power) {
        Array2DRowFieldMatrix<Complex> poweredS = S;

        for (int i = 0; i < nRow; i++) {
            poweredS.setEntry(i, i, S.getEntry(i, i).pow(power));
        }

        return poweredS;
    }
}
