package com.example.administrator.androidbss.audioprocessing.commons;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;

public class Conjugate {
    public static Array2DRowFieldMatrix<Complex> conjugate(Array2DRowFieldMatrix<Complex> z) {

        int row = z.getRowDimension();
        int col = z.getColumnDimension();

        Array2DRowFieldMatrix<Complex> zconj = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), row, col);

        for (int r = 0; r < row; r++) {
            for (int c = 0; c < col; c++) {
                zconj.setEntry(r, c, z.getEntry(r, c).conjugate());
            }
        }

        return zconj;
    }

    public static Array2DRowFieldMatrix<Complex> transjugate(Array2DRowFieldMatrix<Complex> z){

        return (Array2DRowFieldMatrix<Complex>) conjugate(z).transpose();
    }
}
