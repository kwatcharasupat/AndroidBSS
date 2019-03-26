package com.example.administrator.audiorecord.audioprocessing.commons;

import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import java.util.ArrayList;

import static com.example.administrator.audiorecord.audioprocessing.commons.Conjugate.transjugate;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.sqrt;

public class K_Hyperlines {

    private Array2DRowFieldMatrix<Complex> A;

    private final double EPSILON = 1e-8;
    private final double PROGRESS_TOLERANCE = 1e-5;

    public K_Hyperlines(Array2DRowFieldMatrix<Complex> X_bar, int nSrc, int f) {

        boolean isConverged = false;

        int maxEpoch = 100;

        int nChannels = X_bar.getRowDimension();
        int nFrames = X_bar.getColumnDimension();

        Complex[][] randomInit = new Complex[nChannels][nSrc];

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                randomInit[c][s] = new Complex(random(), random());
            }
        }

        A = new Array2DRowFieldMatrix<>(randomInit);

        double[] colNorm = new double[nSrc];

        for (int s = 0; s < nSrc; s++) {

            colNorm[s] = 0.0;

            for (int c = 0; c < nChannels; c++) {
                colNorm[s] += pow(A.getEntry(c, s).abs(), 2.0);
            }

            colNorm[s] = sqrt(colNorm[s] + EPSILON);
        }

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                A.multiplyEntry(c, s, Complex.ONE.divide(colNorm[s]));
            }
        }

        Array2DRowFieldMatrix<Complex> angle;
        Array2DRowRealMatrix D = new Array2DRowRealMatrix(nSrc, nFrames);

        double[] Dmin = new double[nFrames];
        int[] DminIdx = new int[nFrames];
        int vvSize;

        double cost, prevCost = 0.0;
        double scale = 1.0;

        ArrayList<Complex>[] vvTemp = new ArrayList[nChannels];

        for (int c = 0; c < nChannels; c++) {
            vvTemp[c] = new ArrayList<>();
        }

        Array2DRowFieldMatrix<Complex> vv, Rj, eRj;

        for (int epoch = 0; epoch < maxEpoch; epoch++) {

            angle = transjugate(A).multiply(X_bar);

            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    D.setEntry(s, t, 1.0 - pow(angle.getEntry(s, t).abs(), 2.0));
                }
            }

            cost = 0.0;

            for (int t = 0; t < nFrames; t++) {

                DminIdx[t] = 0;

                for (int s = 0; s < nSrc; s++) {
                    if (D.getEntry(s, t) < D.getEntry(DminIdx[t], t)) {
                        DminIdx[t] = s;
                        Dmin[t] = D.getEntry(s, t);
                    }
                }

                cost += Dmin[t];
            }

            cost = cost * scale / nFrames;

            if ((epoch > 0) && (abs(cost - prevCost) < PROGRESS_TOLERANCE)) {

                Log.i("DEBUG", "Bin " + f + " K-hyperlines completed at epoch " + epoch);

                isConverged = true;

                break;
            }

            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    if (s == DminIdx[t]) {
                        for (int c = 0; c < nChannels; c++) {
                            vvTemp[c].add(X_bar.getEntry(c, t));
                        }
                    }
                }

                vvSize = vvTemp[0].size();

                vv = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, vvSize);

                for (int c = 0; c < nChannels; c++) {
                    for (int t = 0; t < vvSize; t++) {
                        vv.setEntry(c, t, vvTemp[c].get(t));
                    }
                }

                Rj = vv.multiply(transjugate(vv));
                eRj = new ComplexSingularValueDecomposition(Rj, true).getU();
                A.setColumnMatrix(s, eRj.getColumnMatrix(0));
            }

            prevCost = cost;
        }

        if (!isConverged) {
            Log.i("DEBUG", "Bin " + f + " K-hyperlines completed after max epoch ");
        }
    }

    public Array2DRowFieldMatrix<Complex> getA() {
        return A;
    }
}
