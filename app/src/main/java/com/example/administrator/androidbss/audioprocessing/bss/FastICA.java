package com.example.administrator.androidbss.audioprocessing.bss;

import com.example.administrator.androidbss.audioprocessing.commons.ComplexSingularValueDecomposition;
import com.example.administrator.androidbss.audioprocessing.commons.PermutationAlignment;
import com.example.administrator.androidbss.audioprocessing.commons.Whitening;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Max;

import java.util.Random;
import java.util.stream.IntStream;

import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.conjugate;
import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.transjugate;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;


public class FastICA {

    private Complex[][][] STFTin;
    private Complex[][][][] STFTout;

    private int nSrc, maxItr;
    private int nChannels, nFrames, nFreqs;

    private final double EPSILON = 1e-8;

    private final double PROGRESS_TOLERANCE = 1e-9;

    private Array2DRowFieldMatrix<Complex>[] STFTwhite;
    private Array2DRowFieldMatrix<Complex>[] STFTper;

    private final double beta = 12.5;

    private final double r = 1.25;

    private double[][][] posteriorProb;

    //Complex[][] randomInit;

    public FastICA(Complex[][][] STFTin, int maxItr) {

        this.STFTin = SerializationUtils.clone(STFTin);

        this.maxItr = maxItr;

        this.nChannels = STFTin.length;
        this.nFrames = STFTin[0].length;
        this.nFreqs = STFTin[0][0].length;

        this.nSrc = nChannels;

        STFTwhite = new Array2DRowFieldMatrix[nFreqs];
        STFTper = new Array2DRowFieldMatrix[nFreqs];

        posteriorProb = new double[nFreqs][nFrames][nSrc];
    }

    public void run() {
        whiten();

        IntStream.range(0, nFreqs).parallel().forEach(this::runThisFreq); //

        //Log.i("DEBUG", "Completed all frequencies");

        permutationAlignment();
    }

    private void permutationAlignment() {

        double[][][] P = PermutationAlignment.initiate()
                .setMaximumIterations(10)
                .setProgressTolerance(1e-3)
                .setFinetuningEnabled(true)
                .setInputData(posteriorProb)
                .run();

        //Log.i("DEBUG", "Permutation aligned.");

        STFTout = new Complex[this.nSrc][nChannels][nFrames][nFreqs];

        for (int s = 0; s < nSrc; s++) {
            for (int c = 0; c < nChannels; c++) {
                for (int t = 0; t < nFrames; t++) {
                    int finalS = s;
                    int finalC = c;
                    int finalT = t;
                    IntStream.range(0, nFreqs).parallel().forEach(f -> {
                        STFTout[finalS][finalC][finalT][f] = STFTin[finalC][finalT][f].multiply(P[f][finalT][finalS]);
                    });
                }
            }
        }
    }

    private void whiten() {

        /*
        IMPORTANT:

        STFTin is nChannels by nFrames by nFreq to match STFT's output

        However,

        STFTout is nFreq by nChannels by nFrames to simplify further calculations
        */

        Whitening whitening = new Whitening(STFTin);
        whitening.run();
        STFTwhite = whitening.getWhitenedMatrixArray(); // nFreqs x (nChannels x nFrames)
    }

    private Array2DRowFieldMatrix<Complex> normalizeX(Array2DRowFieldMatrix<Complex> X) {

        Array2DRowFieldMatrix<Complex> X_bar = (Array2DRowFieldMatrix<Complex>) X.copy();

        double colNorm;

        for (int t = 0; t < nFrames; t++) {
            colNorm = 0.0;

            for (int c = 0; c < nChannels; c++) {
                colNorm += pow(X_bar.getEntry(c, t).abs(), 2);
            }

            colNorm = sqrt(colNorm + EPSILON);

            for (int c = 0; c < nChannels; c++) {
                X_bar.multiplyEntry(c, t, Complex.ONE.divide(colNorm));
            }
        }

        return X_bar;
    }

    private Array2DRowFieldMatrix<Complex> randomInit() {

        Complex[][] randomInit = new Complex[nChannels][nSrc];

        Random random = new Random();

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                randomInit[c][s] = new Complex(random.nextGaussian(), random.nextGaussian());
            }
        }

        return new Array2DRowFieldMatrix<>(randomInit);
    }

    private void runThisFreq(int f) {
        Array2DRowFieldMatrix<Complex> X_bar = normalizeX(STFTwhite[f]);
        Array2DRowFieldMatrix<Complex> A = randomInit();

        Array2DRowFieldMatrix<Complex> cov = (Array2DRowFieldMatrix<Complex>) X_bar.multiply(transjugate(X_bar)).scalarMultiply(Complex.ONE.divide(nFrames));
        Array2DRowFieldMatrix<Complex> pseudoCov = (Array2DRowFieldMatrix<Complex>) X_bar.multiply(X_bar.transpose()).scalarMultiply(Complex.ONE.divide(nFrames));

        Array2DRowFieldMatrix<Complex> oldA, y;

        Array2DRowFieldMatrix<Complex> G = new Array2DRowFieldMatrix<Complex>(ComplexField.getInstance(), nChannels, nFrames);
        Array2DRowFieldMatrix<Complex> dG = new Array2DRowFieldMatrix<Complex>(ComplexField.getInstance(), nChannels, nFrames);
        Array2DRowFieldMatrix<Complex> d2G = new Array2DRowFieldMatrix<Complex>(ComplexField.getInstance(), nChannels, nFrames);

        FieldMatrix<Complex> a1, a2, a3, a4;

        Complex thisY;

        ComplexSingularValueDecomposition csvd;

        for (int epoch = 0; epoch < maxItr; epoch++) {

            oldA = (Array2DRowFieldMatrix<Complex>) A.copy();

            y = transjugate(A).multiply(X_bar);

            for (int c = 0; c < nChannels; c++) {
                for (int t = 0; t < nFrames; t++) {

                    thisY = y.getEntry(c, t);

                    G.setEntry(c, t, thisY.pow(r));
                    dG.setEntry(c, t, G.getEntry(c, t).divide(thisY).multiply(r));
                    d2G.setEntry(c, t, dG.getEntry(c, t).divide(thisY).multiply(r - 1));
                }
            }

            a1 = X_bar.scalarMultiply(Complex.ONE.negate());
            a2 = elementwiseMultiply(conjugate(G), dG).transpose();
            a3 = cov.multiply(A).multiply(diagRowSumAbsSquare(dG));
            a4 = pseudoCov.multiply(conjugate(A)).multiply(diagRowSum(elementwiseMultiply(conjugate(G), d2G)));

            A =  (Array2DRowFieldMatrix<Complex>) a1.multiply(a2).add(a3).add(a4);

            csvd = new ComplexSingularValueDecomposition(A, true);
            A = csvd.getU().multiply(csvd.getVH());

            if (frobeniusNorm(transjugate(oldA).multiply(A), true) < PROGRESS_TOLERANCE) {
                break;
            }
        }

        Array2DRowFieldMatrix<Complex> finalAngle = transjugate(A).multiply(X_bar); // (nSrc by nChannels) * (nChannels by nFrames)

        posteriorProb[f] = calculateMask(finalAngle).getData();
    }

    private Array2DRowRealMatrix calculateMask(Array2DRowFieldMatrix<Complex> finalAngle) {
        Array2DRowRealMatrix cos2H = new Array2DRowRealMatrix(nSrc, nFrames);

        for (int s = 0; s < nSrc; s++) {
            for (int t = 0; t < nFrames; t++) {
                cos2H.setEntry(s, t, pow(finalAngle.getEntry(s, t).abs(), 2));
            }
        }

        double[] maxCos = new double[nFrames];
        Max max = new Max();
        for (int t = 0; t < nFrames; t++) {
            maxCos[t] = max.evaluate(cos2H.getColumn(t));
        }

        Array2DRowRealMatrix mask = new Array2DRowRealMatrix(nFrames, nSrc);
        double[] maskColSum = new double[nFrames];

        for (int t = 0; t < nFrames; t++) {
            maskColSum[t] = 0.0;
            for (int s = 0; s < nSrc; s++) {
                mask.setEntry(t, s, exp(beta * (cos2H.getEntry(s, t) - maxCos[t])));
                maskColSum[t] += mask.getEntry(t, s);
            }
            for (int s = 0; s < nSrc; s++) {
                mask.multiplyEntry(t, s, 1.0 / maskColSum[t]);
            }
        }

        return mask;
    }


    public Complex[][][][] getSourceEstimatesSTFT() {
        return STFTout;
    }

    private Array2DRowFieldMatrix<Complex> elementwiseMultiply(Array2DRowFieldMatrix<Complex> A, Array2DRowFieldMatrix<Complex> B) {

        int row = A.getRowDimension();
        int col = A.getColumnDimension();

        Array2DRowFieldMatrix<Complex> out = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), row, col);

        for (int r = 0; r < row; r++) {
            for (int c = 0; c < col; c++) {
                out.setEntry(r, c, A.getEntry(r, c).multiply(B.getEntry(r, c)));
            }
        }

        return out;
    }

    private Array2DRowFieldMatrix<Complex> diagRowSum(Array2DRowFieldMatrix<Complex> A) {

        int row = A.getRowDimension();

        Array2DRowFieldMatrix<Complex> out = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), row, row);

        Complex[] rowSum = new Complex[row];

        for (int r = 0; r < row; r++) {
            rowSum[r] = Complex.ZERO;
            for (int c = 0; c < row; c++) {
                rowSum[r].add(A.getEntry(r, c));
            }
        }

        for (int r = 0; r < row; r++) {
            for (int c = 0; c < row; c++) {
                if (r == c) {
                    out.setEntry(r, c, rowSum[r]);
                } else {
                    out.setEntry(r, c, Complex.ZERO);
                }

            }
        }

        return out;
    }

    private Array2DRowFieldMatrix<Complex> diagRowSumAbsSquare(Array2DRowFieldMatrix<Complex> A) {

        int row = A.getRowDimension();

        Array2DRowFieldMatrix<Complex> out = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), row, row);

        Complex[] rowSum = new Complex[row];

        for (int r = 0; r < row; r++) {
            rowSum[r] = Complex.ZERO;
            for (int c = 0; c < row; c++) {
                rowSum[r].add(pow(A.getEntry(r, c).abs(), 2));
            }
        }

        for (int r = 0; r < row; r++) {
            for (int c = 0; c < row; c++) {
                if (r == c) {
                    out.setEntry(r, c, rowSum[r]);
                } else {
                    out.setEntry(r, c, Complex.ZERO);
                }

            }
        }

        return out;
    }

    private double frobeniusNorm(Array2DRowFieldMatrix<Complex> A, boolean isMinusIdentity) {

        int row = A.getRowDimension();
        int col = A.getColumnDimension();

        double out = 0.0;

        if (isMinusIdentity) {
            for (int r = 0; r < row; r++) {
                for (int c = 0; c < col; c++) {
                    if (r == c) {
                        out += pow((A.getEntry(r, c)).abs() - 1.0, 2);
                    } else {
                        out += pow((A.getEntry(r, c)).abs(), 2);
                    }
                }
            }
        } else {
            for (int r = 0; r < row; r++) {
                for (int c = 0; c < col; c++) {
                    out += pow(A.getEntry(r, c).abs(), 2);
                }
            }
        }

        return out;
    }
}
