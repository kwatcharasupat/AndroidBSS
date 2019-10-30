package com.example.administrator.androidbss.audioprocessing.bss;

import android.util.Log;

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

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.conjugate;
import static com.example.administrator.androidbss.audioprocessing.commons.Conjugate.transjugate;
import static com.example.administrator.androidbss.audioprocessing.commons.Epsilon.EPSILON;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;


public class FastICA {

    public static final int POWER_CONTRAST = 1004;
    public static final int SQRT_CONTRAST = 1001;
    public static final int LOG_CONTRAST = 1002;
    public static final int KURTOSIS_CONTRAST = 1003;

    private Complex[][][] input;
    private Complex[][][][] STFTout;

    private int nSrc = 2, maxItr;
    private int nChannels, nFrames, nFreqs;

    private double PROGRESS_TOLERANCE = 1e-9;

    private Array2DRowFieldMatrix<Complex>[] STFTwhite;

    private double beta = 12.5;

    private final double POWER_R_DEFAULT = 1.25;

    private double contrastFunctionParams[];

    private int contrastFunction;

    private double[][][] posteriorProb;

    private double r = 1.25;

    //Complex[][] randomInit;


    private FastICA() {

    }

    public static FastICA initiate() {
        return new FastICA();
    }

    public FastICA setMaximumIterations(int maxItr) {
        this.maxItr = maxItr;
        return this;
    }

    public FastICA setMaskSoftness(double beta) {
        this.beta = beta;
        return this;
    }

    public FastICA setProgressTolerance(double progTol) {
        this.PROGRESS_TOLERANCE = progTol;
        return this;
    }

    public FastICA setContrastFunction(int contrastFunction) {
        this.contrastFunction = contrastFunction;
        return this;
    }

    public FastICA setContrasFunctionParameters(double... parameters) {
        switch (contrastFunction) {
            case POWER_CONTRAST:
            case SQRT_CONTRAST:
            case LOG_CONTRAST:
                contrastFunctionParams = new double[1];
                this.contrastFunctionParams[0] = parameters[0];
                break;
            case KURTOSIS_CONTRAST:
                contrastFunctionParams = null;
                break;
        }

        return this;
    }

    public FastICA setInputData(Complex[][][] input) {
        this.nChannels = input.length;
        this.nFrames = input[0].length;
        this.nFreqs = input[0][0].length;

        this.input = new Complex[nChannels][nFrames][nFreqs];

        for (int c = 0; c < nChannels; c++) {
            for (int t = 0; t < nFrames; t++) {
                System.arraycopy(input[c][t], 0, this.input[c][t], 0, nFreqs);
            }
        }

        return this;
    }

    public FastICA otherwiseUseDefault() {
        return this;
    }

    public Complex[][][][] run() {
        whiten();

        posteriorProb = new double[nFreqs][nFrames][nSrc];
        IntStream.range(0, nFreqs).parallel().forEach(this::runThisFreq);
        permutationAlignment();

        return STFTout;
    }

    private void permutationAlignment() {

        double[][][] P = PermutationAlignment.initiate()
                .setMaximumIterations(10)
                .setProgressTolerance(1e-3)
                .setEnabledFinetune(true)
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
                        STFTout[finalS][finalC][finalT][f] = input[finalC][finalT][f].multiply(P[f][finalT][finalS]);
                    });
                }
            }
        }
    }

    private void whiten() {

        /*
        IMPORTANT:

        input is nChannels by nFrames by nFreq to match STFT's output

        However,

        STFTout is nFreq by nChannels by nFrames to simplify further calculations
        */
        STFTwhite = new Array2DRowFieldMatrix[nFreqs];
        Whitening whitening = new Whitening(input);
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

        ComplexSingularValueDecomposition csvd;

        for (int epoch = 0; epoch < maxItr; epoch++) {

            oldA = (Array2DRowFieldMatrix<Complex>) A.copy();

            y = transjugate(A).multiply(X_bar);

            contrastFunction(y,G,dG,d2G);

            a1 = X_bar.scalarMultiply(Complex.ONE.negate());
            a2 = elementwiseMultiply(conjugate(G), dG).transpose();
            a3 = cov.multiply(A).multiply(diagRowSumAbsSquare(dG));
            a4 = pseudoCov.multiply(conjugate(A)).multiply(diagRowSum(elementwiseMultiply(conjugate(G), d2G)));

            A = (Array2DRowFieldMatrix<Complex>) a1.multiply(a2).add(a3).add(a4);

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

        Complex[][] tempA = A.getDataRef();

        double out = 0.0;

        if (isMinusIdentity) {
            for (int r = 0; r < row; r++) {
                for (int c = 0; c < col; c++) {
                    if (r == c) {
                        out += pow(tempA[r][c].abs() - 1.0, 2);
                    } else {
                        out += pow(tempA[r][c].abs(), 2);
                    }
                }
            }
        } else {
            for (int r = 0; r < row; r++) {
                for (int c = 0; c < col; c++) {
                    out += pow(tempA[r][c].abs(), 2);
                }
            }
        }

        return out;
    }

    /**
     * Calculate the contrast function and its derivatives in-place
     *
     * @param y   the demixed data
     * @param G   the contrast function matrix
     * @param dG  the first derivative matrix
     * @param d2G the second derivative matrix
     */
    public void contrastFunction(Array2DRowFieldMatrix<Complex> y, Array2DRowFieldMatrix<Complex> G, Array2DRowFieldMatrix<Complex> dG, Array2DRowFieldMatrix<Complex> d2G) {

        /* Modify underlying array directly to reduce overhead*/
        Complex[][] yArray = y.getDataRef();
        Complex[][] Garray = G.getDataRef();
        Complex[][] dGarray = dG.getDataRef();
        Complex[][] d2Garray = d2G.getDataRef();

        Complex temp;

        switch (contrastFunction) {
            case POWER_CONTRAST:
                temp = null;
                for (int c = 0; c < nChannels; c++) {
                    for (int t = 0; t < nFrames; t++) {
                        Garray[c][t] = yArray[c][t].pow(contrastFunctionParams[0]);
                        dGarray[c][t] = Garray[c][t].multiply(contrastFunctionParams[0]).divide(yArray[c][t]);
                        d2Garray[c][t] = dGarray[c][t].multiply(contrastFunctionParams[0] - 1.0).divide(yArray[c][t]);
                    }
                }
                break;
            case SQRT_CONTRAST:
                for (int c = 0; c < nChannels; c++) {
                    for (int t = 0; t < nFrames; t++) {
                        temp = yArray[c][t].add(contrastFunctionParams[0]);
                        Garray[c][t] = temp.sqrt();
                        dGarray[c][t] = Garray[c][t].reciprocal().multiply(0.5);
                        d2Garray[c][t] = dGarray[c][t].divide(temp).multiply(-0.5);
                    }
                }
                break;
            case LOG_CONTRAST:
                for (int c = 0; c < nChannels; c++) {
                    for (int t = 0; t < nFrames; t++) {
                        temp = yArray[c][t].add(contrastFunctionParams[0]);
                        Garray[c][t] = temp.log();
                        dGarray[c][t] = temp.reciprocal();
                        d2Garray[c][t] = d2Garray[c][t].divide(temp).multiply(-1.0);
                    }
                }
                break;
            case KURTOSIS_CONTRAST:
                temp = null;
                for (int c = 0; c < nChannels; c++) {
                    for (int t = 0; t < nFrames; t++) {
                        Garray[c][t] = yArray[c][t].multiply(yArray[c][t]).multiply(0.5);
                        dGarray[c][t] = yArray[c][t];
                    }
                }

                for (int c = 0; c < nChannels; c++) {
                    Arrays.fill(d2Garray[c], Complex.ONE);
                }
                break;
        }
    }
}
