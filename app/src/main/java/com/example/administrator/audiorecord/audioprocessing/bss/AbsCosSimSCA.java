package com.example.administrator.audiorecord.audioprocessing.bss;

import android.util.Log;

import com.example.administrator.audiorecord.audioprocessing.commons.ComplexSingularValueDecomposition;
import com.example.administrator.audiorecord.audioprocessing.commons.WhiteningNew;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import java.util.Arrays;
import java.util.stream.IntStream;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.random;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;

public class AbsCosSimSCA {

    private Complex[][][] STFTin, STFTout;

    private int nSrc, nItr;
    private int nChannels, nFrames, nFreqs;
    private boolean isProjBack, derivCheck;

    private final double EPSILON = 1e-15;
    private final double PROGRESS_TOLERANCE = 1e-9;

    Array2DRowFieldMatrix<Complex>[] STFTwhite;
    Array2DRowFieldMatrix<Complex>[] Y;

    Array2DRowFieldMatrix<Complex>[] demix;

    Complex eta;

    //Complex[][] randomInit;

    public AbsCosSimSCA(Complex[][][] STFTin, int nItr, boolean isProjBack, int nSrc, double eta, boolean derivCheck) {

        this.STFTin = STFTin;

        this.nItr = nItr;

        this.eta = new Complex(eta);

        this.derivCheck = derivCheck;

        this.nSrc = nSrc;

        this.nChannels = STFTin.length;
        this.nFrames = STFTin[0].length;
        this.nFreqs = STFTin[0][0].length;

        STFTout = new Complex[this.nSrc][nFrames][nFreqs];

        STFTwhite = new Array2DRowFieldMatrix[nFreqs];
        Y = new Array2DRowFieldMatrix[nFreqs];
        demix = new Array2DRowFieldMatrix[nFreqs];

        /*
        randomInit = new Complex[nChannels][this.nSrc];
        */

        /*
        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                if (c == s) {
                    randomInit[c][s] = Complex.ONE;
                } else {
                    randomInit[c][s] = Complex.ZERO;
                }
            }
        }
        */
    }

    public void run() {
        whiten();

        IntStream.range(0, nFreqs).forEach(this::runThisFreq); //parallel()

        Log.i("DEBUG", "COMPLETED ALL FREQUENCIES");
    }

    private void whiten() {

        /*
        IMPORTANT:

        STFTin is nChannels by nFrames by nFreq to match STFT's output

        However,

        STFTout is nFreq by nChannels by nFrames to simplify further calculations
        */

        WhiteningNew whitening = new WhiteningNew(STFTin);
        whitening.run();
        STFTwhite = whitening.getWhitenedMatrixArray();
    }

    private void runThisFreq(int f) {

        Log.i("DEBUG", "bin = " + f);

        /* normalizing the columns of X */

        Array2DRowFieldMatrix<Complex> X_bar = STFTwhite[f];

        double colNorm;

        for (int t = 0; t < nFrames; t++) {

            colNorm = 0.0;

            for (int c = 0; c < nChannels; c++) {
                colNorm += pow(X_bar.getEntry(c, t).abs(), 2.0);
            }

            colNorm = sqrt(colNorm);

            for (int c = 0; c < nChannels; c++) {
                X_bar.multiplyEntry(c, t, Complex.ONE.divide(colNorm));
            }
        }

        /* checking */
        /*
            double checkColNorm;

            for (int t = 0; t < nFrames; t++) {

                checkColNorm = 0.0;

                for (int c = 0; c < nChannels; c++) {
                    checkColNorm += pow(X_bar.getEntry(c, t).abs(), 2.0);
                }

                colNorm = sqrt(checkColNorm);

                Log.i("DEBUG", "frame " + t + " norm = " + checkColNorm);
            }
        */

        /* Initializing mixing matrix */

        Complex[][] randomInit = new Complex[nChannels][nSrc];

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                randomInit[c][s] = new Complex(random(), random());
            }
        }


        /* decouple the rows of A */

        Array2DRowFieldMatrix<Complex> A = new Array2DRowFieldMatrix<>(randomInit);

        ComplexSingularValueDecomposition csvd = new ComplexSingularValueDecomposition(A, true);

        A = csvd.getU().multiply(csvd.getVH());

        PowerMeanContrast contrast = new PowerMeanContrast(-0.5, true);

        if (derivCheck) {
            contrast.computeCost(A, X_bar);

            Array2DRowFieldMatrix<Complex> derivCheckGrad = contrast.getGrad();

            Log.i("DEBUG", "grad = " + derivCheckGrad.toString());

            Complex[][] randDirection = new Complex[nChannels][nSrc];

            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    randDirection[c][s] = new Complex(signum(random() - 0.5));
                }
            }

            Log.i("DEBUG", "rand = " + Arrays.deepToString(randDirection));

            double[] gradVec = new double[2 * nChannels * nSrc];
            double[] directVec = new double[2 * nChannels * nSrc];
            Complex tempGrad;

            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    tempGrad = derivCheckGrad.getEntry(c, s);

                    gradVec[s * nChannels + c] = tempGrad.getReal();
                    gradVec[nChannels * nSrc + s * nChannels + c] = tempGrad.getImaginary();

                    directVec[s * nChannels + c] = randDirection[c][s].getReal();
                    directVec[nChannels * nSrc + s * nChannels + c] = randDirection[c][s].getImaginary();
                }
            }

            double gtd1 = 0.0;

            for (int i = 0; i < 2 * nChannels * nSrc; i++) {
                gtd1 += gradVec[i] * directVec[i];
            }

            double doubleEps = 1e-15;
            Complex eps = new Complex(doubleEps);

            Array2DRowFieldMatrix<Complex> epsRandDirection = (Array2DRowFieldMatrix<Complex>) new Array2DRowFieldMatrix<>(randDirection).scalarMultiply(eps);

            contrast.computeCost(A.add(epsRandDirection), X_bar);
            double costPlus = contrast.getCost();
            Log.i("DEBUG", "costPlus = " + costPlus);

            contrast.computeCost(A.subtract(epsRandDirection), X_bar);
            double costMinus = contrast.getCost();
            Log.i("DEBUG", "costMinus = " + costMinus);

            double gtd2 = abs(costPlus - costMinus) / (2 * doubleEps);

            double relDiff = abs(gtd1 - gtd2) / (abs(gtd1) + abs(gtd2));

            Log.i("DEBUG", "analytical deriv = " + gtd1);
            Log.i("DEBUG", "numerical deriv = " + gtd2);
            Log.i("DEBUG", "Relative difference between analytical and numerical directional-derivative is " + relDiff);
        }

        /* Nesterov's accelerated gradient */

        Array2DRowFieldMatrix<Complex> A_nag = A;
        Array2DRowFieldMatrix<Complex> A_prev;

        int counter = 0;
        Complex gamma;

        double cost = 0.0, prevCost = 0.0;
        Array2DRowFieldMatrix<Complex> grad;

        for (int iter = 0; iter < nItr; iter++) {
            A_prev = A;

            if (iter > 0) {
                prevCost = cost;
            }

            gamma = new Complex((counter - 1.0) / (counter + 2.0));

            contrast.computeCost(A, X_bar);

            cost = contrast.getCost();
            grad = contrast.getGrad();

            A = A_nag.subtract((Array2DRowFieldMatrix<Complex>) grad.scalarMultiply(eta));

            A_nag = (Array2DRowFieldMatrix<Complex>) A.scalarMultiply(Complex.ONE.add(gamma)).subtract(A_prev.scalarMultiply(gamma));

            counter++;

            if ((iter > 0) &&
                    abs(cost - prevCost) < PROGRESS_TOLERANCE){
                break;
            }

            if ((iter > 0) && (prevCost < cost)) {
                counter = 0;
                A_nag = A;
            }
        }

        csvd = new ComplexSingularValueDecomposition(A, true);

        A = csvd.getU().multiply(csvd.getVH());

        demix[f] = conjugate((Array2DRowFieldMatrix<Complex>) A.transpose());
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

    public class PowerMeanContrast {

        double r, inv_r;

        double g = 0.0;
        double[][] dg = new double[nSrc][nFrames];

        double cost, prevCost;

        boolean isRowDecoupling;

        ComplexSingularValueDecomposition csvd;

        Array2DRowFieldMatrix<Complex> grad;

        Array2DRowFieldMatrix<Complex> Sinv, C, Scov;

        PowerMeanContrast(double r, boolean isRowDecoupling) {
            this.r = r;
            this.inv_r = 1.0 / r;
            this.isRowDecoupling = isRowDecoupling;
        }

        void computeCost(Array2DRowFieldMatrix<Complex> A, Array2DRowFieldMatrix<Complex> X) {

            prevCost = cost;

            Complex[][] Aorth = A.getData();

            if (isRowDecoupling) {
                csvd = new ComplexSingularValueDecomposition(A, true);

                Aorth = csvd.getU().multiply(csvd.getVH()).getData();
            }

            double[] colNorm = new double[nSrc];
            Complex[] invColNorm = new Complex[nSrc];

            for (int s = 0; s < nSrc; s++) {

                colNorm[s] = 0.0;

                for (int c = 0; c < nChannels; c++) {
                    colNorm[s] += pow(A.getEntry(c, s).abs(), 2.0);
                }

                colNorm[s] = sqrt(colNorm[s]);

                invColNorm[s] = Complex.ONE.divide(colNorm[s]);

                for (int c = 0; c < nChannels; c++) {
                    A.multiplyEntry(c, s, invColNorm[s]);
                }
            }

            //angle: nSrc by nFrames
            Array2DRowFieldMatrix<Complex> angle = (Array2DRowFieldMatrix<Complex>) conjugate(A).transpose().multiply(X);

            Array2DRowRealMatrix Y2_contra_mat = new Array2DRowRealMatrix(nSrc, nFrames);
            double[][] Y2_contra; // = new double[nFrames][nSrc]
            double[] Y2_contra_min = new double[nFrames];
            double[] Y2_power_mean_contra = new double[nFrames];

            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    Y2_contra_mat.setEntry(s, t, 1.0 - pow(angle.getEntry(s, t).abs(), 2.0));
                }
            }

            Y2_contra = Y2_contra_mat.transpose().getData();

            double colMin;
            for (int t = 0; t < nFrames; t++) {

                colMin = Arrays.stream(Y2_contra_mat.getColumn(t)).min().orElse(EPSILON);

                if (colMin == 0) {
                    colMin = EPSILON;
                }

                Y2_contra_min[t] = colMin;
            }

            double frameSum;

            g = 0.0;

            for (int t = 0; t < nFrames; t++) {

                frameSum = 0.0;

                for (int s = 0; s < nSrc; s++) {
                    frameSum += pow(Y2_contra[t][s] / Y2_contra_min[t], r);
                }

                frameSum = pow(frameSum, inv_r) / pow(nSrc, inv_r);

                Y2_power_mean_contra[t] = frameSum * Y2_contra_min[t];

                g += Y2_power_mean_contra[t];
            }

            for (int t = 0; t < nFrames; t++) {
                for (int s = 0; s < nSrc; s++) {
                    dg[s][t] = -pow(Y2_power_mean_contra[t] / Y2_contra[t][s], 1.0 - r) / nSrc;
                }
            }

            /* scaling */

            double scale = 1.0;
            cost = (scale / nFrames) * g;

            Complex gradMultiplier = new Complex(2.0).multiply(scale / nFrames);

            Array2DRowFieldMatrix<Complex> angle_dg = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nSrc, nFrames);
            Array2DRowFieldMatrix<Complex> temp_grad;

            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    angle_dg.setEntry(s, t, angle.getEntry(s, t).multiply(dg[s][t]));
                }
            }

            temp_grad = (Array2DRowFieldMatrix<Complex>) X.multiply((Array2DRowFieldMatrix<Complex>) angle_dg.transpose()).scalarMultiply(gradMultiplier);

            /* Apply the chain rule through column norm normalization */

            Complex[] Aorth_multiplier = new Complex[nSrc];
            Complex[][] grad_array = new Complex[nChannels][nSrc];

            for (int s = 0; s < nSrc; s++) {

                Aorth_multiplier[s] = Complex.ZERO;

                for (int c = 0; c < nChannels; c++) {
                    Aorth_multiplier[s] = Aorth_multiplier[s].add(Aorth[c][s].conjugate().multiply(temp_grad.getEntry(c, s)));
                }

                Aorth_multiplier[s].divide(pow(colNorm[s], 3.0));
            }

            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    grad_array[c][s] = temp_grad.getEntry(c, s).divide(colNorm[s]).subtract(Aorth[c][s].multiply(Aorth_multiplier[s]));
                }
            }

            grad = new Array2DRowFieldMatrix<>(grad_array);

            if (isRowDecoupling) {

                Sinv = csvd.getPoweredS(-1.0);

                Scov = (Array2DRowFieldMatrix<Complex>) csvd.getColumnVectorS().multiply(csvd.getColumnVectorS().transpose());

                C = Sinv.multiply(csvd.getUH()).multiply(grad).multiply(csvd.getV());

                for (int ci = 0; ci < nChannels; ci++) {
                    for (int cj = 0; cj < nChannels; cj++) {
                        C.multiplyEntry(ci, cj, Complex.ONE.divide(Scov.getEntry(ci, cj).negate()));
                    }
                }

                grad = (Array2DRowFieldMatrix<Complex>) csvd.getU()
                        .multiply(conjugate(C).transpose().add(C))
                        .multiply(csvd.getS())
                        .multiply(csvd.getVH())
                        .add(csvd.getU()
                                .multiply(Sinv)
                                .multiply(csvd.getUH())
                                .multiply(grad)
                        );
            }
        }

        Array2DRowFieldMatrix<Complex> getGrad() {
            return grad;
        }

        double getCost() {
            return cost;
        }
    }

    private void demix() {
        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            Y[f] = demix[f].multiply(STFTwhite[f]);
        });
    }

    private void OutMatrixToOutArray() {
        for (int f = 0; f < nFreqs; f++) {
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    STFTout[s][t][f] = Y[f].getEntry(s, t);
                }
            }
        }
    }

    public Complex[][][] getSourceEstimatesSTFT() {

        demix();
        OutMatrixToOutArray();

        return STFTout;
    }


}
