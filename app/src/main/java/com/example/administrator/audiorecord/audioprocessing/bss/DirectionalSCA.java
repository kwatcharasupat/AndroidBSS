package com.example.administrator.audiorecord.audioprocessing.bss;

import android.util.Log;

import com.example.administrator.audiorecord.audioprocessing.commons.ComplexSingularValueDecomposition;
import com.example.administrator.audiorecord.audioprocessing.commons.PermutationAlignment;
import com.example.administrator.audiorecord.audioprocessing.commons.Whitening;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Max;
import org.apache.commons.math3.stat.descriptive.rank.Min;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

import java.util.Arrays;
import java.util.stream.IntStream;

import static com.example.administrator.audiorecord.audioprocessing.commons.Conjugate.transjugate;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.signum;
import static org.apache.commons.math3.util.FastMath.sqrt;


public class DirectionalSCA {

    private final Complex[][][] STFTin;
    private Complex[][][] STFTout;

    private int nSrc, maxItr;
    private int nChannels, nFrames, nFreqs;
    private boolean derivCheck, isRowDecoupling;

    private final double EPSILON = 1e-8;

    private final double OPTIMUM_TOLERANCE = 1e-5;
    private final double PROGRESS_TOLERANCE = 1e-9;

    private Array2DRowFieldMatrix<Complex>[] STFTwhite;

    private Complex eta;

    private final double beta = 12.5;

    private double[][][] posteriorProb;


    public DirectionalSCA(Complex[][][] STFTin, int maxItr, int nSrc, double eta, boolean derivCheck, boolean isRowDecoupling) {

        this.STFTin = SerializationUtils.clone(STFTin);

        this.maxItr = maxItr;
        this.eta = new Complex(eta);

        this.derivCheck = derivCheck;
        this.isRowDecoupling = isRowDecoupling;

        this.nSrc = nSrc;

        this.nChannels = STFTin.length;
        this.nFrames = STFTin[0].length;
        this.nFreqs = STFTin[0][0].length;

        STFTout = new Complex[this.nSrc][nFrames][nFreqs];

        STFTwhite = new Array2DRowFieldMatrix[nFreqs];

        posteriorProb = new double[nFreqs][nFrames][nSrc];
    }

    public void run() {
        whiten();

        IntStream.range(0, nFreqs).parallel().forEach(this::runThisFreq);

        permutationAlignment();
    }


    private void permutationAlignment() {

        PermutationAlignment permutationAlignment = new PermutationAlignment(posteriorProb, 10, 1e-3, true);

        permutationAlignment.run();
        double[][][] P = permutationAlignment.getP();

        //Log.i("DEBUG", "Permutation aligned.");


        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            for (int t = 0; t < nFrames; t++) {
                for (int s = 0; s < nSrc; s++) {
                    synchronized (STFTin) {
                        STFTout[s][t][f] = STFTin[0][t][f].multiply(P[f][t][s]);
                    }
                }
            }
        });

        //Log.i("DEBUG", "SCA completed");
    }

    private void whiten() {

        /*
        IMPORTANT:
        STFTin is nChannels by nFrames by nFreq to match STFT's output
        However, STFTout is nFreq by nChannels by nFrames to simplify further calculations
        */

        Whitening whitening = new Whitening(STFTin);
        whitening.run();
        STFTwhite = whitening.getWhitenedMatrixArray(); // nFreqs x (nChannels x nFrames)
    }


    private Array2DRowFieldMatrix<Complex> normalizeX(Array2DRowFieldMatrix<Complex> X){

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

    private Array2DRowFieldMatrix<Complex> randomInit(){

        Complex[][] randomInit = new Complex[nChannels][nSrc];

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                randomInit[c][s] = new Complex(random(), random());
            }
        }

        return new Array2DRowFieldMatrix<>(randomInit);
    }

    private Array2DRowFieldMatrix<Complex> normalizeA(Array2DRowFieldMatrix<Complex> Ain){

        Array2DRowFieldMatrix<Complex> Aout = (Array2DRowFieldMatrix<Complex>) Ain.copy();

        double[] AcolNorm = new double[nSrc];

        for (int s = 0; s < nSrc; s++) {
            AcolNorm[s] = 0.0;
            for (int c = 0; c < nChannels; c++) {
                AcolNorm[s] += pow(Ain.getEntry(c, s).abs(), 2);
            }
            AcolNorm[s] = sqrt(AcolNorm[s] + EPSILON);
        }

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                Aout.multiplyEntry(c, s, Complex.ONE.divide(AcolNorm[s]));
            }
        }

        return Aout;
    }




    private void runThisFreq(int f) {
        //Log.i("DEBUG", "bin = " + f);

        /* normalizing the columns of X */

        Array2DRowFieldMatrix<Complex> X_bar = normalizeX(STFTwhite[f]);

        /* Initializing mixing matrix */
        Array2DRowFieldMatrix<Complex> A = randomInit();

        /* decouple the rows of A */
        if (isRowDecoupling) {
            ComplexSingularValueDecomposition csvd = new ComplexSingularValueDecomposition(A, true);
            A = csvd.getU().multiply(csvd.getVH());
        }

        /* normalizing A */
        A = normalizeA(A);

        if (derivCheck) {
            checkDerivative(A, X_bar);
        }

        /* Nesterov's accelerated gradient */
        Array2DRowFieldMatrix<Complex> bestA = optimize(A, X_bar);

        if (isRowDecoupling) {
            ComplexSingularValueDecomposition csvd = new ComplexSingularValueDecomposition(bestA, true);
            bestA = csvd.getU().multiply(csvd.getVH());
        }

        /* Normalizing bestA */
        bestA = normalizeA(bestA);

        Array2DRowFieldMatrix<Complex> finalAngle = transjugate(bestA).multiply(X_bar); // (nSrc by nChannels) * (nChannels by nFrames) = nSrc by nFrames

        posteriorProb[f] = calculateMask(finalAngle).getData();
    }

    private Array2DRowFieldMatrix<Complex> optimize(Array2DRowFieldMatrix<Complex> A, Array2DRowFieldMatrix<Complex> X_bar){
        double cost = 0.0;
        double prevCost = Double.POSITIVE_INFINITY;
        double bestCost = Double.POSITIVE_INFINITY;
        double opt_cond;

        Array2DRowFieldMatrix<Complex> grad;
        Array2DRowFieldMatrix<Complex> A_nag = (Array2DRowFieldMatrix<Complex>) A.copy();
        Array2DRowFieldMatrix<Complex> bestA = (Array2DRowFieldMatrix<Complex>) A.copy();
        Array2DRowFieldMatrix<Complex> A_prev_nag, oldA, momentum;
        Complex gamma, complexRestart;

        int counter = 1;

        PowerMeanContrast contrast = new PowerMeanContrast(-0.5, isRowDecoupling);

        for (int iter = 0; iter < maxItr; iter++) {
            oldA = (Array2DRowFieldMatrix<Complex>) A.copy();
            A_prev_nag = (Array2DRowFieldMatrix<Complex>) A_nag.copy();

            if (iter > 0) {
                prevCost = cost;
            }

            gamma = new Complex(((double) counter - 1.0) / ((double) counter + 2.0));

            contrast.computeCost(A, X_bar);
            cost = contrast.getCost();
            grad = contrast.getGrad();

            A_nag = A.subtract((Array2DRowFieldMatrix<Complex>) grad.scalarMultiply(eta));
            A = (Array2DRowFieldMatrix<Complex>) (A_nag.scalarMultiply(Complex.ONE.add(gamma))).subtract(A_prev_nag.scalarMultiply(gamma));

            momentum = A_nag.subtract(A_prev_nag);
            complexRestart = Complex.ZERO;
            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    complexRestart.add(grad.getEntry(c, s).conjugate().multiply(momentum.getEntry(c, s)));
                }
            }

            /* calculating opt_cond */

            opt_cond = grad.getEntry(0, 0).abs();

            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    if (opt_cond < grad.getEntry(c, s).abs()) {
                        opt_cond = grad.getEntry(c, s).abs();
                    }
                }
            }

            if (cost < bestCost) {
                bestCost = cost;
                bestA = (Array2DRowFieldMatrix<Complex>) oldA.copy();
            }

            if (iter > 0) {
                if (abs(prevCost - cost) < PROGRESS_TOLERANCE) {
                    //Log.i("DEBUG", "freq " + f + " completed at iter " + iter + "(PROG_TOL");
                    break;
                } else if (opt_cond < OPTIMUM_TOLERANCE) {
                    //Log.i("DEBUG", "freq " + f + " completed at iter " + iter + "(OPT_TOL");
                    break;
                }
            }

            if (complexRestart.getReal() > 0) {
                counter = 1;
                A = (Array2DRowFieldMatrix<Complex>) A_nag.copy();

                //Log.i("DEBUG", "freq " + f + " restarts at iter " + iter);
            } else {
                counter++;
            }
        }

        return bestA;
    }

    private Array2DRowRealMatrix calculateMask(Array2DRowFieldMatrix<Complex> finalAngle){
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

        return  mask;
    }

    private void checkDerivative(Array2DRowFieldMatrix<Complex> A, Array2DRowFieldMatrix<Complex> X_bar) {

        PowerMeanContrast contrast = new PowerMeanContrast(-0.5, isRowDecoupling);


            /*Log.i("DEBUG", "Testing contrast function");

            double[][] Y2array = new double[][]{{0.0846, -1.3256}, {0.7291, -1.5246}, {0.9008, 0.5024}};

            Array2DRowRealMatrix Y2test = new Array2DRowRealMatrix(Y2array);

            contrast.contrast(Y2test);

            Log.i("DEBUG", "Expected g = 1.4707, Actual g = " + contrast.g);

            Log.i("DEBUG", "Expected dg = [-0.0444, -0.1285; - 0.2757, - 0.1136; - 1.2431, - 1.2988]");
            Log.i("DEBUG", "Actual dg = " + Arrays.deepToString(contrast.dg));*/

        Log.i("DEBUG", "Testing Derivative");

            /*double[][] realA = new double[][]{{-1.0891, 0.5525, 1.5442}, {0.0326, 1.1006, 0.0859}};
            double[][] imagA = new double[][]{{-1.4916, -1.0616, -0.6156}, {-0.7423, 2.3505, 0.7481}};

            Complex[][] testAarray = new Complex[2][3];

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 3; j++) {
                    testAarray[i][j] = new Complex(realA[i][j], imagA[i][j]);
                }
            }

            Array2DRowFieldMatrix<Complex> testA = new Array2DRowFieldMatrix<>(testAarray);

            Log.i("DEBUG", "test A = " + testA);*/

        contrast.computeCost(A, X_bar);

        Array2DRowFieldMatrix<Complex> derivCheckGrad = contrast.getGrad();

        Complex[][] randDirection = new Complex[nChannels][nSrc];

        for (int c = 0; c < nChannels; c++) {
            for (int s = 0; s < nSrc; s++) {
                randDirection[c][s] = new Complex(signum(random() - 0.5));
            }
        }

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

        Log.i("DEBUG", "grad = " + derivCheckGrad.toString());
        Log.i("DEBUG", "gradVec = " + Arrays.toString(gradVec));
        Log.i("DEBUG", "rand = " + Arrays.deepToString(randDirection));
        Log.i("DEBUG", "directVec = " + Arrays.toString(directVec));

        double gtd1 = 0.0;

        for (int i = 0; i < 2 * nChannels * nSrc; i++) {
            gtd1 += gradVec[i] * directVec[i];
        }

        double doubleEps = 1e-8;
        Complex eps = new Complex(doubleEps);

            /*doubleEps = sqrt(doubleEps);
            eps = eps.sqrt();*/

        Array2DRowFieldMatrix<Complex> epsRandDirection = (Array2DRowFieldMatrix<Complex>) new Array2DRowFieldMatrix<>(randDirection).scalarMultiply(eps);

        contrast.computeCost(A.add(epsRandDirection), X_bar);
        double costPlus = contrast.getCost();
        Log.i("DEBUG", "costPlus = " + costPlus);

        contrast.computeCost(A.subtract(epsRandDirection), X_bar);
        double costMinus = contrast.getCost();
        Log.i("DEBUG", "costMinus = " + costMinus);

        double gtd2 = (costPlus - costMinus) / (2 * doubleEps);

        double relDiff = abs(gtd1 - gtd2) / (abs(gtd1) + abs(gtd2));

        Log.i("DEBUG", "analytical deriv = " + gtd1);
        Log.i("DEBUG", "numerical deriv = " + gtd2);
        Log.i("DEBUG", "Relative difference between analytical and numerical directional-derivative is " + relDiff);

    }

    public class PowerMeanContrast {

        double double_r, double_inv_r, double_1minusr;
        int int_r, int_inv_r, int_1minusr;

        boolean isrIntegerValued, isInvrIntegerValued, isrNegInf = false;

        double g = 0.0;
        double[][] dg;

        double cost;

        boolean isRowDecoupling;

        ComplexSingularValueDecomposition csvd;
        Array2DRowFieldMatrix<Complex> grad;

        Array2DRowFieldMatrix<Complex> Sinv, C;
        Complex[] Sdiag;
        Array2DRowFieldMatrix<Complex> U, UH, V, VH, S, SS;

        int nChannels, nSrc, nFrames;

        PowerMeanContrast(double r, boolean isRowDecoupling) {

            if (r == Double.NEGATIVE_INFINITY) {
                isrNegInf = true;
            } else if (r % 1 == 0) {
                isrIntegerValued = true;
                this.int_r = (int) r;
                this.int_1minusr = (int) (1.0 - r);
            } else {
                isrIntegerValued = false;
                this.double_r = r;
                this.double_1minusr = 1.0 - r;
            }

            if ((1.0 / r) % 1 == 0) {
                isInvrIntegerValued = true;
                this.int_inv_r = (int) (1.0 / r);
                double_inv_r = int_inv_r;
            } else {
                isInvrIntegerValued = false;
                double_inv_r = 1.0 / r;
            }

            this.isRowDecoupling = isRowDecoupling;
        }

        private void computeCost(Array2DRowFieldMatrix<Complex> A, Array2DRowFieldMatrix<Complex> X_bar) {

            nChannels = A.getRowDimension();
            nSrc = A.getColumnDimension();
            nFrames = X_bar.getColumnDimension();

            dg = new double[nSrc][nFrames];

            Array2DRowFieldMatrix<Complex> Aorth;
            Complex[][] AorthArray = A.getData();

            if (isRowDecoupling) {
                csvd = new ComplexSingularValueDecomposition(A, true);

                U = csvd.getU();
                UH = csvd.getUH();

                V = csvd.getV();
                VH = csvd.getVH();

                S = csvd.getS();
                Sinv = csvd.getPoweredS(-1);
                Sdiag = csvd.getSentryComplex(); //size = nChannels

                Aorth = U.multiply(VH);
                AorthArray = Aorth.getData();
            }

            /* Column norm normalization */

            double[] colNorm = new double[nSrc];
            Complex[][] AnormArray = new Complex[nChannels][nSrc];

            for (int s = 0; s < nSrc; s++) {
                colNorm[s] = 0.0;
                for (int c = 0; c < nChannels; c++) {
                    colNorm[s] += pow(AorthArray[c][s].abs(), 2);
                }
                colNorm[s] = sqrt(colNorm[s] + EPSILON);
            }

            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    AnormArray[c][s] = AorthArray[c][s].divide(colNorm[s]);
                }
            }

            Array2DRowFieldMatrix<Complex> Anorm = new Array2DRowFieldMatrix<>(AnormArray);
            Array2DRowFieldMatrix<Complex> angle = transjugate(Anorm).multiply(X_bar); //nSrc by nFrames

            /* calculating contrast function */
            Array2DRowRealMatrix Y2 = new Array2DRowRealMatrix(nSrc, nFrames);
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    Y2.setEntry(s, t, pow(angle.getEntry(s, t).abs(), 2));
                }
            }

            contrast(Y2); //return values go into g and dg

            /* scaling */
            double scale = 1.0;
            cost = (scale / nFrames) * g;

            Array2DRowFieldMatrix<Complex> angle_dg = (Array2DRowFieldMatrix<Complex>) angle.copy();
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    angle_dg.multiplyEntry(s, t, new Complex(dg[s][t]));
                }
            }

            Array2DRowFieldMatrix<Complex> temp_grad;
            temp_grad = (Array2DRowFieldMatrix<Complex>) X_bar.multiply(transjugate(angle_dg)).scalarMultiply(new Complex(2.0 * scale / nFrames));

            /* Apply the chain rule through column norm normalization */

            Complex[] Aorth_multiplier = new Complex[nSrc];
            for (int s = 0; s < nSrc; s++) {
                Aorth_multiplier[s] = Complex.ZERO;
                for (int c = 0; c < nChannels; c++) {
                    Aorth_multiplier[s] = Aorth_multiplier[s].add(AorthArray[c][s].conjugate().multiply(temp_grad.getEntry(c, s)));
                }
                Aorth_multiplier[s] = Aorth_multiplier[s].divide(pow(colNorm[s], 3));
            }

            Array2DRowFieldMatrix<Complex> tempGrad1 = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nSrc);
            Array2DRowFieldMatrix<Complex> tempGrad2 = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nSrc);

            for (int c = 0; c < nChannels; c++) {
                for (int s = 0; s < nSrc; s++) {
                    tempGrad1.setEntry(c, s, temp_grad.getEntry(c, s).divide(colNorm[s]));
                    tempGrad2.setEntry(c, s, AorthArray[c][s].multiply(Aorth_multiplier[s]));
                }
            }

            grad = tempGrad1.subtract(tempGrad2);

            if (isRowDecoupling) {
                SS = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), nChannels, nChannels);
                for (int ci = 0; ci < nChannels; ci++) {
                    for (int cj = 0; cj < nChannels; cj++) {
                        SS.setEntry(ci, cj, Sdiag[ci].add(Sdiag[cj]));
                    }
                }

                C = Sinv.multiply(UH).multiply(grad).multiply(V);
                for (int ci = 0; ci < nChannels; ci++) {
                    for (int cj = 0; cj < nChannels; cj++) {
                        C.multiplyEntry(ci, cj, Complex.ONE.divide((SS.getEntry(ci, cj).negate()).add(EPSILON)));
                    }
                }

                grad = U.multiply(transjugate(C).add(C)).multiply(S).multiply(VH)
                        .add(U.multiply(Sinv).multiply(UH).multiply(grad));
            }
        }

        Array2DRowFieldMatrix<Complex> getGrad() {
            return grad;
        }

        private double getCost() {
            return cost;
        }

        private void contrast(Array2DRowRealMatrix Y2) {
            Array2DRowRealMatrix Y2_contra_mat = new Array2DRowRealMatrix(nSrc, nFrames);
            double[][] Y2_contra; // = new double[nFrames][nSrc]
            double[] Y2_contra_min = new double[nFrames];
            double[] Y2_power_mean_contra = new double[nFrames];

            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    Y2_contra_mat.setEntry(s, t, 1.0 - Y2.getEntry(s, t));
                }
            }

            Y2_contra = Y2_contra_mat.getData();

            double colMin;
            double[] thisCol;
            int[] colMinIdx = new int[nFrames];
            int[] colZeroIdx = new int[nFrames];
            Min min = new Min();

            for (int t = 0; t < nFrames; t++) {
                thisCol = Y2_contra_mat.getColumn(t);
                colMin = min.evaluate(thisCol);
                colZeroIdx[t] = ArrayUtils.indexOf(thisCol, 0.0); //returns -1 if there is no zero

                if (isrNegInf) {
                    colMinIdx[t] = ArrayUtils.indexOf(thisCol, colMin);
                }

                Y2_contra_min[t] = colMin;
            }

            Sum sum = new Sum();

            if (isrNegInf) {
                g = sum.evaluate(Y2_contra_min);
                for (int s = 0; s < nSrc; s++) {
                    for (int t = 0; t < nFrames; t++) {
                        if (colMinIdx[t] > -1) { // if there is a zero
                            dg[colMinIdx[t]][t] = -1;
                        }
                    }
                }
            } else {
                double thisFramePowerSum;
                g = 0.0;

                double nSrcPowered;
                if (isInvrIntegerValued) {
                    nSrcPowered = pow(nSrc, int_inv_r);
                } else {
                    nSrcPowered = pow(nSrc, double_inv_r);
                }

                for (int t = 0; t < nFrames; t++) {
                    if (colZeroIdx[t] > -1) {
                        Y2_power_mean_contra[t] = 0;
                        continue;
                    }

                    thisFramePowerSum = 0.0;

                    if (isrIntegerValued) {
                        for (int s = 0; s < nSrc; s++) {
                            thisFramePowerSum += pow(Y2_contra[s][t] / Y2_contra_min[t], int_r);
                        }
                    } else {
                        for (int s = 0; s < nSrc; s++) {
                            thisFramePowerSum += pow(Y2_contra[s][t] / Y2_contra_min[t], double_r);
                        }
                    }

                    if (isInvrIntegerValued) {
                        thisFramePowerSum = pow(thisFramePowerSum, int_inv_r);
                        Y2_power_mean_contra[t] = thisFramePowerSum * Y2_contra_min[t] / nSrcPowered;
                    } else {
                        thisFramePowerSum = pow(thisFramePowerSum, double_inv_r);
                        Y2_power_mean_contra[t] = thisFramePowerSum * Y2_contra_min[t] / nSrcPowered;
                    }

                    g += Y2_power_mean_contra[t];
                }

                double zeroDg = -1.0 / pow(nSrc, double_inv_r);

                if (isrIntegerValued) {
                    for (int s = 0; s < nSrc; s++) {
                        for (int t = 0; t < nFrames; t++) {
                            if (colZeroIdx[t] == s) {
                                for (int tt = 0; tt < nFrames; tt++) {
                                    dg[s][tt] = zeroDg;
                                }
                                break; //no need to iterate thru the frames anymore
                            }

                            dg[s][t] = -pow(Y2_power_mean_contra[t] / Y2_contra[s][t], int_1minusr) / (double) nSrc;
                        }
                    }
                } else {
                    for (int s = 0; s < nSrc; s++) {
                        for (int t = 0; t < nFrames; t++) {
                            if (colZeroIdx[t] == s) {
                                for (int tt = 0; tt < nFrames; tt++) {
                                    dg[s][tt] = zeroDg;
                                }
                                break; //no need to iterate thru the frames anymore
                            }
                            dg[s][t] = -pow(Y2_power_mean_contra[t] / Y2_contra[s][t], double_1minusr) / (double) nSrc;
                        }
                    }
                }
            }
        }
    }

    public Complex[][][] getSourceEstimatesSTFT() {
        return STFTout;
    }
}
