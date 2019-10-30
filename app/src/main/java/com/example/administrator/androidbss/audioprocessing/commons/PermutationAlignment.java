package com.example.administrator.androidbss.audioprocessing.commons;

import android.util.Log;

import com.example.administrator.androidbss.audioprocessing.commons.lapjv.LAPJV;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

import java.util.Arrays;
import java.util.HashSet;
import java.util.function.Predicate;
import java.util.stream.IntStream;

import static com.example.administrator.androidbss.audioprocessing.commons.Epsilon.EPSILON;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.round;
import static org.apache.commons.math3.util.FastMath.sqrt;

public class PermutationAlignment {

    private double[][][] P, Pnorm;

    private int maxItr = 10;
    private double tolerance = 1e-3;

    private int nFreqs, nFrames, nSrc;
    private boolean fineTune = true;

    private PermutationAlignment() {

    }

    public static PermutationAlignment initiate() {
        return new PermutationAlignment();
    }

    public PermutationAlignment setMaximumIterations(int maxItr) {
        this.maxItr = maxItr;
        return this;
    }

    public PermutationAlignment setProgressTolerance(double tolerance) {
        this.tolerance = tolerance;
        return this;
    }

    public PermutationAlignment setEnabledFinetune(boolean fineTune) {
        this.fineTune = fineTune;
        return this;
    }

    public PermutationAlignment setInputData(double[][][] P) {
        nFreqs = P.length;
        nFrames = P[0].length;
        nSrc = P[0][0].length;

        this.P = new double[nFreqs][nFrames][nSrc];

        IntStream.range(0,nFreqs).parallel().forEach(f->{
            for (int t = 0; t < nFrames; t++) {
                System.arraycopy(P[f][t], 0, this.P[f][t], 0, nSrc);
            }
        });

        return this;
    }

    private void normalizeP() {
        Pnorm = new double[nFreqs][nFrames][nSrc];

        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            double[] PtimeAvg = new double[nSrc];
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    PtimeAvg[s] += P[f][t][s];
                }
                PtimeAvg[s] = PtimeAvg[s] / (double) nFrames;
            }

            for (int t = 0; t < nFrames; t++) {
                for (int s = 0; s < nSrc; s++) {
                    Pnorm[f][t][s] = P[f][t][s] - PtimeAvg[s];
                }
            }
        });
    }

    private double[][] calculateCentroid(int[] groupIdx) {
        double[][] f_group = new double[nSrc][nFrames];

        for (int t = 0; t < nFrames; t++) {
            for (int s = 0; s < nSrc; s++) {
                for (int f : groupIdx) {
                    f_group[s][t] += Pnorm[f][t][s];
                }
                f_group[s][t] = f_group[s][t] / (double) groupIdx.length;
            }
        }

        return f_group;
    }


    private void alignThisSubband(int i, int step) {
        StandardDeviation stdev = new StandardDeviation(true);

        int[] group1idx = IntStream.range(i, i + step / 2).toArray();
        int[] group2idx;

        if ((step > 2) && (i + step == nFreqs - 2)) {
            group2idx = IntStream.rangeClosed(i + step / 2, nFreqs - 1).toArray();
        } else {
            group2idx = IntStream.range(i + step / 2, i + step).toArray();
        }

        double[][] f_group1 = calculateCentroid(group1idx);
        double[][] f_group2 = calculateCentroid(group2idx);

        double[] stdev_group1 = new double[nSrc];
        double[] stdev_group2 = new double[nSrc];

        for (int s = 0; s < nSrc; s++) {
            stdev_group1[s] = stdev.evaluate(f_group1[s]);
            stdev_group2[s] = stdev.evaluate(f_group2[s]);
        }

        for (int s = 0; s < nSrc; s++) {
            for (int t = 0; t < nFrames; t++) {
                f_group1[s][t] = f_group1[s][t] / (stdev_group1[s] + EPSILON);
                f_group2[s][t] = f_group2[s][t] / (stdev_group2[s] + EPSILON);
            }
        }
        Array2DRowRealMatrix f_group1mat = new Array2DRowRealMatrix(f_group1, false); //nSrc by nFrames
        Array2DRowRealMatrix f_group2mat = new Array2DRowRealMatrix(f_group2, false); //nSrc by nFrames
        Array2DRowRealMatrix Qmat = (Array2DRowRealMatrix) f_group1mat.multiply(f_group2mat.transpose()).scalarMultiply(1.0 / (double) nFrames);

        for (int si = 0; si < nSrc; si++) {
            for (int sj = 0; sj < nSrc; sj++) {
                Qmat.setEntry(si, sj, 1.0 - Qmat.getEntry(si, sj));
            }
        }

        rearrange(Qmat, group1idx);
    }

    private void prepareFinetuneFrequencies(int[][] fineTuningFreqArray, int[] fineTuningFreqCount) {
        IntStream.range(1, nFreqs).parallel().forEach(f -> {
            Predicate<Integer> illegalFreq = ff -> (ff < 1);
            illegalFreq = illegalFreq.or(ff -> (ff > nFreqs - 1)).or(ff -> (ff == f));

            HashSet<Integer> fineTuningFreq = new HashSet(); //maximum number of fine tune frequencies
            int underF, overF;

            int multiplier = 2;

            for (int e = 1; e < 6; e++) {
                underF = (int) round((double) f / multiplier);
                fineTuningFreq.add(underF - 1);
                fineTuningFreq.add(underF);
                fineTuningFreq.add(underF + 1);

                overF = multiplier * (f);
                fineTuningFreq.add(overF - 1);
                fineTuningFreq.add(overF);
                fineTuningFreq.add(overF + 1);

                multiplier *= 2;
            }

            fineTuningFreq.add(f - 3);
            fineTuningFreq.add(f - 2);
            fineTuningFreq.add(f - 1);
            fineTuningFreq.add(f + 1);
            fineTuningFreq.add(f + 2);
            fineTuningFreq.add(f + 3);

            fineTuningFreq.removeIf(illegalFreq);

            //Array is faster than HashSet
            fineTuningFreqArray[f] = ArrayUtils.toPrimitive(fineTuningFreq.toArray(new Integer[0]));
            fineTuningFreqCount[f] = fineTuningFreqArray[f].length;
        });

        Log.i("DEBUG", Arrays.deepToString(fineTuningFreqArray));
    }

    private void fineTuneThisFreq(int f, double[][][] oldPnorm, int[] fineTuningFreqArray, int fineTuningFreqCount) {
        StandardDeviation stdev = new StandardDeviation(true);
        double stdev_pf;

        double[][] pf_fineTuning = new double[nSrc][nFrames]; //nSrc by nFrames, for speed when calculating stdev
        for (int t = 0; t < nFrames; t++) {
            for (int s = 0; s < nSrc; s++) {
                for (int ff : fineTuningFreqArray) {
                    pf_fineTuning[s][t] += oldPnorm[ff][t][s];
                }
                pf_fineTuning[s][t] = pf_fineTuning[s][t] / fineTuningFreqCount;
            }
        }

        for (int s = 0; s < nSrc; s++) {
            stdev_pf = stdev.evaluate(pf_fineTuning[s]);
            for (int t = 0; t < nFrames; t++) {
                pf_fineTuning[s][t] = pf_fineTuning[s][t] / (stdev_pf + EPSILON);
            }
        }
        Array2DRowRealMatrix pf_mat = new Array2DRowRealMatrix(oldPnorm[f], false); //nFrames by nSrc
        Array2DRowRealMatrix pf_fineTuning_mat = new Array2DRowRealMatrix(pf_fineTuning, false); //nSrc by nFrames

        Array2DRowRealMatrix thisQmat = (Array2DRowRealMatrix) (pf_fineTuning_mat.multiply(pf_mat)).transpose().scalarMultiply(1.0 / (double) nFrames);

        for (int si = 0; si < nSrc; si++) {
            for (int sj = 0; sj < nSrc; sj++) {
                thisQmat.setEntry(si, sj, 1.0 - thisQmat.getEntry(si, sj));
            }
        }

        rearrangeFinetune(thisQmat, f, oldPnorm);
    }


    public double[][][] run() {

        normalizeP();

        int step = 2;
        while (step < nFreqs) {
            int thisStep = step;
            IntStream.range(0, nFreqs / step).map(i -> thisStep * i + 1).parallel().forEach(i -> { // equivalent to for i = 1 : step : nFreqs - 1
                alignThisSubband(i, thisStep);
            });
            step *= 2;
        }

        if (fineTune) {

            int[][] fineTuningFreqArray = new int[nFreqs][];
            int[] fineTuningFreqCount = new int[nFreqs];
            prepareFinetuneFrequencies(fineTuningFreqArray, fineTuningFreqCount);

            double[][][] oldPnorm = new double[nFreqs][nFrames][nSrc];

            double[] numf = new double[nFreqs], denomf = new double[nFreqs];
            double relError = 0.0, num, denom;
            Sum sum = new Sum();

            for (int iter = 0; iter < maxItr; iter++) {
                //Log.i("DEBUG", "Fine tuning iter " + iter);
                IntStream.range(1, nFreqs).parallel().forEach(f -> {
                    for (int t = 0; t < nFrames; t++) {
                        System.arraycopy(Pnorm[f][t], 0, oldPnorm[f][t], 0, nSrc);
                    }
                });

                IntStream.range(1, nFreqs).parallel().forEach(f -> {
                    fineTuneThisFreq(f, oldPnorm, fineTuningFreqArray[f], fineTuningFreqCount[f]);
                });

                Arrays.fill(numf, 0.0);
                Arrays.fill(denomf, 0.0);

                IntStream.range(0, nFreqs).parallel().forEach(f -> {
                    for (int t = 0; t < nFrames; t++) {
                        for (int s = 0; s < nSrc; s++) {
                            denomf[f] += pow(oldPnorm[f][t][s] + EPSILON, 2);
                            numf[f] += pow(oldPnorm[f][t][s] - Pnorm[f][t][s], 2);
                        }
                    }
                });

                num = sum.evaluate(numf);
                denom = sum.evaluate(denomf);
                relError = sqrt(num / denom);

//                Log.i("DEBUG", "relError = " + relError);

                if ((iter > 0) && (relError < tolerance)) {
                    //Log.i("DEBUG", "Fine tuning: Error below tolerance.");
                    break;
                }
            }
            Log.i("DEBUG", "Fine tuning completed with relError " + relError);
        }

        double[][] pf_ac = new double[nSrc][nFrames];

        IntStream.range(0, nFrames).parallel().forEach(t -> {
            for (int s = 0; s < nSrc; s++) {
                for (int f = 1; f < nFreqs; f++) {
                    pf_ac[s][t] += Pnorm[f][t][s];
                }
                pf_ac[s][t] = pf_ac[s][t] / (nFreqs - 1);
            }
        });

        for (int s = 0; s < nSrc; s++) {
            StandardDeviation stdev = new StandardDeviation(true);

            double stdev_pf_ac = stdev.evaluate(pf_ac[s]);
            int finalS = s;
            IntStream.range(0, nFrames).parallel().forEach(t -> {
                pf_ac[finalS][t] = pf_ac[finalS][t] / (stdev_pf_ac + EPSILON);
            });
        }

        Array2DRowRealMatrix pf_dc_mat = new Array2DRowRealMatrix(Pnorm[0],false); //nFrames by nSrc
        Array2DRowRealMatrix pf_ac_mat = new Array2DRowRealMatrix(pf_ac,false); //nSrc by nFrames

        Array2DRowRealMatrix Qmat = (Array2DRowRealMatrix) (pf_ac_mat.multiply(pf_dc_mat)).transpose().scalarMultiply(1.0 / (double) nFrames);
        // same as pf_dc_mat' * pf_ac_mat'

        for (int si = 0; si < nSrc; si++) {
            for (int sj = 0; sj < nSrc; sj++) {
                Qmat.setEntry(si, sj, 1.0 - Qmat.getEntry(si, sj));
            }
        }

        int[] per = LAPJV.execute(Qmat.getData());

        double[][] tempP = new double[nFrames][nSrc];
        IntStream.range(0, nFrames).parallel().forEach(t -> {
            for (int s = 0; s < nSrc; s++) {
                tempP[t][per[s]] = P[0][t][s];
            }
        });

        IntStream.range(0, nFrames).parallel().forEach(t -> {
            System.arraycopy(tempP[t], 0, P[0][t], 0, nSrc);
        });



        return P;
    }

    private void rearrange(Array2DRowRealMatrix Qmat, int[] indices) {

        int[] per = LAPJV.execute(Qmat.getData());

        double[][] tempP = new double[nFrames][nSrc];
        double[][] tempPnorm = new double[nFrames][nSrc];

        for (int f : indices) {
            for (int t = 0; t < nFrames; t++) {
                for (int s = 0; s < nSrc; s++) {
                    tempP[t][per[s]] = P[f][t][s];
                    tempPnorm[t][per[s]] = Pnorm[f][t][s];
                }
            }

            for (int t = 0; t < nFrames; t++) {
                System.arraycopy(tempP[t], 0, P[f][t], 0, nSrc);
            }

            for (int t = 0; t < nFrames; t++) {
                System.arraycopy(tempPnorm[t], 0, Pnorm[f][t], 0, nSrc);
            }
        }
    }

    private void rearrangeFinetune(Array2DRowRealMatrix Qmat, int f, double[][][] oldPnorm) {

        int[] per = LAPJV.execute(Qmat.getData());

        double[][] tempP = new double[nFrames][nSrc];
        double[][] tempPnorm = new double[nFrames][nSrc];

        for (int t = 0; t < nFrames; t++) {
            for (int s = 0; s < nSrc; s++) {
                tempP[t][per[s]] = P[f][t][s];
                tempPnorm[t][per[s]] = oldPnorm[f][t][s];
            }
        }

        for (int t = 0; t < nFrames; t++) {
            System.arraycopy(tempP[t], 0, P[f][t], 0, nSrc);
        }

        for (int t = 0; t < nFrames; t++) {
            System.arraycopy(tempPnorm[t], 0, Pnorm[f][t], 0, nSrc);
        }
    }

}
