package com.example.administrator.audiorecord.audioprocessing.commons;

import android.util.Log;

import com.example.administrator.audiorecord.audioprocessing.bss.lapjv.LAPJV;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

import java.util.Arrays;
import java.util.HashSet;
import java.util.function.Predicate;
import java.util.stream.IntStream;

import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.round;
import static org.apache.commons.math3.util.FastMath.sqrt;

public class PermutationAlignment {

    private double[][][] P, Pnorm;

    private int maxItr;
    private double tolerance;

    private int nFreqs, nFrames, nSrc;

    private final double EPSILON = 1e-8;

    private boolean fineTune;

    public PermutationAlignment(final double[][][] P, int maxItr, double tolerance, boolean fineTune) {
        this.P = SerializationUtils.clone(P);
        this.maxItr = maxItr;
        this.tolerance = tolerance;
        this.fineTune = fineTune;

        nFreqs = P.length;
        nFrames = P[0].length;
        nSrc = P[0][0].length;
    }

    public void run() {
        StandardDeviation stdev = new StandardDeviation(true);

        double[][] PtimeAvg = new double[nFreqs][nSrc];
        Pnorm = new double[nFreqs][nFrames][nSrc];

        IntStream.range(0, nFreqs).parallel().forEach(f -> {
            for (int s = 0; s < nSrc; s++) {
                for (int t = 0; t < nFrames; t++) {
                    PtimeAvg[f][s] += P[f][t][s];
                }
                PtimeAvg[f][s] = PtimeAvg[f][s] / (double) nFrames;
            }

            for (int t = 0; t < nFrames; t++) {
                for (int s = 0; s < nSrc; s++) {
                    Pnorm[f][t][s] = P[f][t][s] - PtimeAvg[f][s];
                }
            }
        });

        int step = 2;

        while (step < nFreqs) {

//            Log.i("DEBUG", "Step = " + step);

            int thisStep = step;
            IntStream.range(0, nFreqs / step).map(i -> thisStep * i + 1).parallel().forEach(i -> { // equivalent to for i = 1 : step : nFreqs - 1

//                Log.i("DEBUG", "i = " + i);
                //Log.i("DEBUG", "Group 1: from " + i + " to " + (i + thisStep / 2 - 1));

                int[] group1idx = IntStream.range(i, i + thisStep / 2).toArray();
                int[] group2idx;

                if ((thisStep > 2) && (i + thisStep == nFreqs - 2)) {
                    group2idx = IntStream.rangeClosed(i + thisStep / 2, nFreqs - 1).toArray();
                    //Log.i("DEBUG", "Group 2*: from " + (i + thisStep / 2) + " to " + (nFreqs - 1));
                } else {
                    group2idx = IntStream.range(i + thisStep / 2, i + thisStep).toArray();
                    //Log.i("DEBUG", "Group 2: from " + (i + thisStep / 2) + " to " + (i + thisStep - 1));
                }

                double[][] f_group1 = new double[nSrc][nFrames]; //initialize to zeros, nSrc by nFrames for speed when calculating stdev
                double[][] f_group2 = new double[nSrc][nFrames]; //initialize to zeros, nSrc by nFrames for speed when calculating stdev

                for (int t = 0; t < nFrames; t++) {
                    for (int s = 0; s < nSrc; s++) {
                        //f_group1[s][t] = 0.0;
                        for (int f : group1idx) {
                            f_group1[s][t] += Pnorm[f][t][s];
                        }
                        f_group1[s][t] = f_group1[s][t] / (double) group1idx.length;
                    }
                }

                for (int t = 0; t < nFrames; t++) {
                    for (int s = 0; s < nSrc; s++) {
                        //f_group2[s][t] = 0.0;
                            for (int f : group2idx) {
                                f_group2[s][t] += Pnorm[f][t][s];
                            }
                        f_group2[s][t] = f_group2[s][t] / (double) group2idx.length;
                    }
                }

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

                Array2DRowRealMatrix f_group1mat = new Array2DRowRealMatrix(f_group1); //nSrc by nFrames
                Array2DRowRealMatrix f_group2mat = new Array2DRowRealMatrix(f_group2); //nSrc by nFrames

                Array2DRowRealMatrix Qmat = (Array2DRowRealMatrix) f_group1mat.multiply(f_group2mat.transpose()).scalarMultiply(1.0 / (double) nFrames);

                for (int si = 0; si < nSrc; si++) {
                    for (int sj = 0; sj < nSrc; sj++) {
                        Qmat.setEntry(si, sj, 1.0 - Qmat.getEntry(si, sj));
                    }
                }

                rearrange(Qmat, group1idx);
            });
            step *= 2;
        }

        if (fineTune) {
            int[][] fineTuningFreqArray = new int[nFreqs][];
            int[] fineTuningFreqCount = new int[nFreqs];

            double relError = 0.0;

            double[] numf = new double[nFreqs], denomf = new double[nFreqs];
            double num, denom;

            //Prepare finetune frequencies once only
            IntStream.range(1, nFreqs).parallel().forEach(f -> {
                Predicate<Integer> illegalFreq = ff -> (ff < 1);
                illegalFreq = illegalFreq.or(ff -> (ff > nFreqs - 1)).or(ff -> (ff == f));

                HashSet<Integer> fineTuningFreq = new HashSet(30); //maximum number of fine tune frequencies
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

                //Log.i("DEBUG", "finetune " + f + ": " + Arrays.toString(fineTuningFreqArray[f]));
            });

            //int additionalIter = 0;

            for (int iter = 0; iter < maxItr; iter++) {
                //Log.i("DEBUG", "Fine tuning iter " + iter);
                double[][][] oldPnorm = SerializationUtils.clone(Pnorm);

                IntStream.range(1, nFreqs).parallel().forEach(f -> { //.map(f -> nFreqs - f) || no difference since we are using oldPnorm
                    //Log.i("DEBUG", "f = " + f);
                    Array2DRowRealMatrix pf_mat, pf_fineTuning_mat;
                    double stdev_pf;

                    double[][] pf_fineTuning = new double[nSrc][nFrames]; //nSrc by nFrames, for speed when calculating stdev

                    for (int t = 0; t < nFrames; t++) {
                        for (int s = 0; s < nSrc; s++) {
                            synchronized (oldPnorm) {
                                for (int ff : fineTuningFreqArray[f]) {
                                    pf_fineTuning[s][t] += oldPnorm[ff][t][s];
                                }
                            }
                            pf_fineTuning[s][t] = pf_fineTuning[s][t] / fineTuningFreqCount[f];
                        }
                    }

                    for (int s = 0; s < nSrc; s++) {
                        stdev_pf = stdev.evaluate(pf_fineTuning[s]);
                        for (int t = 0; t < nFrames; t++) {
                            pf_fineTuning[s][t] = pf_fineTuning[s][t] / (stdev_pf + EPSILON);
                        }
                    }

                    pf_mat = new Array2DRowRealMatrix(oldPnorm[f]); //nFrames by nSrc
                    pf_fineTuning_mat = new Array2DRowRealMatrix(pf_fineTuning); //nSrc by nFrames

                    Array2DRowRealMatrix thisQmat = (Array2DRowRealMatrix) (pf_fineTuning_mat.multiply(pf_mat)).transpose().scalarMultiply(1.0 / (double) nFrames);
                    //same as pf_mat' * pf_finetuning_mat'

                    for (int si = 0; si < nSrc; si++) {
                        for (int sj = 0; sj < nSrc; sj++) {
                            thisQmat.setEntry(si, sj, 1.0 - thisQmat.getEntry(si, sj));
                        }
                    }

                    rearrangeFinetune(thisQmat, f, oldPnorm);
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

                Sum sum = new Sum();

                num = sum.evaluate(numf);
                denom = sum.evaluate(denomf);

                relError = sqrt(num / denom);

                //Log.i("DEBUG", "relError = " + relError);

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

        IntStream.range(0, nSrc).parallel().forEach(s -> {
            double stdev_pf_ac = stdev.evaluate(pf_ac[s]);
            for (int t = 0; t < nFrames; t++) {
                pf_ac[s][t] = pf_ac[s][t] / (stdev_pf_ac + EPSILON);
            }
        });

        Array2DRowRealMatrix pf_dc_mat = new Array2DRowRealMatrix(Pnorm[0]); //nFrames by nSrc
        Array2DRowRealMatrix pf_ac_mat = new Array2DRowRealMatrix(pf_ac); //nSrc by nFrames

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

        P[0] = SerializationUtils.clone(tempP);
    }

    public double[][][] getP() {
        return SerializationUtils.clone(P);
    }

    private synchronized void rearrange(Array2DRowRealMatrix Qmat, int[] indices) {

        int[] per = LAPJV.execute(Qmat.getData());

        for (int f : indices) {
            double[][] tempP = new double[nFrames][nSrc];
            double[][] tempPnorm = new double[nFrames][nSrc];

            for (int t = 0; t < nFrames; t++) {
                for (int s = 0; s < nSrc; s++) {
                    tempP[t][per[s]] = P[f][t][s];
                    tempPnorm[t][per[s]] = Pnorm[f][t][s];
                }
            }

            P[f] = SerializationUtils.clone(tempP);
            Pnorm[f] = SerializationUtils.clone(tempPnorm);
        }
    }

    private synchronized void rearrangeFinetune(Array2DRowRealMatrix Qmat, int f, double[][][] oldPnorm) {

        int[] per = LAPJV.execute(Qmat.getData());

        double[][] tempP = new double[nFrames][nSrc];
        double[][] tempPnorm = new double[nFrames][nSrc];

        for (int t = 0; t < nFrames; t++) {
            for (int s = 0; s < nSrc; s++) {
                tempP[t][per[s]] = P[f][t][s];
                tempPnorm[t][per[s]] = oldPnorm[f][t][s];
            }
        }

        P[f] = SerializationUtils.clone(tempP);
        Pnorm[f] = SerializationUtils.clone(tempPnorm);
    }

    /*private synchronized int[] performLAPJV (Array2DRowRealMatrix Qmat){

        double[] cc = new double[nSrc*nSrc];
        int[] kk = new int[nSrc*nSrc];
        int[] nInEachRow = new int[nSrc];

        for(int si = 0; si < nSrc; si++){
            for(int sj = 0; sj < nSrc; sj++){
                cc[si*nSrc + sj] = Qmat.getEntry(si,sj);
            }
        }

        for(int ss = 0; ss < nSrc*nSrc; ss++){
            kk[ss] = ss%3;
        }

        for(int s = 0; s < nSrc; s++){
            nInEachRow[s] = nSrc;
        }

        SparseCostMatrix Qcm = new SparseCostMatrix(cc,kk,nInEachRow,nSrc);

        LAPJV lapjv = new LAPJV(Qcm);
        int[] out;
        boolean isSuccess = lapjv.process();

        if(isSuccess){
            out = lapjv.getResult();
            return out;
        } else {
            out = IntStream.range(0,nSrc).toArray(); //keep same order
            return out;
        }
    }*/


}
