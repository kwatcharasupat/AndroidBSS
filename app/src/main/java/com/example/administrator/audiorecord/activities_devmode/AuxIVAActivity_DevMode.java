package com.example.administrator.audiorecord.activities_devmode;

import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.Toast;

import com.example.administrator.audiorecord.R;
import com.example.administrator.audiorecord.audioprocessing.bss.AuxIVA;
import com.example.administrator.audiorecord.audioprocessing.commons.STFT;
import com.example.administrator.audiorecord.audioprocessing.datahandler.AudioFileWriter;
import com.example.administrator.audiorecord.audioprocessing.datahandler.STFTParcel;

import org.apache.commons.math3.complex.Complex;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;

import static android.os.Debug.startMethodTracing;
import static android.os.Environment.getExternalStorageDirectory;
import static com.example.administrator.audiorecord.activities_devmode.STFTActivity_DevMode.STFT_PARCEL;
import static com.example.administrator.audiorecord.activities_devmode.STFTActivity_DevMode.audioData;
import static java.lang.Math.abs;

public class AuxIVAActivity_DevMode extends AppCompatActivity {

    Complex[][][] STFTdata, STFTout;
    int audioDataLength, nChannels, nFrames, nFreqs, winLen, nOverlap;
    String winFunc;

    Button runButton, istftButton, writeToFileButton;

    AuxIVA auxIVA;

    STFT stft;
    double[][] reSig;

    Thread auxIVAThread = null;
    Thread istftThread = null;
    Thread writingToFilesThread[];

    AtomicBoolean isAuxIVArunning = new AtomicBoolean(false);
    AtomicBoolean isAuxIVAcompleted = new AtomicBoolean(false);
    AtomicBoolean isISTFTrunning = new AtomicBoolean(false);
    AtomicBoolean isISTFTcompleted = new AtomicBoolean(false);
    AtomicBoolean isWritingToFiles[];
    AtomicBoolean isWritingToFilesCompleted[];

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.devmodeactivity_auxiva);

        getDataFromIntent();

        Log.i("DEBUG", "STFT channels: " + STFTdata.length);
        Log.i("DEBUG", "STFT frames: " + STFTdata[0].length);
        Log.i("DEBUG", "STFT bins: " + STFTdata[0][0].length);

        isWritingToFiles = new AtomicBoolean[nChannels];
        isWritingToFilesCompleted = new AtomicBoolean[nChannels];
        writingToFilesThread = new Thread[nChannels];

        runButton = findViewById(R.id.bRunAuxIVA);
        runButton.setOnClickListener(v -> {
            istftButton.setEnabled(false);
            writeToFileButton.setEnabled(false);
            Log.i("DEBUG", "Running AuxIVA");
            runAuxIVA();
        });

        istftButton = findViewById(R.id.bReconSig);
        istftButton.setOnClickListener(v -> {
            if (isAuxIVAcompleted.get()) {
                Toast toast = Toast.makeText(getApplicationContext(), "AuxIVA is completed. Running ISTFT", Toast.LENGTH_LONG);
                toast.show();

                auxIVAThread.interrupt();
                auxIVAThread = null;

                runISTFT();
            } else if (isAuxIVArunning.get()) {
                Toast toast = Toast.makeText(getApplicationContext(), "AuxIVA is still running", Toast.LENGTH_LONG);
                toast.show();
            } else {
                Toast toast = Toast.makeText(getApplicationContext(), "AuxIVA has not been run", Toast.LENGTH_LONG);
                toast.show();
            }
        });

        writeToFileButton = findViewById(R.id.bSaveSig);
        writeToFileButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (isISTFTcompleted.get()) {

                    istftThread.interrupt();
                    istftThread = null;

                    saveChannelsToFiles();
                } else if (isISTFTrunning.get()) {
                    Toast toast = Toast.makeText(getApplicationContext(), "ISTFT is still running", Toast.LENGTH_LONG);
                    toast.show();
                } else {
                    Toast toast = Toast.makeText(getApplicationContext(), "ISTFT has not been run", Toast.LENGTH_LONG);
                    toast.show();
                }
            }
        });
    }

    void runAuxIVA() {

        //startMethodTracingSampling("AuxIVAThread", 2147483647, 1);

        isAuxIVArunning.set(true);

        int nItr = 30;
        double[] cFuncParam = new double[2];
        cFuncParam[0] = 1.0;
        cFuncParam[1] = 1.0;

        auxIVA = new AuxIVA(STFTdata, nItr, true, null, "norm", cFuncParam);

        auxIVAThread = new Thread(new AuxIVARunnable(), "AuxIVA Thread");
        auxIVAThread.start();

        istftButton.setEnabled(true);
    }

    void runISTFT() {


        isISTFTrunning.set(true);

        stft = new STFT();

        /*
        stft.istftm(STFTdata, 2048, 1024, "sine", audioDataLength);

        reSig = new double[nChannels][audioDataLength];
        Log.i("DEBUG", "Getting real signal");
        reSig = stft.getRealSigFromInvSTFT();
        Log.i("DEBUG", "Reconstructed sig: = " + Arrays.deepToString(reSig));

        testISTFT();

        */

        istftThread = new Thread(new ISTFTrunnable(), "ISTFT Thread");
        istftThread.start();
        writeToFileButton.setEnabled(true);
    }

    void saveChannelsToFiles() {
        for (int c = 0; c < nChannels; c++) {
            isWritingToFiles[c] = new AtomicBoolean(false);
            isWritingToFilesCompleted[c] = new AtomicBoolean(false);
            writingToFilesThread[c] = new Thread(new WriteToFileRunnable(c), "Write To File Thread: Channel " + (c + 1));
            writingToFilesThread[c].start();
        }
    }

    void getDataFromIntent() {
        Intent intent = getIntent();

        STFTParcel stftParcel = intent.getParcelableExtra(STFT_PARCEL);

        STFTdata = BSSMenuActivity_DevMode.STFTdata;
        nChannels = stftParcel.getnChannels();
        Log.i("DEBUG", "nChannels = " + nChannels);
        nFrames = stftParcel.getnFrames();
        Log.i("DEBUG", "nFrm = " + nFrames);
        nFreqs = stftParcel.getnFreqs();
        Log.i("DEBUG", "nFreq = " + nFreqs);
        winLen = stftParcel.getWinLen();
        Log.i("DEBUG", "winLen = " + winLen);
        nOverlap = stftParcel.getnOverlap();
        Log.i("DEBUG", "nOverlap = " + nOverlap);
        audioDataLength = stftParcel.getAudioDataLength();
        winFunc = stftParcel.getWinFunc();
    }

    private class AuxIVARunnable implements Runnable {

        @Override
        public void run() {
            //auxIVA.runDemixTest();



            auxIVA.run();
            //stopMethodTracing();

            STFTout = new Complex[nChannels][nFrames][nFreqs];
            STFTout = auxIVA.getSourceEstimatesSTFT();

            /*
            for (int s = 0; s < nChannels; s++) {
                for (int t = 0; t < nFrames; t++) {
                    Log.i("DEBUG", "STFTout[" + s + "][" + t + "] : " + Arrays.deepToString(STFTout[s][t]));
                }
            }
            */

            Log.i("DEBUG", "Separated STFT channels: " + STFTout.length);
            Log.i("DEBUG", "Separated STFT frames: " + STFTout[0].length);
            Log.i("DEBUG", "Separated STFT bins: " + STFTout[0][0].length);
            isAuxIVArunning.set(false);
            isAuxIVAcompleted.set(true);



            //testSTFT();
        }
    }

    private class ISTFTrunnable implements Runnable {
        @Override
        public void run() {
            stft.istftm(STFTout, winLen, nOverlap, winFunc, audioDataLength);

            reSig = new double[nChannels][audioDataLength];
            Log.i("DEBUG", "Getting real signal");
            reSig = stft.getRealSigFromInvSTFT();
            Log.i("DEBUG", "Reconstructed sig: = " + Arrays.deepToString(reSig));
            isISTFTrunning.set(false);
            isISTFTcompleted.set(true);
        }
    }

    private class WriteToFileRunnable implements Runnable {

        int chanIdx;

        WriteToFileRunnable(int channel) {
            this.chanIdx = channel;
        }

        @Override
        public void run() {

            short[] shortData = new short[audioDataLength];

            Log.i("DEBUG", "Converting data to shorts");

            Log.i("DEBUG", "Channel " + chanIdx);
            for (int t = 0; t < audioDataLength; t++) {
                shortData[t] = (short) (reSig[chanIdx][t] * 32768.0);
            }

            Log.i("DEBUG", "audioDataLength = " + audioDataLength);
            Log.i("DEBUG", "short length = " + shortData.length);

            Log.i("DEBUG", "Data converted");


            AudioFileWriter audioFileWriter = new AudioFileWriter(1, MainActivity_DevMode.SAMPLING_RATE_IN_HZ, 16);


            Log.i("DEBUG", "Writing chanIdx " + chanIdx);
            File file = new File(getExternalStorageDirectory().getAbsolutePath(), MainActivity_DevMode.fileNameNoExt + "AuxIVA_Channel" + (chanIdx + 1) + ".wav");
            try {
                audioFileWriter.convertShortArrayToFile(shortData, file);
            } catch (IOException e) {
                e.printStackTrace();
            }

            isWritingToFiles[chanIdx].set(false);
            isWritingToFilesCompleted[chanIdx].set(true);
            Log.i("DEBUG", "Demixed audio saved to file successfully");

        }
    }

    private boolean areAllChannelsWrittenToFiles() {
        boolean out = true;

        for (int c = 0; c < nChannels; c++) {
            out = (out && isWritingToFilesCompleted[c].get());
        }

        return out;
    }

    private boolean areSomeChannelsWritingToFiles() {
        boolean out = false;

        for (int c = 0; c < nChannels; c++) {
            out = (out || isWritingToFiles[c].get());
        }

        return out;
    }

    private void testISTFT() {

        Log.i("DEBUG", "Checking inverse against original value");

        double[][] absDiff = new double[nChannels][audioDataLength];
        double[][] ratioDiff = new double[nChannels][audioDataLength];

        for (int i = 0; i < nChannels; i++) {
            for (int j = 0; j < audioDataLength; j++) {
                absDiff[i][j] = abs(audioData[i][j] - reSig[i][j]);
                ratioDiff[i][j] = 100 * reSig[i][j] / audioData[i][j];

                Log.i("DEBUGTest",
                        "Channel = " + (i + 1) + ", sample = " + j
                                + ", x = " + audioData[i][j]
                                + ", x_hat = " + reSig[i][j]
                                + ", diff = " + absDiff[i][j]
                                + ", ratio = " + ratioDiff[i][j]);
            }
        }
    }

    private void testSTFT() {

        Log.i("DEBUG", "Checking inverse against original value");

        //int offcount = 0;

        Complex[][][] absDiff = new Complex[nChannels][nFrames][nFreqs];
        double[][][] ratioDiff = new double[nChannels][nFrames][nFreqs];

        for (int i = 0; i < nChannels; i++) {
            for (int j = 0; j < nFrames; j++) {
                for (int k = 0; k < nFreqs; k++) {
                    absDiff[i][j][k] = STFTout[i][j][k].subtract(STFTdata[i][j][k]);
                    ratioDiff[i][j][k] = 100 * STFTout[i][j][k].abs() / STFTdata[i][j][k].abs();

                /*
                if (abs(ratioDiff[i][j] - 100) > 0.1) {
                    offcount++;
                }
                */

                    Log.i("DEBUGTest",
                            "x = " + STFTdata[i][j][k]
                                    + ", x_hat = " + STFTout[i][j][k]
                                    + ", diff = " + absDiff[i][j][k]
                                    + ", ratio = " + ratioDiff[i][j][k]);
                }
            }

            //Log.i("DEBUGTest", "Off count = " + offcount);
        }
    }

    private void testSTFTtoOriginal() {

        Log.i("DEBUG", "Checking inverse against original value");

        //int offcount = 0;

        Complex[][][] absDiff = new Complex[nChannels][nFrames][nFreqs];
        double[][][] ratioDiff = new double[nChannels][nFrames][nFreqs];

        for (int i = 0; i < nChannels; i++) {
            for (int j = 0; j < nFrames; j++) {
                for (int k = 0; k < nFreqs; k++) {
                    absDiff[i][j][k] = STFTout[i][j][k].subtract(STFTActivity_DevMode.STFToutput[i][j][k]);
                    ratioDiff[i][j][k] = 100 * STFTout[i][j][k].abs() / STFTActivity_DevMode.STFToutput[i][j][k].abs();

                /*
                if (abs(ratioDiff[i][j] - 100) > 0.1) {
                    offcount++;
                }
                */

                    Log.i("DEBUGTest",
                            "x = " + STFTActivity_DevMode.STFToutput[i][j][k]
                                    + ", x_hat = " + STFTout[i][j][k]
                                    + ", diff = " + absDiff[i][j][k]
                                    + ", ratio = " + ratioDiff[i][j][k]);
                }
            }

            //Log.i("DEBUGTest", "Off count = " + offcount);
        }
    }

    public void launchNextActivity(View view) {

        if (areAllChannelsWrittenToFiles()) {
            for (int c = 0; c < nChannels; c++) {
                writingToFilesThread[c].interrupt();
                writingToFilesThread[c] = null;
            }
            //Intent intent = new Intent(this, BSSMenuActivity_UserMode.class);
            //intent.putExtra(STFT_PARCEL, stftParcel);
            //startActivity(intent);
        } else if (areSomeChannelsWritingToFiles()){
            Toast toast = Toast.makeText(getApplicationContext(), "Some channels are still being written to files", Toast.LENGTH_LONG);
            toast.show();
        } else {
            Toast toast = Toast.makeText(getApplicationContext(), "Separated channels have not been written to files", Toast.LENGTH_LONG);
            toast.show();
        }
    }
}
