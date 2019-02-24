package com.example.administrator.audiorecord.activities_devmode;

import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;
import android.widget.Toast;

import com.example.administrator.audiorecord.R;
import com.example.administrator.audiorecord.audioprocessing.commons.STFT;
import com.example.administrator.audiorecord.audioprocessing.datahandler.AudioFileWriter;
import com.example.administrator.audiorecord.audioprocessing.datahandler.STFTParcel;

import org.apache.commons.math3.complex.Complex;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.concurrent.atomic.AtomicBoolean;

import static android.os.Debug.startMethodTracing;
import static android.os.Debug.startMethodTracingSampling;
import static android.os.Debug.stopMethodTracing;
import static android.os.Environment.getExternalStorageDirectory;
import static java.lang.Math.abs;

public class STFTActivity_DevMode extends AppCompatActivity {

    int NUM_CHANNELS = 2;

    File file;
    short[] shortData;
    public static double[][] audioData;
    public static Complex[][][] STFToutput;
    double[][] reSig;

    int testLen;

    TextView textView;
    Button readButton;
    Button testButton;
    Button saveChan1Button, saveChan2Button, saveReconstructedSigButton;

    int audioDataLength, nFrames, nFreqs, winLen, nOverlap;
    String winFunc;

    String fileName, fileNameNoExt;

    public static final String STFT_PARCEL = "com.example.administrator.STFT_PARCEL";

    Thread pcmToSTFTthread;

    AtomicBoolean isDataFormatted = new AtomicBoolean(false);
    AtomicBoolean isSTFTrunning = new AtomicBoolean(false);
    AtomicBoolean isSTFTrunAlready = new AtomicBoolean(false);
    AtomicBoolean isISTFTrunAlready = new AtomicBoolean(false);

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.devmodeactivity_stft);

        Intent intent = getIntent();

        fileName = intent.getStringExtra(MainActivity_DevMode.FILE_NAME);
        fileNameNoExt = intent.getStringExtra(MainActivity_DevMode.FILE_NAME_NO_EXT);

        file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);
        Log.i("DEBUG", file.toString());

        textView = findViewById(R.id.textView);
        textView.setText(fileName);

        readButton = findViewById(R.id.bRead);
        readButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

                readButton.setEnabled(false);

                //Log.i("DEBUG", "Read file button clicked");
                isSTFTrunning.set(true);

                pcmToSTFTthread = new Thread(new PCMtoSTFTRunnable(), "PCM data to STFT Thread");
                pcmToSTFTthread.start();

                readButton.setEnabled(true);
            }
        });

        testButton = findViewById(R.id.bTest);
        testButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                runSTFTtest();
            }
        });


        saveChan1Button = findViewById(R.id.bSaveChan1);
        saveChan2Button = findViewById(R.id.bSaveChan2);
        saveReconstructedSigButton = findViewById(R.id.bSaveReconstructedSig);

        saveChan1Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (isDataFormatted.get()) {
                    Log.i("DEBUG", "Writing chanIdx 1");

                    short[] chan1 = new short[audioDataLength];

                    for (int i = 0; i < audioDataLength; i++) {
                        chan1[i] = (short) (audioData[0][i] * 32768.0);
                    }

                    File Wavfile1 = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_chan1mono.wav");

                    AudioFileWriter rtw1 = new AudioFileWriter(1, 16000, 16);
                    try {
                        rtw1.convertShortArrayToFile(chan1, Wavfile1);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    Log.i("DEBUG", "Writing chanIdx 1 - Done");
                } else {
                    Toast toast = Toast.makeText(getApplicationContext(), "Data has not been read to channels", Toast.LENGTH_LONG);
                    toast.show();
                }

            }
        });

        saveChan2Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (isDataFormatted.get()) {
                    Log.i("DEBUG", "Writing chanIdx 2");
                    short[] chan2 = new short[audioDataLength];

                    for (int i = 0; i < audioDataLength; i++) {
                        chan2[i] = (short) (audioData[1][i] * 32768.0);
                    }

                    File Wavfile2 = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_chan2mono.wav");

                    AudioFileWriter rtw1 = new AudioFileWriter(1, 16000, 16);
                    try {
                        rtw1.convertShortArrayToFile(chan2, Wavfile2);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    Log.i("DEBUG", "Writing chanIdx 2 - Done");
                } else {
                    Toast toast = Toast.makeText(getApplicationContext(), "Data has not been read to channels", Toast.LENGTH_LONG);
                    toast.show();
                }
            }
        });

        saveReconstructedSigButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (isISTFTrunAlready.get()) {
                    Log.i("DEBUG", "Writing reSig");
                    short[] recon = new short[audioDataLength * 2];

                    for (int i = 0; i < audioDataLength; i++) {
                        recon[2 * i] = (short) (reSig[0][i] * 32768.0);
                        recon[2 * i + 1] = (short) (reSig[1][i] * 32768.0);
                    }

                    File WavfileRe = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_reconStereo.wav");

                    AudioFileWriter rtw1 = new AudioFileWriter(2, 16000, 16);
                    try {
                        rtw1.convertShortArrayToFile(recon, WavfileRe);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    Log.i("DEBUG", "Writing reSig Done");
                } else {
                    Toast toast = Toast.makeText(getApplicationContext(), "Inverse STFT has not been run", Toast.LENGTH_LONG);
                    toast.show();
                }
            }
        });
    }

    private void readPCM() {
        Log.i("DEBUG", "Reading PCM to short array");

        int BufferSize;
        InputStream inputStream; //FileInputStream is a subclass of InputStream. Keep the highest type
        DataInputStream dataInputStream;

        try {
            if (file == null) {
                Log.e("DEBUG", "File pointer is null!");
                return;
            } else {
                Log.i("DEBUG", "File pointer is not null");
            }

            BufferSize = (int) (file.length() / 2); // max value of int in java is 2,147,483,647
            Log.i("DEBUG", String.valueOf(BufferSize));

            audioDataLength = BufferSize / 2;
            //Log.i("DEBUG", String.valueOf(audioDataLength));

            inputStream = new FileInputStream(file);
            //Log.i("DEBUG", "inputStream has been created.");

            dataInputStream = new DataInputStream(inputStream);
            //Log.i("DEBUG", "dataInputStream has been created.");

            shortData = new short[BufferSize];  // very impt to allocate memory to the array
            //Log.i("DEBUG", "Memory has been allocated for shortData array.");


            for (int i = 0; i < BufferSize; i++) {
                byte[] bytes = new byte[2];

                bytes[0] = dataInputStream.readByte();
                bytes[1] = dataInputStream.readByte();

                ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);

                short result = buffer.getShort();

                //Log.i("DEBUG", "readShort() returns " + result);

                shortData[i] = result;
                //Log.i("DEBUG", "and passed " + shortData[i] + " to shortData array");
            }


            Log.i("DEBUG", "Expected size of shortData = " + BufferSize);
            Log.i("DEBUG", "Actual size of shortData = " + shortData.length);


        } catch (FileNotFoundException e) {
            Log.e("File not found", "" + e);
        } catch (IOException e) {
            e.printStackTrace();
        }

        Log.i("DEBUG", "PCM successfully saved to short array");
    }

    private void rawDataToChannels() {

        Log.i("DEBUG", "Formatting short array to channels");

        audioData = new double[NUM_CHANNELS][audioDataLength];
        for (int c = 0; c < NUM_CHANNELS; c++) {
            for (int i = 0; i < audioDataLength; i++) {
                double value = shortData[NUM_CHANNELS * i + c] / 32768.0;   // Channels are interleaved
                audioData[c][i] = value;
            }
        }
        Log.i("DEBUG", "Formatting completed. PCM successfully saved to double array");
    }

    private void runSTFT() {
        Log.i("DEBUG", "Running STFT");

        STFT stft = new STFT();

        winLen = 2048;
        nOverlap = 1024;
        winFunc = "sine";

        reSig = new double[NUM_CHANNELS][audioDataLength];
        startMethodTracingSampling("STFT", 2147483647,1);
        stft.stftm(audioData, winLen, nOverlap, winFunc);
        stopMethodTracing();
        nFrames = stft.get_nFrames();
        nFreqs = stft.get_nFreqs();

        STFToutput = new Complex[NUM_CHANNELS][nFrames][nFreqs];
        STFToutput = stft.getSTFT();

        /*
        for (int s = 0; s < NUM_CHANNELS; s++) {
            for (int t = 0; t < nFrames; t++) {
                Log.i("DEBUG", "STFTout[" + s + "][" + t + "] : " + Arrays.deepToString(STFToutput[s][t]));
            }
        }
        */
        Log.i("DEBUG", "STFT channels: " + STFToutput.length);
        Log.i("DEBUG", "STFT frames: " + STFToutput[0].length);
        Log.i("DEBUG", "STFT bins: " + STFToutput[0][0].length);

        Log.i("DEBUG", "STFT successfully ran");
    }

    private void runSTFTtest() {

        Log.i("DEBUG", "Running STFT Test");

        STFT stft = new STFT();

        int winLen = 2048, nOverlap = 1024;
        String winFunc = "sine";

        testLen = audioDataLength;

        reSig = new double[NUM_CHANNELS][testLen];

        stft.stftm(audioData, winLen, nOverlap, winFunc);

        nFrames = stft.get_nFrames();
        nFreqs = stft.get_nFreqs();

        STFToutput = new Complex[NUM_CHANNELS][nFrames][nFreqs];

        STFToutput = stft.getSTFT();

        stft.istftm(STFToutput, winLen, nOverlap, winFunc, testLen);

        reSig = stft.getRealSigFromInvSTFT();

        isISTFTrunAlready.set(true);

        testSTFT();
    }

    private void testSTFT() {

        Log.i("DEBUG", "Checking inverse against original value");

        //int offcount = 0;

        double[][] absDiff = new double[NUM_CHANNELS][testLen];
        double[][] ratioDiff = new double[NUM_CHANNELS][testLen];

        for (int i = 0; i < NUM_CHANNELS; i++) {
            for (int j = 0; j < testLen; j++) {

                absDiff[i][j] = abs(audioData[i][j] - reSig[i][j]);
                ratioDiff[i][j] = 100 * reSig[i][j] / audioData[i][j];

                /*
                if (abs(ratioDiff[i][j] - 100) > 0.1) {
                    offcount++;
                }
                */

                Log.i("DEBUGTest",
                        "Channel = " + (i + 1) + ", sample = " + j
                                + ", x = " + audioData[i][j]
                                + ", x_hat = " + reSig[i][j]
                                + ", diff = " + absDiff[i][j]
                                + ", ratio = " + ratioDiff[i][j]);
            }
        }

        //Log.i("DEBUGTest", "Off count = " + offcount);
    }

    private class PCMtoSTFTRunnable implements Runnable {

        @Override
        public void run() {
            readPCM();
            Log.i("DEBUG", "PCM is read");

            rawDataToChannels();
            isDataFormatted.set(true);
            runSTFT();

            isSTFTrunning.set(false);
            isSTFTrunAlready.set(true);
        }
    }

    public void launchBSSMenuActivity(View view) {

        if (isSTFTrunAlready.get()) {
            Log.i("DEBUG", "Launching next activity");

            pcmToSTFTthread.interrupt();
            pcmToSTFTthread = null;

            STFTParcel stftParcel = new STFTParcel(NUM_CHANNELS, nFrames, nFreqs, audioDataLength, winLen, nOverlap, winFunc);

            Intent intent = new Intent(this, BSSMenuActivity_DevMode.class);

            intent.putExtra(STFT_PARCEL, stftParcel);

            startActivity(intent);
        } else if (isSTFTrunning.get()) {
            Toast toast = Toast.makeText(getApplicationContext(), "STFT is still running", Toast.LENGTH_LONG);
            toast.show();
        } else {
            Toast toast = Toast.makeText(getApplicationContext(), "STFT has not been run", Toast.LENGTH_LONG);
            toast.show();
        }
    }
}





