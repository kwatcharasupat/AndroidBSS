package com.example.administrator.audiorecord;

import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import org.apache.commons.math3.complex.Complex;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import static android.os.Environment.getExternalStorageDirectory;
import static java.lang.Math.abs;


// package com.google.corp.productivity.specialprojects.android.fft

public class STFTActivityComplex extends AppCompatActivity {

    int NUM_CHANNELS = 2;

    File file;
    short[] shortData;
    double[][] audioData;
    Complex[][][] STFToutput;
    double[][] reSig;

    int testLen;

    TextView textView;
    Button readButton;
    Button testButton;

    int audioDataLength;

    double[][] paddedInv;

    private Thread stftThread = null;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_stftcomplex);

        Intent intent = getIntent();

        String fileName = intent.getStringExtra(MainActivity.FILE_NAME);
        file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);

        textView = findViewById(R.id.textView);
        textView.setText(fileName);

        readButton = findViewById(R.id.bRead);
        readButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

                readButton.setEnabled(false);

                //Log.i("DEBUG", "Read file button clicked");
                readPCM();

                rawDataToChannels();

                runSTFT();

                readButton.setEnabled(true);
            }
        });

        testButton = findViewById(R.id.bTest);
        testButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                //runFFTtest();
                runSTFTtest();
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

            audioDataLength = BufferSize / 2;

            inputStream = new FileInputStream(file);
            //Log.i("DEBUG", "inputStream has been created.");

            dataInputStream = new DataInputStream(inputStream);
            //Log.i("DEBUG", "dataInputStream has been created.");

            shortData = new short[BufferSize];  // very impt to allocate memory to the array
            //Log.i("DEBUG", "Memory has been allocated for shortData array.");


            for (int i = 0; i < BufferSize; i++) {
                short result = dataInputStream.readShort();
                //Log.i("DEBUG", "readShort() returns " + result);

                shortData[i] = result;
                //Log.i("DEBUG", "and passed " + shortData[i] + " to shortData array");
            }


            //Log.i("DEBUG", "Expected size of shortData = " + BufferSize);
            //Log.i("DEBUG", "Actual size of shortData = " + shortData.length);


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

        for (int j = 0; j < NUM_CHANNELS; j++) {
            for (int i = 0; i < audioDataLength; i++) {

                double value = shortData[NUM_CHANNELS * i + j] / 32768.0;

                audioData[j][i] = value;
            }
        }


        //Log.d("DEBUG", "Expected size of audioData = " + audioDataLength);
        //Log.d("DEBUG", "Actual size of audioData channel 1 = " + audioData[0].length);
        //Log.d("DEBUG", "Actual size of audioData channel 2 = " + audioData[0].length);

        Log.i("DEBUG", "Formatting completed. PCM successfully saved to double array");
    }

    private void runSTFT() {
        stftThread = new Thread(new STFTRunnable(), "STFT Thread");
    }

    private void runFFTtest() {

        STFT stft = new STFT();
        testLen = 16;

        double[] x = new double[testLen];
        double[] y = new double[testLen];

        Arrays.fill(x, 1);

        stft.fft(x, y);

        Log.i("DEBUG", "x_stft_re: " + Arrays.toString(x));
        Log.i("DEBUG", "x_stft_im: " + Arrays.toString(y));

        Arrays.fill(x, 1);
        Arrays.fill(y, 0);

        stft.ifft(x, y);

        Log.i("DEBUG", "x_hat: " + Arrays.toString(x));
        Log.i("DEBUG", "y_hat: " + Arrays.toString(y));


    }

    private void runSTFTtest() {

        Log.i("DEBUG", "Running STFT Test");

        STFTComplex stft = new STFTComplex();

        int winLen = 2048, nOverlap = 1024;
        String winFunc = "sine";

        testLen = audioDataLength;

        reSig = new double[NUM_CHANNELS][testLen];

        stft.stftm(audioData, winLen, nOverlap, winFunc);

        STFToutput = new Complex[NUM_CHANNELS][stft.get_nFrames()][stft.get_nFreq()];

        STFToutput = stft.getSTFT();

        stft.istftm(STFToutput, winLen, nOverlap, winFunc, testLen);

        reSig = stft.getRealSigFromInvSTFT();

        /*
        Log.i("DEBUG", "Displaying reSig");
        for (int c = 0; c < NUM_CHANNELS; c++) {
                Log.i("DEBUG", Arrays.toString(reSig[c]));
        }
        */

        //Log.i("DEBUG", "Inverse STFT ran");

        /*
        paddedInv = new double[NUM_CHANNELS][stft.getPaddedLength()];
        paddedInv = stft.getInvPaddedOutput();

        for (int c = 0; c < NUM_CHANNELS; c++){
            for(int i = 0; i < stft.getPaddedLength(); i++){
                if(paddedInv[c][i] < 1e-3){
                    Log.i("DEBUG", "paddedInv[" + c + "][" + i + "] = 0");
                }
            }
        }
        */

        testSTFT();
        Log.i("DEBUG", "STFT Test completed");
    }

    private class STFTRunnable implements Runnable {

        @Override
        public void run() {

            Log.i("DEBUG", "Running STFT");

            STFTComplex stft = new STFTComplex();

            int winLen = 2048, nOverlap = 1024;
            String winFunc = "sine";

            reSig = new double[NUM_CHANNELS][audioDataLength];

            stft.stftm(audioData, winLen, nOverlap, winFunc);

            STFToutput = new Complex[NUM_CHANNELS][stft.get_nFrames()][stft.get_nFreq()];
            STFToutput = stft.getSTFT();

            Log.i("DEBUG", "STFT successfully ran");
        }
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
}





