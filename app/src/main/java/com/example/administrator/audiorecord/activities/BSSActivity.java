package com.example.administrator.audiorecord.activities;

import android.content.Intent;
import android.media.AudioAttributes;
import android.media.AudioFormat;
import android.media.AudioTrack;
import android.os.Bundle;
import android.support.design.widget.FloatingActionButton;
import android.support.design.widget.Snackbar;
import android.support.v4.content.ContextCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.view.WindowManager;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Spinner;
import android.widget.Switch;
import android.widget.TextView;

import com.example.administrator.audiorecord.R;
import com.example.administrator.audiorecord.audioprocessing.bss.AuxIVA;
import com.example.administrator.audiorecord.audioprocessing.bss.DirectionalSCA;
import com.example.administrator.audiorecord.audioprocessing.bss.FastICA;
import com.example.administrator.audiorecord.audioprocessing.commons.STFT;
import com.example.administrator.audiorecord.audioprocessing.datahandler.AudioFileWriter;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.descriptive.rank.Max;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

import static android.os.Environment.getExternalStorageDirectory;
import static android.support.design.widget.Snackbar.LENGTH_SHORT;
import static org.apache.commons.math3.util.FastMath.abs;

public class BSSActivity extends AppCompatActivity implements AdapterView.OnItemSelectedListener {


    final int AuxIVAIndex = 0;
    final int SCAIndex = 1;
    final int FastICAIndex = 2;

    int NUM_CHANNELS_DEFAULT = 2, SAMPLING_RATE_DEFAULT = 16000, BIT_DEPTH = 16;
    int samplingRate;

    Spinner spinnerBSS, spinnerPlaybackSrc, spinnerBssSrc;
    TextView txtBssSrc;
    int bssType, playbackSrc;

    FloatingActionButton fabBSS, fabSave, fabPlaySrc, fabStopPlaySrc;
    Snackbar progressBar = null;
    View progressBarView = null;

    Thread bssThread = null, wavThread = null, playbackThread = null;

    String fileName, fileNameNoExt;

    SimpleDateFormat dateFormat;
    String datetimeSuffix;


    double averageTime;

    /* for STFT */

    STFT stft;

    File file;
    short[] shortData;
    double[][] audioData;
    Complex[][][] obsSTFT, demixedSTFT;
    double[][] demixedSig;

    int audioDataLength, nChannels, nFrames, nFreqs, winLen, nOverlap, nSrc;
    String winFunc;

    int testLen;

    /* for BSS */

    AuxIVA auxIVA = null;
    DirectionalSCA sca = null;
    FastICA fastIca = null;

    /* for data handling */

    AudioFileWriter audioFileWriter = null;
    File[] pcmFiles, wavFiles;
    File pcmMultichannel, wavMultichannel;
    String pcmMultichannelFileName;

    /* for playback */

    AudioTrack audioTrack;

    /* debugging */

    AtomicBoolean isSTFTtestMode = new AtomicBoolean(false);
    AtomicBoolean isBssRunning = new AtomicBoolean(false);
    AtomicBoolean isBssCompleted = new AtomicBoolean(false);
    AtomicBoolean isPlayingBack = new AtomicBoolean(false);

    public static final String FILE_NAME_NO_EXT = "com.example.administrator.FILE_NAME_NO_EXT";
    public static final String NUM_CHANNELS = "com.example.administrator.NUM_CHANNELS";
    public static final String SAMPLING_RATE = "com.example.administrator.SAMPLING_RATE";

    String bssString;

    String[] Mixture2Src, Mixture3Src, Mixture4Src;

    /*------------*/

    long startTime, endTime;
    long[] bssTime, stftTime, istftTime;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.usermodeactivity_bssmenu);

        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);

        Intent intent = getIntent();

        fileName = intent.getStringExtra(MainActivity.FILE_NAME);
        fileNameNoExt = intent.getStringExtra(MainActivity.FILE_NAME_NO_EXT);
        file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);

        nChannels = intent.getIntExtra(MainActivity.NUM_CHANNELS, NUM_CHANNELS_DEFAULT);
        Log.i("DEBUG", "nChannels = " + nChannels);

        nSrc = nChannels;

        samplingRate = intent.getIntExtra(MainActivity.SAMPLING_RATE, SAMPLING_RATE_DEFAULT);
        Log.i("DEBUG", "samplingRate = " + samplingRate);

        setBSSspinner();
        setBssSrcSpinnerInitial();
        txtBssSrc = findViewById(R.id.txtBssSrc);
        txtBssSrc.setVisibility(View.INVISIBLE);

        dateFormat = new SimpleDateFormat("yyMMddHHmmss", Locale.getDefault());

        prepFileNames();

        fabBSS = findViewById(R.id.fabBSS);
        fabBSS.setOnClickListener(v -> {
            if (isBssRunning.get()) {
                progressBar.setText("BSS is already running");
            } else {

                bssTime = new long[18];
                stftTime = new long[18];
                istftTime = new long[18];
                //datetimeSuffix = dateFormat.format(new Date());

                switch (bssType) {
                    case AuxIVAIndex: //AuxIVA

                        progressBar = Snackbar.make(v, "AuxIVA running", Snackbar.LENGTH_SHORT)
                                .setAction("Action", null);
                        progressBarView = progressBar.getView();
                        progressBar.show();

                        if (bssThread != null) {
                            bssThread.interrupt();
                            bssThread = null;
                        }

                        bssThread = new Thread(new AuxIvaThread(), "AuxIVA Thread");
                        bssThread.start();
                        isBssRunning.set(true);
                        fabSave.setEnabled(true);
                        break;
                    case SCAIndex: //SCA
                        progressBar = Snackbar.make(v, "SCA running", Snackbar.LENGTH_SHORT)
                                .setAction("Action", null);
                        progressBarView = progressBar.getView();
                        progressBar.show();

                        if (bssThread != null) {
                            bssThread.interrupt();
                            bssThread = null;
                        }

                        bssThread = new Thread(new ScaThread(), "SCA Thread");
                        bssThread.start();
                        isBssRunning.set(true);
                        fabSave.setEnabled(true);
                        break;
                    case FastICAIndex: //SCA
                        progressBar = Snackbar.make(v, "FastICA running", Snackbar.LENGTH_SHORT)
                                .setAction("Action", null);
                        progressBarView = progressBar.getView();
                        progressBar.show();

                        if (bssThread != null) {
                            bssThread.interrupt();
                            bssThread = null;
                        }

                        bssThread = new Thread(new IcaThread(), "FastICA Thread");
                        bssThread.start();
                        isBssRunning.set(true);
                        fabSave.setEnabled(true);
                        break;
                }

                setPlaybackSrcSpinner();
            }
        });


        fabSave = findViewById(R.id.fabSaveBSS);
        fabSave.setEnabled(false);
        fabSave.setOnClickListener(v -> {
            if (isBssCompleted.get()) {
                progressBar.setText("Saving audio files").show();
                wavThread = new Thread(new WavThread(), "PCM to Wav Thread");
                wavThread.start();
            } else {
                progressBar.setText("AuxIVA is not completed yet").show();
            }
        });

        Switch switchTestMode = findViewById(R.id.switchSTFT);

        switchTestMode.setOnCheckedChangeListener((buttonView, isChecked) -> {
            if (isChecked) {
                isSTFTtestMode.set(true);
            } else {
                isSTFTtestMode.set(false);
            }
        });

        fabPlaySrc = findViewById(R.id.fabPlaySrc);

        fabPlaySrc.setOnClickListener(v -> {
            if (isBssCompleted.get()) {
                if (isPlayingBack.get()) {
                    Snackbar.make(v, "Playback paused", LENGTH_SHORT)
                            .setAction("Action", null).show();
                    fabPlaySrc.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
                    pausePlayback();

                } else {
                    progressBar.setText("Playing back source " + (playbackSrc + 1)).show();
                    fabPlaySrc.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_pause_24px));

                    String filePath = pcmFiles[playbackSrc].getAbsolutePath();
                    playbackThread = new Thread(new PlayRunnable(filePath), "Playback Thread");
                    playbackThread.start();

                    isPlayingBack.set(true);

                    fabStopPlaySrc.setEnabled(true);
                }
            } else {
                progressBar.setText("AuxIVA is not completed yet").show();
            }
        });


        fabStopPlaySrc = findViewById(R.id.fabStopPlaySrc);
        fabStopPlaySrc.setEnabled(false);
        fabStopPlaySrc.setOnClickListener(v -> {
            stopPlayback();
        });

    }

    private void prepFileNames() {

        int[] dist = {50, 75, 100};

        Mixture2Src = new String[18];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                Mixture2Src[2 * j + i] = "mixture_otd_2src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
                Mixture2Src[2 * j + i + 6] = "mixture_ofc_2src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
                Mixture2Src[2 * j + i + 12] = "mixture_lth_2src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
            }
        }

        Mixture3Src = new String[18];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                Mixture3Src[2 * j + i] = "mixture_otd_3src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
                Mixture3Src[2 * j + i + 6] = "mixture_ofc_3src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
                Mixture3Src[2 * j + i + 12] = "mixture_lth_3src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
            }
        }

        Mixture4Src = new String[18];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                Mixture4Src[2 * j + i] = "mixture_otd_4src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
                Mixture4Src[2 * j + i + 6] = "mixture_ofc_4src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
                Mixture4Src[2 * j + i + 12] = "mixture_lth_4src_combi" + Integer.toString(i + 1) + "_" + dist[j] + "cm";
            }
        }
    }

    private void processTime(long[] timein, String type) {

        double[] doubleTime = Arrays.stream(timein).parallel().mapToDouble(i -> i / 1000.0).toArray();

        saveTimeLog(doubleTime, type);
    }

    void saveTimeLog(double[] print, String type) {
        String filename = "timeLog_" + type + "_" + dateFormat.format(new Date()) + ".csv";

        File file = new File(getExternalStorageDirectory().getAbsolutePath() + File.separator + filename);

        FileOutputStream outputStream;

        try {
            outputStream = new FileOutputStream(file);

            for (double i : print) {
                outputStream.write(Double.toString(i).getBytes());
                outputStream.write(", ".getBytes());
                outputStream.write(Objects.requireNonNull(System.getProperty("line.separator")).getBytes());
            }
            outputStream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /* Threads */
    private class AuxIvaThread implements Runnable {
        @Override
        public void run() {

            Log.i("DEBUG", "AuxIVA");

            if (nSrc != nChannels) {
                nSrc = nChannels;
            }

            for (int i = 0; i < 18; i++) {

                runOnUiThread(() -> {
                    progressBar.setText("Reading file").show();
                });

                Log.i("DEBUG", "Separating " + Mixture2Src[i]);
                fileName = Mixture2Src[i] + ".pcm";
                fileNameNoExt = Mixture2Src[i];
                file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);

                readPCM();
                rawDataToChannels();

                if (isSTFTtestMode.get()) {
                    runSTFTtest();
                } else {
                    datetimeSuffix = "";
                    auxIVA = null;

                    runOnUiThread(() -> {
                        progressBar.setText("Running STFT").show();
                    });

                    System.gc();
                    startTime = System.currentTimeMillis();
                    runSTFT(false);
                    endTime = System.currentTimeMillis();
                    System.gc();
                    stftTime[i] = endTime - startTime;

                    runOnUiThread(() -> {
                        progressBar.setDuration(Snackbar.LENGTH_INDEFINITE).setText("Separating sources").show();
                    });
                    System.gc();
                    startTime = System.currentTimeMillis();
                    runAuxIVA();
                    endTime = System.currentTimeMillis();
                    System.gc();
                    bssTime[i] = endTime - startTime;

                    runOnUiThread(() -> {
                        progressBar.setDuration(Snackbar.LENGTH_SHORT).setText("Reconstructing sources").show();
                    });

                    System.gc();
                    startTime = System.currentTimeMillis();
                    runISTFT(false);
                    endTime = System.currentTimeMillis();
                    System.gc();
                    istftTime[i] = endTime - startTime;

                    sourcesToPCM();
                    pcmToWav(); //for testing only

                }
            }

            processTime(stftTime, "stft");
            processTime(istftTime, "istft");
            processTime(bssTime, "auxiva");

            isBssRunning.set(false);
            isBssCompleted.set(true);

            runOnUiThread(() -> {
                progressBar.setText("AuxIVA completed").show();
            });

        }
    }

    private void runAuxIVA() {

        int nItr = 100;
        double[] cFuncParam = new double[2];
        cFuncParam[0] = 1.0;
        cFuncParam[1] = 1.0;

        demixedSTFT = new Complex[nChannels][nFrames][nFreqs];

        auxIVA = new AuxIVA(obsSTFT, nItr, true, null, "norm", cFuncParam);

        auxIVA.run();

        demixedSTFT = auxIVA.getSourceEstimatesSTFT();
    }

    private class ScaThread implements Runnable {
        @Override
        public void run() {

            for (int i = 0; i < 18; i++) {

                runOnUiThread(() -> {
                    progressBar.setText("Reading file").show();
                });

                Log.i("DEBUG", "Separating " + Mixture3Src[i]);
                fileName = Mixture3Src[i] + ".pcm";
                fileNameNoExt = Mixture3Src[i];
                file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);

                readPCM();
                rawDataToChannels();

                if (isSTFTtestMode.get()) {
                    runSTFTtest();
                } else {

                    datetimeSuffix = "";
                    sca = null;

                    runOnUiThread(() -> {
                        progressBar.setText("Running STFT").show();
                    });

                    System.gc();
                    startTime = System.currentTimeMillis();
                    runSTFT(false);
                    endTime = System.currentTimeMillis();
                    System.gc();
                    stftTime[i] = endTime - startTime;

                    runOnUiThread(() -> {
                        progressBar.setDuration(Snackbar.LENGTH_INDEFINITE).setText("Separating sources").show();
                    });

                    //startMethodTracingSampling("SCAnew2", 2147483647, 1);
                    System.gc();
                    startTime = System.currentTimeMillis();
                    runSCA();
                    endTime = System.currentTimeMillis();
                    System.gc();
                    bssTime[i] = endTime - startTime;
                    //stopMethodTracing();

                    Log.i("DEBUG", "SCA Completed");

                    runOnUiThread(() -> {
                        progressBar.setDuration(Snackbar.LENGTH_SHORT).setText("Reconstructing sources").show();
                    });

                    System.gc();
                    startTime = System.currentTimeMillis();
                    runISTFT(false);
                    endTime = System.currentTimeMillis();
                    System.gc();
                    istftTime[i] = endTime - startTime;

                    sourcesToPCM();
                    pcmToWav(); //for testing only
                }
            }

            processTime(stftTime, "stft");
            processTime(istftTime, "istft");
            processTime(bssTime, "sca");

            isBssRunning.set(false);
            isBssCompleted.set(true);

            runOnUiThread(() -> {
                progressBar.setText("SCA completed").show();
                getWindow().clearFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
            });
        }

    }

    private void runSCA() {
        int maxItr = 100;
        double eta = 10;
        boolean isDerivCheck = false;
        boolean isRowDecoupling = false;

        demixedSTFT = new Complex[nChannels][nFrames][nFreqs];

        /*boolean isDebugging = true;

        if (isDebugging) {

            //[nChannels][nFrames][nFreqs]
            Complex[][][] testSTFTarray = new Complex[2][5][1];

            double[][] realX = new double[][]{{-0.1022, 0.3192, -0.8649, -0.1649, 1.0933}, {-0.2414, 0.3129, -0.0301, 0.6277, 1.1093}};
            double[][] imagX = new double[][]{{-0.8637, -1.2141, -0.0068, -0.7697, -0.2256}, {0.0774, -1.1135, 1.5326, 0.3714, 1.1174}};

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 5; j++) {
                    testSTFTarray[i][j][0] = new Complex(realX[i][j], imagX[i][j]);
                }
            }

            sca = new DirectionalSCA(testSTFTarray, nItr, true, 3, 1e-2, true, false);
        } else {*/
        sca = new DirectionalSCA(obsSTFT, maxItr, nSrc, eta, isDerivCheck, isRowDecoupling);

        sca.run();

        demixedSTFT = sca.getSourceEstimatesSTFT();
    }

    private class IcaThread implements Runnable {
        @Override
        public void run() {

            for (int i = 0; i < 18; i++) {

                Log.i("DEBUG", "Separating " + Mixture2Src[i]);

                runOnUiThread(() -> {
                    progressBar.setText("Reading file").show();
                });

                fileName = Mixture2Src[i] + ".pcm";
                fileNameNoExt = Mixture2Src[i];
                file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);

                readPCM();
                rawDataToChannels();

                if (isSTFTtestMode.get()) {
                    runSTFTtest();
                } else {
                    datetimeSuffix = "";
                    fastIca = null;

                    runOnUiThread(() -> {
                        progressBar.setText("Running STFT").show();
                    });

                    System.gc();
                    startTime = System.currentTimeMillis();
                    runSTFT(false);
                    endTime = System.currentTimeMillis();
                    System.gc();
                    stftTime[i] = endTime - startTime;

                    runOnUiThread(() -> {
                        progressBar.setDuration(Snackbar.LENGTH_INDEFINITE).setText("Separating sources").show();
                    });

                    //startMethodTracingSampling("SCAnew2", 2147483647, 1);
                    startTime = System.currentTimeMillis();
                    runFastICA();
                    endTime = System.currentTimeMillis();
                    bssTime[i] = endTime - startTime;
                    //stopMethodTracing();

                    Log.i("DEBUG", "FastICA Completed");

                    runOnUiThread(() -> {
                        progressBar.setDuration(Snackbar.LENGTH_SHORT).setText("Reconstructing sources").show();
                    });

                    System.gc();
                    startTime = System.currentTimeMillis();
                    runISTFT(false);
                    endTime = System.currentTimeMillis();
                    System.gc();
                    istftTime[i] = endTime - startTime;
                    sourcesToPCM();
                    pcmToWav(); //for testing only


                }
            }

            processTime(stftTime, "stft");
            processTime(istftTime, "istft");
            processTime(bssTime, "ica");

            isBssRunning.set(false);
            isBssCompleted.set(true);

            runOnUiThread(() -> {
                progressBar.setText("SCA completed").show();
                getWindow().clearFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
            });

        }
    }

    private void runFastICA() {
        int maxItr = 100;

        demixedSTFT = new Complex[nChannels][nFrames][nFreqs];

        fastIca = new FastICA(obsSTFT, maxItr);

        fastIca.run();

        demixedSTFT = fastIca.getSourceEstimatesSTFT();
    }

    /* STFT */

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

            BufferSize = (int) (file.length() / 2);
            // divided by two since each short = 2 bytes
            // max value of int in java is 2,147,483,647
            Log.i("DEBUG", String.valueOf(BufferSize));

            audioDataLength = BufferSize / 2;
            // divided by two since the input is always stereo
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

//            Log.i("DEBUG", "Expected size of shortData = " + BufferSize);
//            Log.i("DEBUG", "Actual size of shortData = " + shortData.length);


        } catch (FileNotFoundException e) {
            Log.e("File not found", "" + e);
        } catch (IOException e) {
            e.printStackTrace();
        }

//        Log.i("DEBUG", "PCM successfully saved to short array");
    }

    private void rawDataToChannels() {

//        Log.i("DEBUG", "Formatting short array to channels");

        audioData = new double[nChannels][audioDataLength];
        for (int c = 0; c < nChannels; c++) {
            for (int i = 0; i < audioDataLength; i++) {
                double value = ((double) shortData[nChannels * i + c]) / 32768.0;   // Channels are interleaved
                audioData[c][i] = value;
            }
        }
//        Log.i("DEBUG", "Formatting completed. PCM successfully saved to double array");
    }

    private void runSTFT(boolean normalize) {
        Log.i("DEBUG", "Running STFT");

        stft = new STFT();

        winLen = 2048;
        nOverlap = 1024;
        winFunc = "sine";

//        Log.i("DEBUG", "audioData.length = " + audioData.length);
//        Log.i("DEBUG", "audioData[0].length = " + audioData[0].length);

        demixedSig = new double[nChannels][audioDataLength];
        //startMethodTracingSampling("STFT", 2147483647, 1);

        if (normalize) {
            IntStream.range(0, nChannels).parallel().forEach(c -> {

                double max = new Max().evaluate(audioData[c]);

                for (int t = 0; t < audioDataLength; t++) {
                    audioData[c][t] /= max;
                }
            });
        }

        stft.stftm(audioData, winLen, nOverlap, winFunc);
        //stopMethodTracing();
        nFrames = stft.get_nFrames();
        nFreqs = stft.get_nFreqs();

        obsSTFT = new Complex[nChannels][nFrames][nFreqs];
        obsSTFT = stft.getSTFT();

        /*
        for (int s = 0; s < NUM_CHANNELS; s++) {
            for (int t = 0; t < nFrames; t++) {
                Log.i("DEBUG", "STFTout[" + s + "][" + t + "] : " + Arrays.deepToString(obsSTFT[s][t]));
            }
        }
        */
//        Log.i("DEBUG", "STFT channels: " + obsSTFT.length);
//        Log.i("DEBUG", "STFT frames: " + obsSTFT[0].length);
//        Log.i("DEBUG", "STFT bins: " + obsSTFT[0][0].length);

//        Log.i("DEBUG", "STFT successfully ran");
    }

    private void runSTFTtest() {

        Log.i("DEBUG", "Running STFT Test");

        STFT stft = new STFT();

        int winLen = 2048, nOverlap = 1024;
        String winFunc = "sine";

        testLen = audioDataLength;

        demixedSig = new double[nChannels][testLen];

        stft.stftm(audioData, winLen, nOverlap, winFunc);

        nFrames = stft.get_nFrames();
        nFreqs = stft.get_nFreqs();

        obsSTFT = new Complex[nChannels][nFrames][nFreqs];

        obsSTFT = stft.getSTFT();

        stft.istftm(obsSTFT, winLen, nOverlap, winFunc, testLen);

        demixedSig = stft.getRealSigFromInvSTFT();

        testSTFT();
    }

    private void testSTFT() {

        Log.i("DEBUG", "Checking inverse against original value");

        //int offcount = 0;

        double[][] absDiff = new double[nChannels][testLen];
        double[][] ratioDiff = new double[nChannels][testLen];

        for (int i = 0; i < nChannels; i++) {
            for (int j = 0; j < testLen; j++) {

                absDiff[i][j] = abs(audioData[i][j] - demixedSig[i][j]);
                ratioDiff[i][j] = 100 * demixedSig[i][j] / audioData[i][j];

                /*
                if (abs(ratioDiff[i][j] - 100) > 0.1) {
                    offcount++;
                }
                */

                Log.i("DEBUGTest",
                        "Channel = " + (i + 1) + ", sample = " + j
                                + ", x = " + audioData[i][j]
                                + ", x_hat = " + demixedSig[i][j]
                                + ", diff = " + absDiff[i][j]
                                + ", ratio = " + ratioDiff[i][j]);
            }
        }

        //Log.i("DEBUGTest", "Off count = " + offcount);
    }

    private void runISTFT(boolean normalize) {
        stft.istftm(demixedSTFT, winLen, nOverlap, winFunc, audioDataLength);

        demixedSig = new double[nSrc][audioDataLength];

        demixedSig = stft.getRealSigFromInvSTFT();

        Log.i("DEBUG", "Real signal retrieved");

        if (normalize) {
            IntStream.range(0, nSrc).parallel().forEach(s -> {

                double max = new Max().evaluate(demixedSig[s]);

                for (int t = 0; t < audioDataLength; t++) {
                    demixedSig[s][t] /= max;
                }
            });
        }

        Log.i("DEBUG", "Normalizing");
    }

    /* Saving */

    private class WavThread implements Runnable {
        @Override
        public void run() {
            pcmToWav();
            runOnUiThread(() -> {
                progressBar.setText("Audio files saved").show();
            });

        }
    }

    private void sourcesToPCM() {

        pcmFiles = new File[nSrc];

        for (int s = 0; s < nSrc; s++) {
            pcmFiles[s] = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_" + bssString + "_Source" + (s + 1) + "_" + datetimeSuffix + ".pcm");
        }

        audioFileWriter = new AudioFileWriter(1, samplingRate, BIT_DEPTH);

        IntStream.range(0, nSrc).parallel().forEach(s -> {
            try {
                audioFileWriter.doubleArrayToPCM(demixedSig[s], pcmFiles[s], audioDataLength);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
    }

    private void pcmToWav() {

        wavFiles = new File[nSrc];

        for (int s = 0; s < nSrc; s++) {
            wavFiles[s] = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_" + bssString + "_Source" + (s + 1) + "_" + datetimeSuffix + ".wav");
        }

        pcmMultichannelFileName = fileNameNoExt + "_" + bssString + "_AllSources_" + datetimeSuffix + ".pcm";
        pcmMultichannel = new File(getExternalStorageDirectory().getAbsolutePath(), pcmMultichannelFileName);
        wavMultichannel = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_" + bssString + "_AllSources_" + datetimeSuffix + ".wav");

        if (audioFileWriter == null) {
            audioFileWriter = new AudioFileWriter(1, samplingRate, BIT_DEPTH);
        }

        IntStream.range(0, nSrc).parallel().forEach(s -> {
            try {
                audioFileWriter.convertPcmToWav(pcmFiles[s], wavFiles[s]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

        audioFileWriter = new AudioFileWriter(nSrc, samplingRate, BIT_DEPTH);

        try {
            audioFileWriter.doubleMultichannelArrayToPCM(demixedSig, pcmMultichannel, audioDataLength);
            audioFileWriter.convertPcmToWav(pcmMultichannel, wavMultichannel);
        } catch (IOException e) {
            e.printStackTrace();
        }


        Log.i("DEBUG", "Saving completed.");
    }

    /* Playback */

    private class PlayRunnable implements Runnable {

        String filePath;

        PlayRunnable(String filePath) {
            this.filePath = filePath;
        }

        @Override
        public void run() {
            Log.i("DEBUG", "Now playing");

            if (filePath == null) {
                Log.i("DEBUG", "filePath is null!");
                return;
            }


            int bufferSize = android.media.AudioTrack.getMinBufferSize(
                    samplingRate,
                    AudioFormat.CHANNEL_OUT_MONO,
                    AudioFormat.ENCODING_PCM_16BIT);

            Log.i("DEBUG", "buffer size = " + bufferSize);

            audioTrack = new AudioTrack.Builder()
                    .setAudioAttributes(new AudioAttributes.Builder()
                            .setUsage(AudioAttributes.USAGE_MEDIA)
                            .setContentType(AudioAttributes.CONTENT_TYPE_SPEECH)
                            .build())
                    .setAudioFormat(new AudioFormat.Builder()
                            .setEncoding(AudioFormat.ENCODING_PCM_16BIT)
                            .setSampleRate(samplingRate)
                            .setChannelMask(AudioFormat.CHANNEL_OUT_MONO)
                            .build())
                    .setBufferSizeInBytes(bufferSize)
                    .build();

            Log.i("DEBUG", "Sampling rate = " + audioTrack.getSampleRate());

            int count = 512 * 1024; // 512 kb

            byte[] byteData;
            File file;
            file = new File(filePath);

            byteData = new byte[count];
            FileInputStream in = null;
            try {
                in = new FileInputStream(file);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            int bytesread = 0, ret = 0;
            int size = (int) file.length();
            audioTrack.play();
            while (bytesread < size) {
                try {
                    assert in != null;
                    ret = in.read(byteData, 0, count);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                if (ret != -1) { // Write the byte array to the track
                    audioTrack.write(byteData, 0, ret);
                    bytesread += ret;
                } else break;
            }
            try {
                assert in != null;
                in.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            audioTrack.pause();
            audioTrack.flush();
            audioTrack.release();

            isPlayingBack.set(false);
            fabPlaySrc.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));

            runOnUiThread(() -> fabStopPlaySrc.setEnabled(false));
            Log.i("DEBUG", "Playback ended");
        }

    }

    private void pausePlayback() {
        audioTrack.pause();

        isPlayingBack.set(false);
    }

    private void stopPlayback() {
        audioTrack.pause();
        audioTrack.flush();
        audioTrack.release();
        playbackThread.interrupt();
        playbackThread = null;

        isPlayingBack.set(false);
    }

    /* Spinner */

    void setPlaybackSrcSpinner() {
        spinnerPlaybackSrc = findViewById(R.id.spinnerPlaybackSrc);
        spinnerPlaybackSrc.setOnItemSelectedListener(this);

        List<Integer> sourceArray = new ArrayList<>(nSrc);
        for (int s = 0; s < nSrc; s++) {
            sourceArray.add(s + 1);
        }

        ArrayAdapter<Integer> adapter = new ArrayAdapter<>(this, android.R.layout.simple_spinner_item, sourceArray);

        adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerPlaybackSrc.setAdapter(adapter);

    }

    void setBSSspinner() {
        spinnerBSS = findViewById(R.id.spinnerBSS);
        spinnerBSS.setOnItemSelectedListener(this);

        ArrayAdapter<CharSequence> adapter = ArrayAdapter.createFromResource(this,
                R.array.bss_options, android.R.layout.simple_spinner_item);
        adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerBSS.setAdapter(adapter);
    }

    void setBssSrcSpinnerInitial() {
        spinnerBssSrc = findViewById(R.id.spinnerBssSrc);
        spinnerBssSrc.setOnItemSelectedListener(this);
        spinnerBssSrc.setVisibility(View.INVISIBLE);
    }

    void setBssSrcSpinnerVisible() {
        int capacity = 3;
        List<Integer> sourceArray = new ArrayList<>(capacity);
        for (int i = 0; i < capacity; i++) {
            sourceArray.add(nChannels + i);
        }

        ArrayAdapter<Integer> adapter = new ArrayAdapter<>(this, android.R.layout.simple_spinner_item, sourceArray);

        adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerBssSrc.setAdapter(adapter);
        spinnerBssSrc.setVisibility(View.VISIBLE);
        txtBssSrc.setVisibility(View.VISIBLE);
    }

    void setBssSrcSpinnerInvisible() {
        spinnerBssSrc.setVisibility(View.INVISIBLE);
        txtBssSrc.setVisibility(View.INVISIBLE);
        nSrc = nChannels;
    }

    @Override
    public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
        switch (parent.getId()) {
            case R.id.spinnerBSS:
                bssType = position;

                bssString = parent.getItemAtPosition(position).toString();

                Snackbar.make(view, parent.getItemAtPosition(position).toString() + " selected", LENGTH_SHORT)
                        .setAction("Action", null).show();

                switch (bssType) {
                    case SCAIndex:
                        setBssSrcSpinnerVisible();
                        break;
                    default:
                        setBssSrcSpinnerInvisible();
                        break;
                }

                break;
            case R.id.spinnerPlaybackSrc:
                playbackSrc = position;
                if (isPlayingBack.get()) {
                    stopPlayback();
                    fabPlaySrc.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
                }

                if (isBssCompleted.get()) {
                    Snackbar.make(view, "Source " + parent.getItemAtPosition(position).toString() + " selected", LENGTH_SHORT)
                            .setAction("Action", null).show();
                }

                break;
            case R.id.spinnerBssSrc:
                if (!isBssRunning.get()) {
                    nSrc = Integer.parseInt(parent.getItemAtPosition(position).toString());
                    Snackbar.make(view, "Number of sources is " + parent.getItemAtPosition(position).toString() + ".", LENGTH_SHORT)
                            .setAction("Action", null).show();
                } else {
                    Snackbar.make(view, "Separation has already started.", LENGTH_SHORT)
                            .setAction("Action", null).show();
                }

                break;
        }
    }

    @Override
    public void onNothingSelected(AdapterView<?> parent) {

    }

    public void launchAsrActivity(View view) {

        Log.i("DEBUG", "Launching next activity");

        Intent intent = new Intent(this, ASRActivity.class);

        intent.putExtra(FILE_NAME_NO_EXT, fileNameNoExt);
        intent.putExtra(NUM_CHANNELS, nChannels);
        intent.putExtra(SAMPLING_RATE, samplingRate);

        startActivity(intent);
    }
}
