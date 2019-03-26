package com.example.administrator.audiorecord.activities;

import android.Manifest;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.media.AudioAttributes;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.AudioTrack;
import android.media.MediaRecorder;
import android.os.Bundle;
import android.support.annotation.NonNull;
import android.support.design.widget.FloatingActionButton;
import android.support.design.widget.Snackbar;
import android.support.v4.app.ActivityCompat;
import android.support.v4.content.ContextCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Switch;

import com.example.administrator.audiorecord.R;
import com.example.administrator.audiorecord.audioprocessing.datahandler.AudioFileWriter;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicBoolean;

import static android.os.Environment.getExternalStorageDirectory;

public class MainActivity extends AppCompatActivity {

    private static final int SAMPLING_RATE_IN_HZ = 16000;   //44100 Hz is supported on all devices
    private static final int CHANNEL_CONFIG = AudioFormat.CHANNEL_IN_STEREO;
    private static final int AUDIO_FORMAT = AudioFormat.ENCODING_PCM_16BIT;
    private static final int BUFFER_SIZE_FACTOR = 2; // preemptively allocated space
    private static final int BUFFER_SIZE = AudioRecord.getMinBufferSize(SAMPLING_RATE_IN_HZ, CHANNEL_CONFIG, AUDIO_FORMAT) * BUFFER_SIZE_FACTOR;


    private int nChannels = 2;

    private AudioRecord recorder = null;
    private AudioTrack audioTrack = null;

    private Thread recordingThread = null;
    private Thread playbackThread = null;
    private Thread saveThread = null;

    private final AtomicBoolean isRecording = new AtomicBoolean(false);
    private final AtomicBoolean isRecordingCompleted = new AtomicBoolean(false);
    private final AtomicBoolean isPlayingBack = new AtomicBoolean(false);

    private final AtomicBoolean isTestMode = new AtomicBoolean(false);

    String fileName, fileNameNoExt;

    FloatingActionButton fabRecord, fabPlay, fabStop, fabSave;

    Snackbar saveBar;

    public static final String FILE_NAME = "com.example.administrator.FILE_NAME";
    public static final String FILE_NAME_NO_EXT = "com.example.administrator.FILE_NAME_NO_EXT";
    public static final String NUM_CHANNELS = "com.example.administrator.NUM_CHANNELS";
    public static final String SAMPLING_RATE = "com.example.administrator.SAMPLING_RATE";

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.usermodeactivity_main);

        requestRecordAudioPermission();
        requestWriteExternalPermission();
        requestInternetPermission();
        requestNetworkStatePermission();
        requestWakeLockPermission();

        Switch switchTestMode = findViewById(R.id.switchTestMode);

        switchTestMode.setOnCheckedChangeListener((buttonView, isChecked) -> {
            if (isChecked) {
                isTestMode.set(true);
                Snackbar.make(buttonView, "Test mode active", Snackbar.LENGTH_SHORT).show();
            } else {
                Snackbar.make(buttonView, "Test mode disabled", Snackbar.LENGTH_SHORT).show();
                isTestMode.set(false);
            }
        });

        fabRecord = findViewById(R.id.fabRecord);
        fabPlay = findViewById(R.id.fabPlay);
        fabStop = findViewById(R.id.fabStop);
        fabSave = findViewById(R.id.fabSave);

        fabPlay.setEnabled(false);
        fabStop.setEnabled(false);
        fabSave.setEnabled(false);

        fabRecord.setOnClickListener(view -> {
            if (isRecording.get()) {
                Snackbar.make(view, "Recording stopped", Snackbar.LENGTH_SHORT).show();
                fabRecord.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_mic_24px));
                stopRecording();
                fabPlay.setEnabled(true);
                fabSave.setEnabled(true);
            } else {
                Snackbar.make(view, "Recording", Snackbar.LENGTH_SHORT).show();
                fabRecord.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_stop_24px));
                startRecording();
            }
        });


        fabPlay.setOnClickListener(view -> {
            if (isPlayingBack.get()) {
                Snackbar.make(view, "Playback paused", Snackbar.LENGTH_SHORT).show();
                fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
                pausePlayback();
            } else {
                Snackbar.make(view, "Playing back", Snackbar.LENGTH_SHORT).show();
                fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_pause_24px));
                startPlayback();
                fabStop.setEnabled(true);
            }
        });

        fabStop.setOnClickListener(view -> {
            Snackbar.make(view, "Playback stopped", Snackbar.LENGTH_SHORT).show();
            fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
            stopPlayback();
            fabStop.setEnabled(false);
        });

        fabSave = findViewById(R.id.fabSave);
        fabSave.setOnClickListener(view -> {
            saveBar = Snackbar.make(view, "Saving file", Snackbar.LENGTH_SHORT);
            saveBar.show();
            saveWavFile();
        });
    }

    private void saveWavFile() {
        saveThread = new Thread(new SaveRunnable(), "Saving WAV File Thread");
        saveThread.run();
    }

    private void startRecording() {

        recorder = new AudioRecord(
                MediaRecorder.AudioSource.DEFAULT,
                SAMPLING_RATE_IN_HZ,
                CHANNEL_CONFIG,
                AUDIO_FORMAT,
                BUFFER_SIZE);

        nChannels = recorder.getChannelCount();

        recorder.startRecording();
        isRecording.set(true);
        recordingThread = new Thread(new RecordingRunnable(), "Recording Thread");
        recordingThread.start();
    }

    private void stopRecording() {
        if (null == recorder) {
            return;
        }

        isRecording.set(false);
        isRecordingCompleted.set(true);
        recorder.stop();
        recorder.release();
        recorder = null;
        recordingThread = null;

        Log.i("DEBUG", "Recording stopped");
    }

    private void startPlayback() {
        String filePath = getExternalStorageDirectory().getAbsolutePath() + File.separator + fileName;
        playbackThread = new Thread(new PlayRunnable(filePath), "Playback Thread");
        playbackThread.start();
    }

    private void pausePlayback() {
        audioTrack.pause();
    }

    private void stopPlayback() {

        if (audioTrack.getPlayState() != AudioTrack.PLAYSTATE_STOPPED) {
            audioTrack.stop();
        }

        playbackThread.interrupt();
        playbackThread = null;
    }

    private class RecordingRunnable implements Runnable {

        @Override
        public void run() {

            Log.i("DEBUG", "Now recording");

            SimpleDateFormat formatter = new SimpleDateFormat("yyMMddHHmmss", Locale.getDefault());
            Date dateTimeNow = new Date();
            fileNameNoExt = "recording" + formatter.format(dateTimeNow);
            fileName = fileNameNoExt + ".pcm";

            final File file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);
            final ByteBuffer buffer = ByteBuffer.allocateDirect(BUFFER_SIZE);

            try (final FileOutputStream outStream = new FileOutputStream(file)) {
                while (isRecording.get()) {
                    int result = recorder.read(buffer, BUFFER_SIZE);
                    if (result < 0) {
                        throw new RuntimeException("Reading of audio buffer failed: " +
                                getBufferReadFailureReason(result));
                    }
                    outStream.write(buffer.array(), 0, BUFFER_SIZE);
                    buffer.clear();
                }
            } catch (IOException e) {
                throw new RuntimeException("Writing of recorded audio failed", e);
            }
        }

        private String getBufferReadFailureReason(int errorCode) {
            switch (errorCode) {
                case AudioRecord.ERROR_INVALID_OPERATION:
                    return "ERROR_INVALID_OPERATION";
                case AudioRecord.ERROR_BAD_VALUE:
                    return "ERROR_BAD_VALUE";
                case AudioRecord.ERROR_DEAD_OBJECT:
                    return "ERROR_DEAD_OBJECT";
                case AudioRecord.ERROR:
                    return "ERROR";
                default:
                    return "Unknown (" + errorCode + ")";
            }
        }
    }

    private class PlayRunnable implements Runnable {

        String filePath;

        PlayRunnable(String filePath) {
            this.filePath = filePath;
        }

        @Override
        public void run() {
            Log.i("DEBUG", "Now playing");

            isPlayingBack.set(true);

            if (filePath == null) {
                Log.i("DEBUG", "filePath is null!");
                return;
            }


            int intSize = android.media.AudioTrack.getMinBufferSize(
                    SAMPLING_RATE_IN_HZ,
                    AudioFormat.CHANNEL_OUT_STEREO,
                    AudioFormat.ENCODING_PCM_16BIT);

            audioTrack = new AudioTrack.Builder()
                    .setAudioAttributes(new AudioAttributes.Builder()
                            .setUsage(AudioAttributes.USAGE_MEDIA)
                            .setContentType(AudioAttributes.CONTENT_TYPE_SPEECH)
                            .build())
                    .setAudioFormat(new AudioFormat.Builder()
                            .setEncoding(AudioFormat.ENCODING_PCM_16BIT)
                            .setSampleRate(SAMPLING_RATE_IN_HZ)
                            .setChannelMask(AudioFormat.CHANNEL_OUT_STEREO)
                            .build())
                    .setBufferSizeInBytes(intSize)
                    .build();

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
            fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
            runOnUiThread(() -> fabStop.setEnabled(false));
            Log.i("DEBUG", "Playback ended");
        }
    }

    private class SaveRunnable implements Runnable {
        @Override
        public void run() {
            File PCMfile = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);
            File Wavfile = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_stereo.wav");

            AudioFileWriter rtw = new AudioFileWriter(2, SAMPLING_RATE_IN_HZ, 16);
            try {
                rtw.convertPcmToWav(PCMfile, Wavfile);
            } catch (IOException e) {
                e.printStackTrace();
            }

            runOnUiThread(() -> saveBar.setText("Recording saved to file").show());
        }
    }

    private void requestRecordAudioPermission() {
        //check API version, do nothing if API version < 23!
        int currentapiVersion = android.os.Build.VERSION.SDK_INT;
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP) {

            if (ContextCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO) != PackageManager.PERMISSION_GRANTED) {

                // Should we show an explanation?
                if (ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.RECORD_AUDIO)) {

                    // Show an expanation to the user *asynchronously* -- don't block
                    // this thread waiting for the user's response! After the user
                    // sees the explanation, try again to request the permission.

                } else {

                    // No explanation needed, we can request the permission.

                    ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.RECORD_AUDIO}, 1);
                }
            }
        }
    }

    private void requestWriteExternalPermission() {
        //check API version, do nothing if API version < 23!
        int currentapiVersion = android.os.Build.VERSION.SDK_INT;
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP) {

            if (ContextCompat.checkSelfPermission(this, Manifest.permission.WRITE_EXTERNAL_STORAGE) != PackageManager.PERMISSION_GRANTED) {

                // Should we show an explanation?
                if (ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.WRITE_EXTERNAL_STORAGE)) {

                    // Show an expanation to the user *asynchronously* -- don't block
                    // this thread waiting for the user's response! After the user
                    // sees the explanation, try again to request the permission.

                } else {

                    // No explanation needed, we can request the permission.

                    ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.WRITE_EXTERNAL_STORAGE}, 1);
                }
            }
        }
    }

    private void requestInternetPermission() {
        //check API version, do nothing if API version < 23!
        int currentapiVersion = android.os.Build.VERSION.SDK_INT;
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP) {

            if (ContextCompat.checkSelfPermission(this, Manifest.permission.INTERNET) != PackageManager.PERMISSION_GRANTED) {

                // Should we show an explanation?
                if (ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.INTERNET)) {

                    // Show an expanation to the user *asynchronously* -- don't block
                    // this thread waiting for the user's response! After the user
                    // sees the explanation, try again to request the permission.

                } else {

                    // No explanation needed, we can request the permission.

                    ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.INTERNET}, 1);
                }
            }
        }
    }

    private void requestNetworkStatePermission() {
        //check API version, do nothing if API version < 23!
        int currentapiVersion = android.os.Build.VERSION.SDK_INT;
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP) {

            if (ContextCompat.checkSelfPermission(this, Manifest.permission.ACCESS_NETWORK_STATE) != PackageManager.PERMISSION_GRANTED) {

                // Should we show an explanation?
                if (ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.ACCESS_NETWORK_STATE)) {

                    // Show an expanation to the user *asynchronously* -- don't block
                    // this thread waiting for the user's response! After the user
                    // sees the explanation, try again to request the permission.

                } else {

                    // No explanation needed, we can request the permission.

                    ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.ACCESS_NETWORK_STATE}, 1);
                }
            }
        }
    }

    private void requestWakeLockPermission() {
        //check API version, do nothing if API version < 23!
        int currentapiVersion = android.os.Build.VERSION.SDK_INT;
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP) {

            if (ContextCompat.checkSelfPermission(this, Manifest.permission.WAKE_LOCK) != PackageManager.PERMISSION_GRANTED) {

                // Should we show an explanation?
                if (ActivityCompat.shouldShowRequestPermissionRationale(this, Manifest.permission.WAKE_LOCK)) {

                    // Show an expanation to the user *asynchronously* -- don't block
                    // this thread waiting for the user's response! After the user
                    // sees the explanation, try again to request the permission.

                } else {

                    // No explanation needed, we can request the permission.

                    ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.WAKE_LOCK}, 1);
                }
            }
        }
    }

    @Override
    public void onRequestPermissionsResult(int requestCode, @NonNull String permissions[], @NonNull int[] grantResults) {
        switch (requestCode) {
            case 1: {
                // If request is cancelled, the result arrays are empty.
                if (grantResults.length > 0 && grantResults[0] == PackageManager.PERMISSION_GRANTED) {

                    // permission was granted, yay! Do the
                    // contacts-related task you need to do.
                    Log.d("Activity", "Granted!");

                } else {

                    // permission denied, boo! Disable the
                    // functionality that depends on this permission.
                    Log.d("Activity", "Denied!");
                    finish();
                }
            }

            // other 'case' lines to check for other
            // permissions this app might request
        }
    }

    public void launchBssActivity(View view) {

        Log.i("DEBUG", "Launching next activity");

        Intent intent = new Intent(this, BSSActivity.class);

        if (isTestMode.get()) {

//            fileNameNoExt = "newTest";
//            fileNameNoExt = "dev1_male2_inst_mix";

//            List of test files:
//
//            3 sources:
//            fileNameNoExt = "dev1_male3_inst_mix";
            fileNameNoExt = "dev1_male3_liverec_130ms_1m_mix";
//            fileNameNoExt = "dev1_male3_liverec_130ms_5cm_mix";
//            fileNameNoExt = "dev1_male3_liverec_250ms_1m_mix";
//            fileNameNoExt = "dev1_male3_liverec_250ms_5cm_mix";
//            fileNameNoExt = "dev1_female3_inst_mix";
//            fileNameNoExt = "dev1_female3_liverec_130ms_1m_mix";
//            fileNameNoExt = "dev1_female3_liverec_130ms_5cm_mix";
//            fileNameNoExt = "dev1_female3_liverec_250ms_1m_mix";
//            fileNameNoExt = "dev1_female3_liverec_250ms_5cm_mix";
//
//            4 sources:
//            fileNameNoExt = "dev1_male4_inst_mix";
//            fileNameNoExt = "dev1_male4_liverec_130ms_1m_mix";
//            fileNameNoExt = "dev1_male4_liverec_130ms_5cm_mix";
//            fileNameNoExt = "dev1_male4_liverec_130ms_1m_mix";
//            fileNameNoExt = "dev1_male4_liverec_130ms_5cm_mix";
//            fileNameNoExt = "dev1_female4_inst_mix";
//            fileNameNoExt = "dev1_female4_liverec_130ms_1m_mix";
//            fileNameNoExt = "dev1_female4_liverec_130ms_5cm_mix";
//            fileNameNoExt = "dev1_female4_liverec_130ms_1m_mix";
//            fileNameNoExt = "dev1_female4_liverec_130ms_5cm_mix";

            fileName = fileNameNoExt + ".pcm";
        }

        intent.putExtra(FILE_NAME, fileName);
        intent.putExtra(FILE_NAME_NO_EXT, fileNameNoExt);
        intent.putExtra(NUM_CHANNELS, nChannels);
        intent.putExtra(SAMPLING_RATE, SAMPLING_RATE_IN_HZ);

        startActivity(intent);
    }
}
