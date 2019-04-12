package com.example.administrator.androidbss.activities;

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

import com.example.administrator.androidbss.R;
import com.example.administrator.androidbss.audioprocessing.datahandler.AudioFileWriter;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
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
    private final AtomicBoolean isSaved = new AtomicBoolean(false);

    private final AtomicBoolean isTestMode = new AtomicBoolean(false);

    String fileName = null, fileNameNoExt;

    FloatingActionButton fabRecord, fabPlay, fabStop, fabSave;

    Snackbar saveBar;

    public static final String FILE_NAME = "com.example.administrator.FILE_NAME";
    public static final String FILE_NAME_NO_EXT = "com.example.administrator.FILE_NAME_NO_EXT";
    public static final String NUM_CHANNELS = "com.example.administrator.NUM_SOURCES";
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

//        switchTestMode.setVisibility(View.INVISIBLE);

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
            if(audioTrack == null || audioTrack.getPlayState() == AudioTrack.PLAYSTATE_STOPPED){
                Snackbar.make(view, "Playing back", Snackbar.LENGTH_SHORT).show();
                fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_pause_24px));
                startPlayback();
                fabStop.setEnabled(true);
            } else if (audioTrack.getPlayState() == AudioTrack.PLAYSTATE_PAUSED){
                Snackbar.make(view, "Playing back", Snackbar.LENGTH_SHORT).show();
                fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_pause_24px));
                resumePlayback();
                fabStop.setEnabled(true);
            } else if (audioTrack.getPlayState() == AudioTrack.PLAYSTATE_PLAYING){
                Snackbar.make(view, "Playback paused", Snackbar.LENGTH_SHORT).show();
                fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
                pausePlayback();
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
        isSaved.set(false);

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

    private void resumePlayback() {
        audioTrack.play();
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

            File file;
            file = new File(filePath);

            int intSize = (int) file.length();

                    android.media.AudioTrack.getMinBufferSize(
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
                    .setTransferMode(AudioTrack.MODE_STATIC)
                    .setBufferSizeInBytes(intSize)
                    .build();

            int count = 512 * 1024; // 512 kb

            try {
                if (file == null) {
                    Log.e("DEBUG", "File pointer is null!");
                    return;
                } else {
                    Log.i("DEBUG", "File pointer is not null");
                }

                int sizeInShorts = (int) (file.length() / 2);
                FileInputStream inputStream = new FileInputStream(file);
                //Log.i("DEBUG", "inputStream has been created.");

                DataInputStream dataInputStream = new DataInputStream(inputStream);
                //Log.i("DEBUG", "dataInputStream has been created.");

                short[] shortData = new short[sizeInShorts];  // very impt to allocate memory to the array
                //Log.i("DEBUG", "Memory has been allocated for shortData array.");

                for (int i = 0; i < sizeInShorts; i++) {
                    byte[] bytes = new byte[2];
                    bytes[0] = dataInputStream.readByte();
                    bytes[1] = dataInputStream.readByte();
                    ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
                    short result = buffer.getShort();
                    shortData[i] = result;
                }

                audioTrack.write(shortData, 0, sizeInShorts);
                audioTrack.setNotificationMarkerPosition(sizeInShorts/2);

            } catch (FileNotFoundException e) {
                Log.e("File not found", "" + e);
            } catch (IOException e) {
                e.printStackTrace();
            }


            audioTrack.setPlaybackPositionUpdateListener(new AudioTrack.OnPlaybackPositionUpdateListener(){

                @Override
                public void onMarkerReached(AudioTrack arg0) {
                    runOnUiThread(() -> {
                        fabPlay.setImageDrawable(ContextCompat.getDrawable(getApplicationContext(), R.drawable.ic_round_play_arrow_24px));
                        fabStop.setEnabled(false);
                    });
                }

                @Override
                public void onPeriodicNotification(AudioTrack arg0) {

                }

            });

            Log.i("DEBUG", "time before = " + System.currentTimeMillis());
            audioTrack.play();
            Log.i("DEBUG", "time after = " + System.currentTimeMillis());

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

            isSaved.set(true);

            //runOnUiThread(() -> saveBar.setText("Recording saved to file").show());
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

        if(!isSaved.get() && !isTestMode.get()){
            saveThread = new Thread(new SaveRunnable(), "save thread");
            saveThread.start();
        }


        if(fileName == null){
            Snackbar.make(findViewById(R.id.fabRecord), "No recording has been made yet", Snackbar.LENGTH_SHORT).show();
        }



        Log.i("DEBUG", "Launching next activity");

        Intent intent = new Intent(this, BSSActivity.class);

        if (isTestMode.get()) {

            fileNameNoExt = "mixture_otd_2src_config2";
            fileName = fileNameNoExt + ".pcm";
        }

        intent.putExtra(FILE_NAME, fileName);
        intent.putExtra(FILE_NAME_NO_EXT, fileNameNoExt);
        intent.putExtra(NUM_CHANNELS, nChannels);
        intent.putExtra(SAMPLING_RATE, SAMPLING_RATE_IN_HZ);

        startActivity(intent);
    }
}
