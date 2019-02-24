
package com.example.administrator.audiorecord.activities_devmode;

import android.Manifest;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioRecord;
import android.media.AudioTrack;
import android.media.MediaRecorder;
import android.os.Bundle;
import android.support.v4.app.ActivityCompat;
import android.support.v4.content.ContextCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.Toast;

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
import java.util.concurrent.atomic.AtomicBoolean;

import static android.os.Environment.getExternalStorageDirectory;

public class MainActivity_DevMode extends AppCompatActivity {

    public static final int SAMPLING_RATE_IN_HZ = 16000;   //44100 Hz is supported on all devices

    private static final int CHANNEL_CONFIG = AudioFormat.CHANNEL_IN_STEREO;

    private static final int AUDIO_FORMAT = AudioFormat.ENCODING_PCM_16BIT;

    private static final int BUFFER_SIZE_FACTOR = 2; // preemptively allocated space

    private static final int BUFFER_SIZE = AudioRecord.getMinBufferSize(
            SAMPLING_RATE_IN_HZ,
            CHANNEL_CONFIG,
            AUDIO_FORMAT) * BUFFER_SIZE_FACTOR;

    private final AtomicBoolean isRecording = new AtomicBoolean(false);
    private final AtomicBoolean isRecordingCompleted = new AtomicBoolean(false);
    private final AtomicBoolean isPlayingBack = new AtomicBoolean(false);

    private AudioRecord recorder = null;
    private Thread recordingThread = null;
    private Thread playThread = null;

    private Button startRecordButton;
    private Button stopRecordButton;
    private Button startPlayButton;
    private Button stopPlayButton;
    private Button pcmToWavButton;

    public static final String FILE_NAME = "com.example.administrator.FILE_NAME";
    public static final String FILE_NAME_NO_EXT = "com.example.administrator.FILE_NAME_NO_EXT";
    public static String fileName, fileNameNoExt;

    AudioTrack audioTrack;

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.devmodeactivity_main);



        requestRecordAudioPermission();
        requestWriteExternalPermission();

        Log.i("DEBUG", getExternalStorageDirectory().getAbsolutePath());

        startRecordButton = findViewById(R.id.bStartRecord);
        startRecordButton.setOnClickListener(v -> {
            startRecording();
            startRecordButton.setEnabled(false);
            stopRecordButton.setEnabled(true);
            startPlayButton.setEnabled(false);
        });

        stopRecordButton = findViewById(R.id.bStopRecord);
        stopRecordButton.setEnabled(false);
        stopRecordButton.setOnClickListener(v -> {
            stopRecording();
            startRecordButton.setEnabled(true);
            stopRecordButton.setEnabled(false);
            startPlayButton.setEnabled(true);
        });

        startPlayButton = findViewById(R.id.bPlayRecord);
        startPlayButton.setEnabled(false);
        startPlayButton.setOnClickListener(v -> {
            if (isPlayingBack.get()) {
                Toast toast = Toast.makeText(getApplicationContext(), "Current recording is still playing", Toast.LENGTH_LONG);
                toast.show();
            } else if (!isRecordingCompleted.get()) {
                Toast toast = Toast.makeText(getApplicationContext(), "No recording has been made", Toast.LENGTH_LONG);
                toast.show();
            } else {
                Log.i("DEBUG", "Creating play Thread");
                playRecording();
                stopPlayButton.setEnabled(true);
            }
            stopRecordButton.setEnabled(false);
            startRecordButton.setEnabled(true);
        });

        stopPlayButton = findViewById(R.id.bStopPlay);
        stopPlayButton.setEnabled(false);
        stopPlayButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (isPlayingBack.get()) {
                    audioTrack.pause();
                    audioTrack.flush();
                    audioTrack.release();
                    playThread.interrupt();
                    playThread = null;
                    stopPlayButton.setEnabled(false);
                }
            }
        });

        pcmToWavButton = findViewById(R.id.bPCM2WAV);
        pcmToWavButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (isRecordingCompleted.get()) {
                    Log.i("DEBUG", "Preparing to save PCM to WAV");
                    File PCMfile = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);
                    File Wavfile = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_stereo.wav");

                    Log.i("DEBUG", "Saving");
                    AudioFileWriter rtw = new AudioFileWriter(2, SAMPLING_RATE_IN_HZ, 16);
                    try {
                        rtw.convertPcmToWav(PCMfile, Wavfile);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    Log.i("DEBUG", "Done");
                } else {
                    Toast toast = Toast.makeText(getApplicationContext(), "No recording has been completed", Toast.LENGTH_LONG);
                    toast.show();
                }

            }
        });
    }

    @Override
    protected void onResume() {
        super.onResume();

        startRecordButton.setEnabled(true);
        stopRecordButton.setEnabled(false);
    }

    @Override
    protected void onPause() {
        super.onPause();
        stopRecording();
    }

    private void startRecording() {

        recorder = new AudioRecord(
                MediaRecorder.AudioSource.DEFAULT,
                SAMPLING_RATE_IN_HZ,
                CHANNEL_CONFIG,
                AUDIO_FORMAT,
                BUFFER_SIZE);

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

    private class RecordingRunnable implements Runnable {

        @Override
        public void run() {

            Log.i("DEBUG", "Now recording");

            SimpleDateFormat formatter = new SimpleDateFormat("yyMMddHHmmss");
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

        public PlayRunnable(String filePath) {
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
                    AudioFormat.CHANNEL_IN_STEREO,
                    AudioFormat.ENCODING_PCM_16BIT);

            audioTrack = new AudioTrack(
                    AudioManager.STREAM_MUSIC,
                    SAMPLING_RATE_IN_HZ,
                    AudioFormat.CHANNEL_IN_STEREO,
                    AudioFormat.ENCODING_PCM_16BIT,
                    intSize,
                    AudioTrack.MODE_STREAM);

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

            Log.i("DEBUG", "Playback ended");
        }
    }
    private void playRecording(/*String filePath*/) {
        String filePath = getExternalStorageDirectory().getAbsolutePath() + "/" + fileName;
        playThread = new Thread(new PlayRunnable(filePath), "Playback Thread");
        playThread.start();
    }

    public void launchSTFTactivityTestMode(View view) {

        Log.i("DEBUG", "Launching next activity");

        Intent intent = new Intent(this, STFTActivity_DevMode.class);

        //fileNameNoExt = "recording190202004045";
        //fileName = fileNameNoExt + ".pcm";

        fileNameNoExt = "testAudio";
        fileName = fileNameNoExt + ".pcm";

        intent.putExtra(FILE_NAME, fileName);
        intent.putExtra(FILE_NAME_NO_EXT, fileNameNoExt);

        startActivity(intent);
    }

    public void launchSTFTactivity(View view) {

        if (isPlayingBack.get()) {
            audioTrack.pause();
            audioTrack.flush();
            audioTrack.release();
            audioTrack = null;
            playThread.interrupt();
            playThread = null;
        }

        if (isRecordingCompleted.get()) {

            Log.i("DEBUG", "Launching next activity");

            Intent intent = new Intent(this, STFTActivity_DevMode.class);

            intent.putExtra(FILE_NAME, fileName);
            intent.putExtra(FILE_NAME_NO_EXT, fileNameNoExt);

            startActivity(intent);
        } else {
            Toast toast = Toast.makeText(getApplicationContext(), "No recording has been made", Toast.LENGTH_LONG);
            toast.show();
        }
    }

    private void requestRecordAudioPermission() {
        //check API version, do nothing if API version < 23!
        int currentapiVersion = android.os.Build.VERSION.SDK_INT;
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP){

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
        if (currentapiVersion > android.os.Build.VERSION_CODES.LOLLIPOP){

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

    @Override
    public void onRequestPermissionsResult(int requestCode, String permissions[], int[] grantResults) {
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
                return;
            }

            // other 'case' lines to check for other
            // permissions this app might request
        }
    }
}


