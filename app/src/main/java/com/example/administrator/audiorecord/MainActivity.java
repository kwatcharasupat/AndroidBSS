
package com.example.administrator.audiorecord;

import android.content.Intent;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioRecord;
import android.media.AudioTrack;
import android.media.MediaRecorder;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;

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

public class MainActivity extends AppCompatActivity {

    private static final int SAMPLING_RATE_IN_HZ = 16000;   //44100 Hz is supported on all devices

    private static final int CHANNEL_CONFIG = AudioFormat.CHANNEL_IN_STEREO;

    private static final int AUDIO_FORMAT = AudioFormat.ENCODING_PCM_16BIT;

    private static final int BUFFER_SIZE_FACTOR = 2; // preemptively allocated space

    private static final int BUFFER_SIZE = AudioRecord.getMinBufferSize(
            SAMPLING_RATE_IN_HZ,
            CHANNEL_CONFIG,
            AUDIO_FORMAT) * BUFFER_SIZE_FACTOR;



    private final AtomicBoolean recordingInProgress = new AtomicBoolean(false);
    private final AtomicBoolean playingInProgress = new AtomicBoolean(false);

    private AudioRecord recorder = null;
    private Thread recordingThread = null;

    private Button startButton;
    private Button stopButton;
    private Button playButton;

    public static final String FILE_NAME = "com.example.administrator.FILE_NAME";
    private String fileName;

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        startButton = findViewById(R.id.bStartRecord);
        startButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                startRecording();
                startButton.setEnabled(false);
                stopButton.setEnabled(true);
                playButton.setEnabled(false);
            }
        });

        stopButton = findViewById(R.id.bStopRecord);
        stopButton.setEnabled(false);
        stopButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                stopRecording();
                startButton.setEnabled(true);
                stopButton.setEnabled(false);
                playButton.setEnabled(true);
            }
        });

        playButton = findViewById(R.id.bPlayRecord);
        playButton.setEnabled(false);
        playButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

                String filePath = getExternalStorageDirectory().getAbsolutePath() + "/" + fileName;
                playRecording(filePath);
                playButton.setEnabled(false);
                stopButton.setEnabled(false);
                startButton.setEnabled(true);
            }
        });
    }

    @Override
    protected void onResume() {
        super.onResume();

        startButton.setEnabled(true);
        stopButton.setEnabled(false);
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
        recordingInProgress.set(true);
        recordingThread = new Thread(new RecordingRunnable(), "Recording Thread");
        recordingThread.start();
    }

    private void stopRecording() {
        if (null == recorder) {
            return;
        }

        recordingInProgress.set(false);
        recorder.stop();
        recorder.release();
        recorder = null;
        recordingThread = null;

        Log.i("DEBUG","Recording stopped");
    }


    private class RecordingRunnable implements Runnable {

        @Override
        public void run() {

            Log.i("DEBUG","Now recording");

            SimpleDateFormat formatter = new SimpleDateFormat("yyMMddHHmmss");
            Date dateTimeNow = new Date();
            fileName = "recording" + formatter.format(dateTimeNow) + ".pcm";

            final File file = new File(getExternalStorageDirectory().getAbsolutePath(), fileName);
            final ByteBuffer buffer = ByteBuffer.allocateDirect(BUFFER_SIZE);

            try (final FileOutputStream outStream = new FileOutputStream(file)) {
                while (recordingInProgress.get()) {
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

    private void playRecording(String filePath) {

        Log.i("DEBUG","Now playing");

        playingInProgress.set(true);

        if (filePath == null)
            return;

        int intSize = android.media.AudioTrack.getMinBufferSize(
                SAMPLING_RATE_IN_HZ,
                AudioFormat.CHANNEL_IN_STEREO,
                AudioFormat.ENCODING_PCM_16BIT);

        AudioTrack audioTrack = new AudioTrack(
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
        audioTrack.stop();
        audioTrack.release();

        playingInProgress.set(false);

        Log.i("DEBUG","Playback ended");
    }

    public void launchSTFTcomplexActivity(View view) {
        Intent intent = new Intent(this, STFTActivity.class);

        intent.putExtra(FILE_NAME, fileName);

        startActivity(intent);
    }

}


