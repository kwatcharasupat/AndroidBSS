package com.example.administrator.audiorecord.activities;

import android.annotation.SuppressLint;
import android.content.Context;
import android.content.Intent;
import android.net.ConnectivityManager;
import android.net.NetworkInfo;
import android.os.AsyncTask;
import android.os.Bundle;
import android.support.design.widget.FloatingActionButton;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Spinner;
import android.widget.TextView;

import com.example.administrator.audiorecord.R;
import com.google.api.gax.core.FixedCredentialsProvider;
import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.speech.v1.RecognitionAudio;
import com.google.cloud.speech.v1.RecognitionConfig;
import com.google.cloud.speech.v1.RecognizeResponse;
import com.google.cloud.speech.v1.SpeechClient;
import com.google.cloud.speech.v1.SpeechRecognitionAlternative;
import com.google.cloud.speech.v1.SpeechRecognitionResult;
import com.google.cloud.speech.v1.SpeechSettings;
import com.google.protobuf.ByteString;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import static android.os.Environment.getExternalStorageDirectory;
import static com.example.administrator.audiorecord.activities.BSSActivity.FILE_NAME_NO_EXT;
import static com.example.administrator.audiorecord.activities.BSSActivity.NUM_CHANNELS;
import static com.example.administrator.audiorecord.activities.BSSActivity.SAMPLING_RATE;

public class ASRActivity extends AppCompatActivity implements AdapterView.OnItemSelectedListener {

    final int ALL_SOURCES_POSITION = 0;

    String[] transcriptions;
    StringBuilder stringBuilder;

    TextView tvTranscriptions;

    Spinner spinnerSrcTranscribe;
    boolean readAllSources;
    int currentSource;

    AtomicBoolean[] isSourceTranscibed;
    AtomicBoolean hasTranscriptionStarted = new AtomicBoolean(false);
    AtomicBoolean isAllSourcesTranscribed = new AtomicBoolean(false);

    int NUM_CHANNELS_DEFAULT = 2, SAMPLING_RATE_DEFAULT = 16000, BIT_DEPTH = 16;

    String fileNameNoExt;

    int nSrc, samplingRate;

    File[] wavFiles;
    File pcmMultichannel, wavMultichannel;

    FloatingActionButton fabASR;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.usermodeactivity_asr);

        Intent intent = getIntent();


        tvTranscriptions = findViewById(R.id.textviewTranscription);

        fileNameNoExt = intent.getStringExtra(FILE_NAME_NO_EXT);
        nSrc = intent.getIntExtra(NUM_CHANNELS, NUM_CHANNELS_DEFAULT);
        samplingRate = intent.getIntExtra(SAMPLING_RATE, SAMPLING_RATE_DEFAULT);

        setSourceSpinner();

        transcriptions = new String[nSrc];
        isSourceTranscibed = new AtomicBoolean[nSrc];

        for (int s = 0; s < nSrc; s++) {
            isSourceTranscibed[s] = new AtomicBoolean(false);
        }

        wavFiles = new File[nSrc];

        for (int s = 0; s < nSrc; s++) {
            wavFiles[s] = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_Source" + (s + 1) + ".wav");
        }

        wavMultichannel = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_AllSources.wav");

        pcmMultichannel = new File(getExternalStorageDirectory().getAbsolutePath(), fileNameNoExt + "_AllSources.pcm");

        tvTranscriptions.setText(getString(R.string.press_to_transcribe));

        fabASR = findViewById(R.id.fabASR);
        fabASR.setOnClickListener(v -> {
            if (isNetworkAvailable()) {
                tvTranscriptions.setText(getString(R.string.transcription_in_progress));
                spinnerSrcTranscribe.setSelection(ALL_SOURCES_POSITION);
                hasTranscriptionStarted.set(true);
                new SpeechTask().execute();
            } else {
                tvTranscriptions.setText(getString(R.string.no_internet));
            }
        });
    }


    private boolean isNetworkAvailable() {
        ConnectivityManager connectivityManager
                = (ConnectivityManager) getSystemService(Context.CONNECTIVITY_SERVICE);
        NetworkInfo activeNetworkInfo = connectivityManager.getActiveNetworkInfo();
        return activeNetworkInfo != null && activeNetworkInfo.isConnected();
    }

    @SuppressLint("StaticFieldLeak")
    private class SpeechTask extends AsyncTask<Void, Void, Void> {

        @Override
        protected Void doInBackground(Void... voids) {


            for (int s = 0; s < nSrc; s++) {

                stringBuilder = new StringBuilder();

                try {
                    transcribe(wavFiles[s].getAbsolutePath(), s);
                } catch (Exception e) {
                    e.printStackTrace();
                }

                isSourceTranscibed[s].set(true);
            }

            isAllSourcesTranscribed.set(true);

            runOnUiThread(ASRActivity.this::readTranscription);

            return null;
        }

    }

    public void transcribe(String fileName, int source) throws Exception {

        Log.i("DEBUG", "Transcibing");

        Path path = Paths.get(fileName);
        byte[] content = Files.readAllBytes(path);

        try (SpeechClient speechClient = SpeechClient.create(getSpeedSettingsWithCredential())) {

            Log.i("DEBUG", "Trying");

            // Get the contents of the local audio file
            RecognitionAudio recognitionAudio =
                    RecognitionAudio.newBuilder().setContent(ByteString.copyFrom(content)).build();

            // Configure request to enable multiple channels
            RecognitionConfig config =
                    RecognitionConfig.newBuilder()
                            .setEncoding(RecognitionConfig.AudioEncoding.LINEAR16)
                            .setLanguageCode("en-US")
                            .build();

            Log.i("DEBUG", "Config done");

            // Perform the transcription request
            RecognizeResponse recognizeResponse = speechClient.recognize(config, recognitionAudio);

            Log.i("DEBUG", "requested");


            // Print out the results
            for (SpeechRecognitionResult result : recognizeResponse.getResultsList()) {
                Log.i("DEBUG", "for loop");
                // There can be several alternative transcripts for a given chunk of speech. Just use the
                // first (most likely) one here.
                SpeechRecognitionAlternative alternative = result.getAlternatives(0);
                stringBuilder.append(alternative.getTranscript());
                //Log.i("Transcription:", String.valueOf(result.getChannelTag()) + ": " + alternative.getTranscript());
            }

            stringBuilder.append("\n\n");
            transcriptions[source] = stringBuilder.toString();

            Log.i("DEBUG", "done");
        }
    }


    public SpeechSettings getSpeedSettingsWithCredential() {

        SpeechSettings speechSettings =
                null;
        try {
            InputStream stream = getResources().openRawResource(R.raw.credential);
            FixedCredentialsProvider credentialsProvider = FixedCredentialsProvider.create(GoogleCredentials.fromStream(stream));
            speechSettings = SpeechSettings.newBuilder()
                    .setCredentialsProvider(credentialsProvider)
                    .build();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return speechSettings;
    }


    void setSourceSpinner() {
        spinnerSrcTranscribe = findViewById(R.id.spinnerSrcTranscribe);
        spinnerSrcTranscribe.setOnItemSelectedListener(this);

        List<CharSequence> sourceArray = new ArrayList<>(nSrc + 1);
        sourceArray.add("All");

        for (int s = 0; s < nSrc; s++) {
            Log.i("DEBUG", "entered loop");
            sourceArray.add("Source " + (s + 1));
            Log.i("DEBUG", "added source" + s + "into the spinner");
        }

        ArrayAdapter<CharSequence> adapter = new ArrayAdapter<>(this, android.R.layout.simple_spinner_item, sourceArray);

        adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerSrcTranscribe.setAdapter(adapter);
    }


    @Override
    public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
        if (hasTranscriptionStarted.get()) {
            if (position > 0) {
                readAllSources = false;
                currentSource = position - 1;
                if (isSourceTranscibed[currentSource].get()) {
                    readTranscription();
                } else {
                    tvTranscriptions.setText(getString(R.string.transcription_in_progress));
                }


            } else {
                readAllSources = true;
                if (isAllSourcesTranscribed.get()) {
                    readTranscription();
                } else if (hasTranscriptionStarted.get()) {
                    tvTranscriptions.setText(getString(R.string.transcription_in_progress));
                }
            }
        } else {
            tvTranscriptions.setText(getString(R.string.press_to_transcribe));
        }
    }

    @Override
    public void onNothingSelected(AdapterView<?> parent) {

    }

    void readTranscription() {
        if (readAllSources) {
            stringBuilder = new StringBuilder();
            for (int s = 0; s < nSrc; s++) {
                stringBuilder.append("Source ").append(s + 1).append(":\n");
                stringBuilder.append(transcriptions[s]);
                stringBuilder.append("\n\n");
            }
            tvTranscriptions.setText(stringBuilder.toString());
        } else {
            tvTranscriptions.setText(transcriptions[currentSource]);
        }
    }
}
