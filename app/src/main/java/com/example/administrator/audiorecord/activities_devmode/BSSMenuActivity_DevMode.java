package com.example.administrator.audiorecord.activities_devmode;

import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.AdapterView.OnItemSelectedListener;
import android.widget.Spinner;
import android.widget.Toast;

import com.example.administrator.audiorecord.R;
import com.example.administrator.audiorecord.audioprocessing.datahandler.STFTParcel;

import org.apache.commons.math3.complex.Complex;

import static com.example.administrator.audiorecord.activities_devmode.STFTActivity_DevMode.STFT_PARCEL;

public class BSSMenuActivity_DevMode extends AppCompatActivity implements OnItemSelectedListener{

    public static Complex[][][] STFTdata;
    int audioDataLength, nChannels, nFrames, nFreqs, winLen, nOverlap;
    String winFunc;

    String BSStype;

    Spinner spinnerBSS;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.devmodeactivity_bssmenu);

        getDataFromIntent();

        setSpinnerBSS();
    }



    void getDataFromIntent() {
        Intent intent = getIntent();

        STFTParcel stftParcel = intent.getParcelableExtra(STFT_PARCEL);

        STFTdata = STFTActivity_DevMode.STFToutput;

        nChannels = stftParcel.getnChannels();

        nFrames = stftParcel.getnFrames();

        nFreqs = stftParcel.getnFreqs();

        winLen = stftParcel.getWinLen();

        nOverlap = stftParcel.getnOverlap();

        audioDataLength = stftParcel.getAudioDataLength();

        winFunc = stftParcel.getWinFunc();
    }


    void setSpinnerBSS() {
        spinnerBSS = findViewById(R.id.spinnerBSS);
        spinnerBSS.setOnItemSelectedListener(this);

        ArrayAdapter<CharSequence> adapter = ArrayAdapter.createFromResource(this,
                R.array.bss_options, android.R.layout.simple_spinner_item);
        adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerBSS.setAdapter(adapter);

    }

    @Override
    public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
        BSStype = parent.getItemAtPosition(position).toString();
        Toast.makeText(parent.getContext(), "Selected: " + BSStype, Toast.LENGTH_LONG).show();
    }

    @Override
    public void onNothingSelected(AdapterView<?> arg0) {    }

    public void launchBSSActivity(View view) {

        Log.i("DEBUG", "Launching next activity");

        STFTParcel stftParcel = new STFTParcel(nChannels, nFrames, nFreqs, audioDataLength, winLen, nOverlap, winFunc);

        Intent intent;

        switch (BSStype){
            case "AuxIVA - Apache Commons":
                intent = new Intent(this, AuxIVAActivity_DevMode.class);
                break;
            case "AuxIVA - Parallel":
                intent = new Intent(this, AuxIVAParallelActivity_DevMode.class);
                break;
            default:
                return;
        }

        intent.putExtra(STFT_PARCEL, stftParcel);

        startActivity(intent);
    }

}
