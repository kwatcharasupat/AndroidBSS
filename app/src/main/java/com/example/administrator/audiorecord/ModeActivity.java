package com.example.administrator.audiorecord;

import android.content.Intent;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;

import com.example.administrator.audiorecord.activities_devmode.MainActivity_DevMode;
import com.example.administrator.audiorecord.activities_usermode.MainActivity_UserMode;

public class ModeActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_mode);
    }

    public void launchUserMode(View view){
        Intent intent = new Intent(this, MainActivity_UserMode.class);
        startActivity(intent);
    }

    public void launchDevMode(View view){
        Intent intent = new Intent(this, MainActivity_DevMode.class);
        startActivity(intent);
    }
}
