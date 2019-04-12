package com.example.administrator.androidbss;

import android.content.Intent;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;

import com.example.administrator.androidbss.activities.MainActivity;

public class ModeActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_mode);
    }

    public void launchUserMode(View view){
        Intent intent = new Intent(this, MainActivity.class);
        startActivity(intent);
    }

    public void launchDevMode(View view){

    }
}
