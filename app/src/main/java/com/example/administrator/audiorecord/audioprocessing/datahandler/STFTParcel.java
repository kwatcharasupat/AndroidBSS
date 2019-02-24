package com.example.administrator.audiorecord.audioprocessing.datahandler;

import android.os.Parcel;
import android.os.Parcelable;

import org.apache.commons.math3.complex.Complex;

public class STFTParcel implements Parcelable {

    private int nChannels, nFrames, nFreqs, audioDataLength, winLen, nOverlap;
    private String winFunc;

    public STFTParcel(int nChannels, int nFrames, int nFreqs, int audioDataLength, int winLen, int nOverlap, String winFunc) {
        this.nChannels = nChannels;
        this.nFrames = nFrames;
        this.nFreqs = nFreqs;
        this.audioDataLength = audioDataLength;
        this.winLen = winLen;
        this.nOverlap = nOverlap;
        this.winFunc = winFunc;
    }

    public STFTParcel(Parcel in) {
        nChannels = in.readInt();
        nFrames = in.readInt();
        nFreqs = in.readInt();
        audioDataLength = in.readInt();
        winLen = in.readInt();
        nOverlap = in.readInt();
        winFunc = in.readString();
    }

    public static final Parcelable.Creator<STFTParcel> CREATOR
            = new Parcelable.Creator<STFTParcel>() {
        public STFTParcel createFromParcel(Parcel in) {
            return new STFTParcel(in);
        }

        public STFTParcel[] newArray(int size) {
            return new STFTParcel[size];
        }
    };

    @Override
    public int describeContents() {
        return 0;
    }

    @Override
    public void writeToParcel(Parcel dest, int flags) {
        dest.writeInt(nChannels);
        dest.writeInt(nFrames);
        dest.writeInt(nFreqs);
        dest.writeInt(audioDataLength);
        dest.writeInt(winLen);
        dest.writeInt(nOverlap);
        dest.writeString(winFunc);
    }

    public int getnChannels(){
        return nChannels;
    }

    public int getAudioDataLength() {
        return audioDataLength;
    }

    public int getnFrames() {
        return nFrames;
    }

    public int getnFreqs() {
        return nFreqs;
    }

    public int getnOverlap() {
        return nOverlap;
    }

    public int getWinLen() {
        return winLen;
    }

    public String getWinFunc() {
        return winFunc;
    }
}

