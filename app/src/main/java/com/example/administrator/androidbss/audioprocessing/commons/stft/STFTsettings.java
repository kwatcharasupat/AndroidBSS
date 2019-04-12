package com.example.administrator.androidbss.audioprocessing.commons.stft;

import static com.example.administrator.androidbss.audioprocessing.commons.stft.WindowFunctions.PERIODIC_HAMMING_WINDOW;
import static org.apache.commons.math3.util.FastMath.round;

public class STFTsettings {

    private static final int DEFAULT_WINDOW_LENGTH = 2048;

    int winLen = -1;
    int nOverlap = -1;
    int winFunc = -1;

    public STFTsettings() {

    }

    public STFTsettings setWindowLength(int winLen) {
        this.winLen = winLen;
        return this;
    }

    public STFTsettings setWindowFunction(int winFunc) {
        this.winFunc = winFunc;
        return this;
    }

    public STFTsettings setOverlap(int nOverlap) {
        this.nOverlap = nOverlap;
        return this;
    }

    public STFTsettings setOverlapToDefault(){
        return this;
    }

    public STFTsettings prepare() {
        if (winLen < 0) {
            winLen = DEFAULT_WINDOW_LENGTH;
        }

        if (winFunc < 0) {
            winFunc = PERIODIC_HAMMING_WINDOW;
        }

        if (nOverlap < 0) {
            nOverlap = (int) round(WindowFunctions.getDefaultOverlap(winFunc) * winLen);
        }

        return this;
    }

}
