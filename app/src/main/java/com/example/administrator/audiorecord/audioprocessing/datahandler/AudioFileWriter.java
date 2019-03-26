package com.example.administrator.audiorecord.audioprocessing.datahandler;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import static org.apache.commons.math3.util.FastMath.round;

public class AudioFileWriter {

    /*
    References
    ----------
    http://soundfile.sapp.org/doc/WaveFormat/
    https://stackoverflow.com/a/37436599/10699730
     */

    private int nChannels, sampleRate, byteRate, blockAlign, bitsPerSample;

    public AudioFileWriter(int nChannels, int sampleRate, int bitsPerSample) {
        this.nChannels = nChannels;
        this.sampleRate = sampleRate;
        this.bitsPerSample = bitsPerSample;

        byteRate = sampleRate * nChannels * bitsPerSample / 8;

        blockAlign = nChannels * bitsPerSample / 8;
    }

    public void convertPcmToWav(final File rawFile, final File waveFile) throws IOException {

        byte[] rawData = new byte[(int) rawFile.length()];

        try (DataInputStream input = new DataInputStream(new FileInputStream(rawFile))) {
            input.read(rawData);
        }

        try (DataOutputStream output = new DataOutputStream(new FileOutputStream(waveFile))) {
            // WAVE header
            writeString(output, "RIFF"); // chunk id
            writeInt(output, 36 + rawData.length); // chunk size
            writeString(output, "WAVE"); // format
            writeString(output, "fmt "); // subchunk 1 id
            writeInt(output, 16); // subchunk 1 size
            writeShort(output, (short) 1); // audio format (1 = PCM)
            writeShort(output, (short) nChannels); // number of channels
            writeInt(output, sampleRate); // sample rate
            writeInt(output, byteRate); // byte rate
            writeShort(output, (short) blockAlign); // block align
            writeShort(output, (short) bitsPerSample); // bits per sample
            writeString(output, "data"); // subchunk 2 id
            writeInt(output, rawData.length); // subchunk 2 size

            // Audio data (conversion big endian -> little endian)
            short[] shorts = new short[rawData.length / 2];
            ByteBuffer.wrap(rawData).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer().get(shorts);
            ByteBuffer bytes = ByteBuffer.allocate(shorts.length * 2);
            for (short s : shorts) {
                bytes.putShort(s);
            }

            output.write(fullyReadFileToBytes(rawFile));
        }
    }

    public void convertShortArrayToFile(final short[] shorts, final File waveFile) throws IOException {

        try (DataOutputStream output = new DataOutputStream(new FileOutputStream(waveFile))) {
            // WAVE header
            writeString(output, "RIFF"); // chunk id
            writeInt(output, 36 + shorts.length * 2); // chunk size
            writeString(output, "WAVE"); // format
            writeString(output, "fmt "); // subchunk 1 id
            writeInt(output, 16); // subchunk 1 size
            writeShort(output, (short) 1); // audio format (1 = PCM)
            writeShort(output, (short) nChannels); // number of channels
            writeInt(output, sampleRate); // sample rate
            writeInt(output, byteRate); // byte rate
            writeShort(output, (short) blockAlign); // block align
            writeShort(output, (short) bitsPerSample); // bits per sample
            writeString(output, "data"); // subchunk 2 id
            writeInt(output, shorts.length * 2); // subchunk 2 size

            for (short aShort : shorts) {
                writeShort(output, aShort);
            }

        }
    }

    private byte[] fullyReadFileToBytes(File f) throws IOException {
        int size = (int) f.length();
        byte bytes[] = new byte[size];
        byte tmpBuff[] = new byte[size];
        try (FileInputStream fis = new FileInputStream(f)) {

            int read = fis.read(bytes, 0, size);
            if (read < size) {
                int remain = size - read;
                while (remain > 0) {
                    read = fis.read(tmpBuff, 0, remain);
                    System.arraycopy(tmpBuff, 0, bytes, size - remain, read);
                    remain -= read;
                }
            }
        }

        return bytes;
    }

    public void doubleArrayToPCM(final double[] doubles, final File waveFile, int sampleLength) throws IOException {

        short[] shorts = new short[sampleLength];

        for (int i = 0; i < sampleLength; i++) {
            shorts[i] = (short) round(doubles[i] * 32768.0);
        }

        try (DataOutputStream output = new DataOutputStream(new FileOutputStream(waveFile))) {

            for (short aShort : shorts) {
                writeShort(output, aShort);
            }

        }
    }

    public void doubleMultichannelArrayToPCM(final double[][] doubles, final File pcmFile, int sampleLength) throws IOException {

        short[] shorts = new short[sampleLength * nChannels];


        for (int i = 0; i < sampleLength; i++) {
            for (int c = 0; c < nChannels; c++) {
                shorts[i * nChannels + c] = (short) round(doubles[c][i] * 32768.0);
            }
        }


        try (DataOutputStream output = new DataOutputStream(new FileOutputStream(pcmFile))) {

            for (short aShort : shorts) {
                writeShort(output, aShort);
            }

        }
    }

    private void writeInt(final DataOutputStream output, final int value) throws IOException {
        output.write(value);
        output.write(value >> 8);
        output.write(value >> 16);
        output.write(value >> 24);
    }

    private void writeShort(final DataOutputStream output, final short value) throws
            IOException {
        output.write(value);
        output.write(value >> 8);
    }

    private void writeString(final DataOutputStream output, final String value) throws
            IOException {
        for (int i = 0; i < value.length(); i++) {
            output.write(value.charAt(i));
        }
    }
}

