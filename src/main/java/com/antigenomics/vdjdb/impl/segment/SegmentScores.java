package com.antigenomics.vdjdb.impl.segment;

public class SegmentScores {
    private final float vScore, cdr1Score, cdr2Score, jScore;

    public SegmentScores(float vScore, float cdr1Score, float cdr2Score, float jScore) {
        this.vScore = vScore;
        this.cdr1Score = cdr1Score;
        this.cdr2Score = cdr2Score;
        this.jScore = jScore;
    }

    public float getvScore() {
        return vScore;
    }

    public float getCdr1Score() {
        return cdr1Score;
    }

    public float getCdr2Score() {
        return cdr2Score;
    }

    public float getjScore() {
        return jScore;
    }
}
