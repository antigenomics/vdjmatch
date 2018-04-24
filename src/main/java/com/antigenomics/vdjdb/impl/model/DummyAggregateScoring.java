package com.antigenomics.vdjdb.impl.model;

public class DummyAggregateScoring implements AggregateScoring {
    public static final DummyAggregateScoring INSTANCE = new DummyAggregateScoring();

    private DummyAggregateScoring() {

    }

    @Override
    public float getPrecomputedZeroFullScore() {
        return 1.0f;
    }

    @Override
    public float computeFullScore(float vScore, float cdr1Score, float cdr2Score, float cdr3Score, float jScore) {
        return vScore + cdr1Score + cdr2Score + cdr3Score + jScore;
    }
}
