package com.antigenomics.vdjdb.impl.model;


public class CloglogAggregateScoring extends LinearAggregateScoring {
    private final float precomputedZeroFullScore;

    public CloglogAggregateScoring(float intercept, float cc1, float cc2, float cc3, float cv, float cj) {
        super(intercept, cc1, cc2, cc3, cv, cj);
        this.precomputedZeroFullScore = this.computeFullScore(0, 0, 0, 0, 0);
    }

    @Override
    public float getPrecomputedZeroFullScore() {
        return this.precomputedZeroFullScore;
    }

    @Override
    public float computeFullScore(float vScore, float cdr1Score, float cdr2Score, float cdr3Score, float jScore) {
        //System.out.println(super.computeFullScore(vScore, cdr1Score, cdr2Score, cdr3Score, jScore));
        return 1.0f - (float) Math.exp(-Math.exp(super.computeFullScore(vScore, cdr1Score, cdr2Score, cdr3Score, jScore)));
    }
}
