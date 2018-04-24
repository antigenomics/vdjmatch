package com.antigenomics.vdjdb.impl.model;


public class LinearAggregateScoring implements AggregateScoring {
    private final float intercept, cc1, cc2, cc3, cv, cj;
    private final float precomputedZeroFullScore;

    public LinearAggregateScoring(float intercept,
                                  float cc1, float cc2, float cc3,
                                  float cv, float cj) {
        this.intercept = intercept;
        this.cc1 = cc1;
        this.cc2 = cc2;
        this.cc3 = cc3;
        this.cv = cv;
        this.cj = cj;
        this.precomputedZeroFullScore = this.computeFullScore(0, 0, 0, 0, 0);
    }

    public float getIntercept() {
        return intercept;
    }

    public float getCc1() {
        return cc1;
    }

    public float getCc2() {
        return cc2;
    }

    public float getCc3() {
        return cc3;
    }

    public float getCv() {
        return cv;
    }

    public float getCj() {
        return cj;
    }

    @Override
    public float getPrecomputedZeroFullScore() {
        return precomputedZeroFullScore;
    }

    @Override
    public float computeFullScore(float vScore, float cdr1Score, float cdr2Score, float cdr3Score, float jScore) {
        return intercept + cv * vScore + cc1 * cdr1Score + cc2 * cdr2Score + cc3 * cdr3Score + cj * jScore;
    }
}
