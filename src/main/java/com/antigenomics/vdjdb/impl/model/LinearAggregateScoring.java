package com.antigenomics.vdjdb.impl.model;


public class LinearAggregateScoring implements AggregateScoring {
    private final float intercept, cc1, cc2, cv, cj;

    public LinearAggregateScoring(float intercept,
                                  float cc1, float cc2, float cv, float cj) {
        this.intercept = intercept;
        this.cc1 = cc1;
        this.cc2 = cc2;
        this.cv = cv;
        this.cj = cj;
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

    public float getCv() {
        return cv;
    }

    public float getCj() {
        return cj;
    }

    @Override
    public float computeFullScore(float vScore, float cdr1Score, float cdr2Score, float cdr3Score, float jScore) {
        return intercept + cv * vScore + cc1 * cdr1Score + cc2 * cdr2Score + cdr3Score + cj * jScore;
    }
}
