package com.antigenomics.vdjdb.impl.model;


public class CloglogAggregateScoring implements AggregateScoring {
    private final float Cc1, Cc2, Cv, Cj;

    public CloglogAggregateScoring(float cc1, float cc2, float cv, float cj) {
        Cc1 = cc1;
        Cc2 = cc2;
        Cv = cv;
        Cj = cj;
    }

    public float getCc1() {
        return Cc1;
    }

    public float getCc2() {
        return Cc2;
    }

    public float getCv() {
        return Cv;
    }

    public float getCj() {
        return Cj;
    }

    @Override
    public float computeFullScore(float vScore, float cdr1Score, float cdr2Score, float cdr3Score, float jScore) {
        float x = Cv * vScore + Cc1 * cdr1Score + Cc2 * cdr2Score + cdr3Score + Cj * jScore;
        return 1.0f - (float) Math.exp(-Math.exp(x));
    }
}
