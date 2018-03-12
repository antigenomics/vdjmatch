package com.antigenomics.vdjdb.impl.model;

public interface AggregateScoring {
    float computeFullScore(float vScore, float cdr1Score, float cdr2Score, float cdr3Score, float jScore);
}
