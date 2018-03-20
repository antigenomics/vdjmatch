package com.antigenomics.vdjdb.impl.weights;

public interface WeightFunction {
    float computeWeight(String v, String j, String cdr3aa);

    int size();
}
