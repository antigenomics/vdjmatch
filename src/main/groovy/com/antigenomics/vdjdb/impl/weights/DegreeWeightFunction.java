package com.antigenomics.vdjdb.impl.weights;

import java.util.Map;

public class DegreeWeightFunction implements WeightFunction {
    private final Map<String, Float> weights;

    public DegreeWeightFunction(Map<String, Float> weights) {
        this.weights = weights;
    }

    @Override
    public float computeWeight(String v, String j, String cdr3aa) {
        return weights.get(cdr3aa);
    }

    @Override
    public int size() {
        return weights.size();
    }
}
