package com.antigenomics.vdjdb.impl.weights;

public class DummyWeightFunction implements WeightFunction {
    public static final DummyWeightFunction INSTANCE = new DummyWeightFunction();

    private DummyWeightFunction() {
    }

    @Override
    public float computeWeight(String v, String j, String cdr3aa) {
        return 1.0f;
    }
}
