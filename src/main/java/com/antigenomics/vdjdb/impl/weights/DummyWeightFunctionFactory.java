package com.antigenomics.vdjdb.impl.weights;

import com.antigenomics.vdjdb.impl.ClonotypeDatabase;

public class DummyWeightFunctionFactory implements WeightFunctionFactory {
    public static final DummyWeightFunctionFactory INSTANCE = new DummyWeightFunctionFactory();

    private DummyWeightFunctionFactory() {

    }

    @Override
    public WeightFunction create(ClonotypeDatabase clonotypeDatabase) {
        return DummyWeightFunction.INSTANCE;
    }
}
