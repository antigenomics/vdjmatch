package com.antigenomics.vdjdb.impl.weights;

import com.antigenomics.vdjdb.impl.ClonotypeDatabase;

public interface WeightFunctionFactory {
    WeightFunction create(ClonotypeDatabase clonotypeDatabase);
}
