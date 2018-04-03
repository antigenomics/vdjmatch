package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.WeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.sample.Sample

class ClonotypeEmbeddingCalculator {
    final ClonotypeDistanceCalculator clonotypeDistanceCalculator

    ClonotypeEmbeddingCalculator(Sample sample,
                                 SearchScope searchScope = SearchScope.EXACT,
                                 ScoringBundle scoringBundle = ScoringBundle.DUMMY,
                                 WeightFunctionFactory weightFunctionFactory = DummyWeightFunctionFactory.INSTANCE) {
        this.clonotypeDistanceCalculator = new ClonotypeDistanceCalculator(sample,
                searchScope, scoringBundle, weightFunctionFactory)
    }

    List<ClonotypeEmbedding> computeIsoMap() {

    }
}
