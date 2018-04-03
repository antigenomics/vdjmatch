package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.WeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.sample.Sample

class ClonotypeDistanceCalculator {
    final ClonotypeDatabase clonotypeDatabase
    final Sample target

    ClonotypeDistanceCalculator(Sample target,
                                SearchScope searchScope = SearchScope.EXACT,
                                ScoringBundle scoringBundle = ScoringBundle.DUMMY,
                                WeightFunctionFactory weightFunctionFactory = DummyWeightFunctionFactory.INSTANCE) {
        this.target = target
        this.clonotypeDatabase = VdjdbInstance.fromSample(target,
                searchScope, scoringBundle, weightFunctionFactory)
    }

    List<ClonotypeDistance> computeDistances(Sample sample) {
        def results = new ArrayList<ClonotypeDistance>()
        clonotypeDatabase.search(sample).each { entry ->
            entry.value.each { result ->
                int idInTarget = result.row[VdjdbInstance.CLONOTYPE_SAMPLE_ID_COL].value.toInteger()
                results.add(new ClonotypeDistance(
                        result.id,
                        idInTarget,
                        entry.key,
                        target[idInTarget],
                        result.score,
                        result.weight,
                        result.weightedScore
                ))
            }
        }
        results
    }

    List<ClonotypeDistance> computeDistances() {
        computeDistances(target)
    }
}
