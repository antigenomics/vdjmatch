package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.WeightFunctionFactory
import com.antigenomics.vdjdb.sequence.ScoringType
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.sample.Sample

/**
 * Compute pairwise distance (dissimilarity) between clonotypes in the same sample / different pair of sample
 */
class ClonotypeDistanceCalculator {
    final ClonotypeDatabase clonotypeDatabase
    final Sample sample
    final boolean probabilisticScoring

    ClonotypeDistanceCalculator(Sample sample,
                                SearchScope searchScope = SearchScope.EXACT,
                                ScoringBundle scoringBundle = ScoringBundle.DUMMY,
                                WeightFunctionFactory weightFunctionFactory = DummyWeightFunctionFactory.INSTANCE) {
        this.sample = sample
        this.clonotypeDatabase = VdjdbInstance.fromSample(sample,
                searchScope, scoringBundle, weightFunctionFactory)
        this.probabilisticScoring = scoringBundle.alignmentScoring.scoringType == ScoringType.Probabilistic
    }

    List<ClonotypeDistance> computeDistancesTo(Sample otherSample) {
        def results = new ArrayList<ClonotypeDistance>()
        clonotypeDatabase.search(otherSample).each { entry ->
            entry.value.each { result ->
                int idInSample = result.row[VdjdbInstance.CLONOTYPE_SAMPLE_ID_COL].value.toInteger()
                results.add(new ClonotypeDistance(
                        idInSample,
                        result.id, // id in other sample
                        sample[idInSample],
                        entry.key, // clonotype in other sample
                        result.score,
                        result.weight,
                        probabilisticScoring ? (1.0d - result.score) : (-result.score)
                ))
            }
        }
        results
    }

    /**
     * Computes distances between clonotypes in the sample and converts them to clonotype graph.
     * @return clonotype graph object
     */
    ClonotypeGraph computeClonotypeGraph() {
        new ClonotypeGraph(sample, computeDistancesTo(sample))
    }
}
