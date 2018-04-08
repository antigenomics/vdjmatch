package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.ScoringProvider
import com.antigenomics.vdjdb.impl.filter.DummyResultFilter
import com.antigenomics.vdjdb.impl.filter.MaxScoreResultFilter
import com.antigenomics.vdjdb.impl.weights.DegreeWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE
import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE2
import static com.antigenomics.vdjdb.TestUtil.TEST_VDJDB_SAMPLE

class ClonotypeDistanceTest {
    @Test
    void sampleTest() {
        def distanceCalc = new ClonotypeDistanceCalculator(new SearchScope(1, 0, 1),
                ScoringBundle.DUMMY, DummyWeightFunctionFactory.INSTANCE, DummyResultFilter.INSTANCE,
                false, false)

        int n = TEST_SAMPLE.findAll { it.coding }.size()
        def distances = distanceCalc.computeDistances(TEST_SAMPLE, TEST_SAMPLE2)
        assert distances.findAll { it.score == 0 && it.from == it.to }.size() == n
        int n0 = distances.findAll { it.score == 0 }.size()
        assert n0 > n
        assert distances.size() > n0
    }

    @Test
    void sampleTest2() {
        def distanceCalc = new ClonotypeDistanceCalculator(
                new SearchScope(2, 0, 2),
                ScoringProvider.loadScoringBundle("HomoSapiens", "TRB",
                        false),
                DegreeWeightFunctionFactory.DEFAULT,
                new MaxScoreResultFilter(0.1f),
                false, false)

        def distances = distanceCalc.computeDistances(TEST_SAMPLE)

        assert !distances.isEmpty()
    }

    @Test
    void sampleTest3() {
        def distanceCalc = new ClonotypeDistanceCalculator(
                new SearchScope(2, 0, 2),
                ScoringProvider.loadScoringBundle("HomoSapiens", "TRB",
                        false),
                DegreeWeightFunctionFactory.DEFAULT,
                new MaxScoreResultFilter(0.1f),
                false, false)

        def distances = distanceCalc.computeDistances(TEST_VDJDB_SAMPLE)

        assert !distances.isEmpty()
    }
}
