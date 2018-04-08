package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.filter.DummyResultFilter
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE
import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE2

class ClonotypeGraphTest {
    @Test
    void sampleTest() {
        def distanceCalc = new ClonotypeDistanceCalculator(new SearchScope(1, 0, 1),
                ScoringBundle.DUMMY, DummyWeightFunctionFactory.INSTANCE, DummyResultFilter.INSTANCE,
                false, false)
        def graph = new ClonotypeGraph(TEST_SAMPLE, distanceCalc.computeDistances(TEST_SAMPLE, TEST_SAMPLE2))

        assert graph.graph.numVertices == TEST_SAMPLE.diversity
    }

    @Test
    void simpleConnectedComponentsTest() {
        def distanceCalc = new ClonotypeDistanceCalculator(new SearchScope(1, 0, 1),
                ScoringBundle.DUMMY, DummyWeightFunctionFactory.INSTANCE, DummyResultFilter.INSTANCE,
                false, false)
        def graph = new ClonotypeGraph(TEST_SAMPLE, distanceCalc.computeDistances(TEST_SAMPLE, TEST_SAMPLE2))

        assert graph.connectedComponents.collect { it.index.length }.sum() == TEST_SAMPLE.diversity
    }
}
