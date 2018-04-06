package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.sequence.SearchScope
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE

class ClonotypeGraphTest {
    @Test
    void sampleTest() {
        def distanceCalc = new ClonotypeDistanceCalculator(new SearchScope(1, 0, 1))
        def graph = new ClonotypeGraph(TEST_SAMPLE, distanceCalc.computeDistances(TEST_SAMPLE))

        assert graph.graph.numVertices == TEST_SAMPLE.diversity
    }

    @Test
    void simpleConnectedComponentsTest() {
        def distanceCalc = new ClonotypeDistanceCalculator(new SearchScope(1, 0, 1))
        def graph = new ClonotypeGraph(TEST_SAMPLE, distanceCalc.computeDistances(TEST_SAMPLE))

        assert graph.connectedComponents.collect { it.index.length }.sum() == TEST_SAMPLE.diversity
    }
}
