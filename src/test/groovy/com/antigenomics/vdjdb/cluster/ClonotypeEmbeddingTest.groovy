package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.sequence.SearchScope
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE

class ClonotypeEmbeddingTest {
    @Test
    void simpleConnectedComponentsTest() {
        def embeddingCalc = new ClonotypeEmbeddingCalculator(TEST_SAMPLE,
                new SearchScope(1, 0, 1))
        def distanceCalc = new ClonotypeDistanceCalculator(TEST_SAMPLE,
                new SearchScope(1, 0, 1))
        def graph = distanceCalc.computeClonotypeGraph()

        int minCompSize = 5
        def cc1 = graph.connectedComponents
                .findAll { it.index.length >= minCompSize }
                .collect { it.index.length }
                .sum()

        assert embeddingCalc.isoMap(minCompSize).findAll { it.coordinates[1] != 0 }.size() == cc1
    }
}
