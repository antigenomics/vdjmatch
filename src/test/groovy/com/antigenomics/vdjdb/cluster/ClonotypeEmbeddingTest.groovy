package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.sequence.SearchScope
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.TEST_SAMPLE

class ClonotypeEmbeddingTest {
    @Test
    void simpleConnectedComponentsTest() {
        def distanceCalc = new ClonotypeDistanceCalculator(new SearchScope(1, 0, 1))
        def graph = new ClonotypeGraph(TEST_SAMPLE, distanceCalc.computeDistances(TEST_SAMPLE))

        int minCompSize = 5
        def cc1 = graph.connectedComponents
                .findAll { it.index.length >= minCompSize }
                .collect { it.index.length }
                .sum()

        //assert ClonotypeEmbeddingCalculator.isoMap(graph,
        //        minCompSize).findAll { it.coordinates[1] != 0d || it.coordinates[2] != 0d }.size() == cc1
    }
}
