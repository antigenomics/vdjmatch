package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.io.InputStreamFactory
import com.antigenomics.vdjtools.io.SampleStreamConnection
import org.junit.Test

import java.util.zip.GZIPInputStream

import static com.antigenomics.vdjdb.Util.resourceAsStream

class ClonotypeDistanceTest {
    @Test
    void clonotypesFromSampleTest() {
        def sample = SampleStreamConnection.load([
                create: {
                    new GZIPInputStream(resourceAsStream("sergey_anatolyevich.gz"))
                },
                getId : { "sergey_anatolyevich.gz" }
        ] as InputStreamFactory)

        def distanceCalc = new ClonotypeDistanceCalculator(sample,
                new SearchScope(1, 0, 1))

        int n = sample.findAll { it.coding }.size()
        def distances = distanceCalc.computeDistances()
        assert distances.findAll { it.score == 0 && it.idTarget == it.idQuery }.size() == n
        int n0 = distances.findAll { it.score == 0 }.size()
        assert n0 > n
        assert distances.size() > n0
    }
}
