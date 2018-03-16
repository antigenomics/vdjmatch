package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.impl.filter.DummyResultFilter
import com.antigenomics.vdjdb.impl.filter.MaxScoreResultFilter
import com.antigenomics.vdjdb.impl.filter.ScoreThresholdResultFilter
import com.antigenomics.vdjdb.impl.filter.TopNResultFilter
import org.junit.Test

class FilterTest {
    @Test
    void dummyFilterTest() {
        def filter = new DummyResultFilter()

        def results = [
                new ClonotypeSearchResult(null, null, -1, 0.1f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.1f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.1f, 1)
        ]

        assert filter.filter(results).size() == 3
    }

    @Test
    void thresholdFilterTest() {
        def filter = new ScoreThresholdResultFilter(0.2f)

        def results = [
                new ClonotypeSearchResult(null, null, -1, 0.5f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.3f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.1f, 1)
        ]

        assert filter.filter(results).size() == 2
    }

    @Test
    void maxFilterTest() {
        def filter = new MaxScoreResultFilter(0.2f)

        def results = [
                new ClonotypeSearchResult(null, null, -1, 0.5f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.3f, 10),
                new ClonotypeSearchResult(null, null, -1, 0.1f, 1)
        ]

        assert filter.filter(results).first().weightedScore == 3f
    }

    @Test
    void topNFilterTest() {
        def filter = new TopNResultFilter(0.11f, 3)

        def results = [
                new ClonotypeSearchResult(null, null, -1, 0.5f, 10),
                new ClonotypeSearchResult(null, null, -1, 0.3f, 10),
                new ClonotypeSearchResult(null, null, -1, 0.1f, 10),
                new ClonotypeSearchResult(null, null, -1, 0.3f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.15f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.01f, 1),
                new ClonotypeSearchResult(null, null, -1, 0.5f, 1)
        ]

        assert filter.filter(results).size() == 3
    }
}
