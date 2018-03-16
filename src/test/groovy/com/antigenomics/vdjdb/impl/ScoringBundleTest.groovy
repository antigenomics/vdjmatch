package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.impl.model.LinearAggregateScoring
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjdb.sequence.SequenceFilter

import static com.antigenomics.vdjdb.sequence.ExampleSequenceColumn.*
import org.junit.Test

class ScoringBundleTest {
    @Test
    void loadTest() {
        def scoringBundle = ScoringProvider.loadScoringBundle("HomoSapiens",
                "TRB", false,
                ["score_coef_legacy.txt", "segm_score_legacy.txt", "vdjam_legacy.txt"] as String[])
    }

    @Test
    void aggregateScoringTest() {
        def scoringBundle = ScoringProvider.loadScoringBundle("HomoSapiens",
                "TRB", false,
                ["score_coef_legacy.txt", "segm_score_legacy.txt", "vdjam_legacy.txt"] as String[])

        def segmentScores = scoringBundle.segmentScoring.computeScores("TRBV12-3*01",
                "TRBV28*01",
                "TRBJ1-3*01",
                "TRBJ2-5*01")

        //println segmentScores

        //def hit = ExampleSequenceColumn.SC

        // TRBV12-3*01 TRBJ1-3*01 CASSLGANTIYF
        // TRBV28*01 TRBJ2-5*01 CASSLGRETQYF
        // 0.06058031

        def filter = new SequenceFilter("sc", "CASSLGRETQYF",
                new SearchScope(3, 0, 0, 3),
                scoringBundle.alignmentScoring)

        def alignmentScore = SC.search(filter).values().first().alignmentScore


        assert scoringBundle.aggregateScoring.computeFullScore(segmentScores.vScore,
                segmentScores.cdr1Score, segmentScores.cdr2Score, alignmentScore, segmentScores.jScore) == 0.060580313f
    }
}
