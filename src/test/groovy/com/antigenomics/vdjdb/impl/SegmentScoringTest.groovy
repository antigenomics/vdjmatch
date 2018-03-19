package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.impl.segment.DummySegmentScoring
import com.antigenomics.vdjdb.impl.segment.PrecomputedSegmentScoring
import org.junit.Test

class SegmentScoringTest {
    @Test
    void dummyTest() {
        def scoring = DummySegmentScoring.INSTANCE
        def scores = scoring.computeScores("", "", "", "")
        assert scores.vScore == 0
        assert scores.jScore == 0
        assert scores.cdr1Score == 0
        assert scores.cdr2Score == 0
    }

    @Test
    void precompScoringTest() {
        def scoring = new PrecomputedSegmentScoring([("a"): [("b"): [1, 2, 3] as float[],
                                                             ("c"): [0, 0, 0] as float[]]],
                [("a"): [("b"): 4f, ("j"): 0f], ("e"): [("f"): 0f]])
        def scores = scoring.computeScores("a", "b", "a", "b")
        assert scores.vScore == 1f
        assert scores.jScore == 4f
        assert scores.cdr1Score == 2f
        assert scores.cdr2Score == 3f

        scores = scoring.computeScores("a", "d", "a", "d")
        assert scores.vScore == 1f / 2
        assert scores.jScore == 4f / 2
        assert scores.cdr1Score == 2f / 2
        assert scores.cdr2Score == 3f / 2

        scores = scoring.computeScores("a", "d", "d", "d")
        assert scores.jScore == (float) (4.0 / 3.0)
    }

    @Test
    void vdjmatchScoringTest() {
        def scoring = ScoringProvider.loadSegmentScoring("HomoSapiens", "TRB",
                "segm_score_legacy.txt")
        def scores = scoring.computeScores("TRBV7-2", "TRBV7-2*01", "TRBJ2-3", "TRBJ2-3*01")
        assert scores.vScore == 0
        assert scores.jScore == 0
        assert scores.cdr1Score == 0
        assert scores.cdr2Score == 0

        scores = scoring.computeScores("TRBV7-2", "TRBV7-3*01", "TRBJ2-3", "TRBJ2-4*01")
        assert scores.vScore == -102
        assert scores.jScore == -12
        assert scores.cdr1Score == 0
        assert scores.cdr2Score == -4.705516f
    }
}
