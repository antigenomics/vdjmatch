package com.antigenomics.vdjdb.scoring

import com.antigenomics.vdjdb.sequence.AlignmentScoringProvider
import com.milaboratory.core.alignment.Aligner
import com.milaboratory.core.alignment.BLASTMatrix
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class ScoringTest {
    @Test
    void loadScoringTest() {
        AlignmentScoringProvider.loadScoring() // asserts inside
    }

    @Test
    void scoringTest() {
        def scoring = AlignmentScoringProvider.loadScoring()

        def seq1 = new AminoAcidSequence("CASSLAPGATNEKLFF"),
            seq2 = new AminoAcidSequence("CASSLAPGATNEKLFF")

        def score = {
            def aln = Aligner.alignGlobalLinear(LinearGapAlignmentScoring.getAminoAcidBLASTScoring(BLASTMatrix.BLOSUM45),
                    seq1, seq2)
            scoring.computeScore(aln)
        }

        assert score() >= scoring.scoreThreshold

        seq2 = new AminoAcidSequence("CASSLAPATNEKLFF")

        assert score() >= scoring.scoreThreshold

        seq2 = new AminoAcidSequence("CASSLAPGATNEFF")

        assert score() >= scoring.scoreThreshold

        seq2 = new AminoAcidSequence("CASSVRLNTGEKLFF")

        assert score() < scoring.scoreThreshold

        seq2 = new AminoAcidSequence("CASSLAGATEKLFF")

        assert score() < scoring.scoreThreshold
    }
}
