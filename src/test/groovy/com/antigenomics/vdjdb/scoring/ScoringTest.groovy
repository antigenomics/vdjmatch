package com.antigenomics.vdjdb.scoring

import com.milaboratory.core.alignment.Aligner
import com.milaboratory.core.alignment.BLASTMatrix
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidAlphabet
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class ScoringTest {
    @Test
    void scoringTest() {
        def scoring = SequenceSearcherPreset.byPrecision(0.8f).scoring

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

    @Test
    void ioTest() {
        def scoring1 = new VdjdbAlignmentScoring(new LinearGapAlignmentScoring(AminoAcidSequence.ALPHABET,
                1, -4, -4), [1.0f] * 11 as float[], 0),
            scoring2 = new VdjdbAlignmentScoring(new LinearGapAlignmentScoring(AminoAcidSequence.ALPHABET,
                    2, -4, -4), [1.0f] * 11 as float[], 1)

        def tempFileName = "temp_scoring_test.txt"

        AlignmentScoringProvider.saveScoring(["1": scoring1, "2": scoring2], tempFileName)

        def scoring11 = AlignmentScoringProvider.loadScoring("1", false, tempFileName)

        assert scoring11.scoring.getScore(AminoAcidAlphabet.H, AminoAcidAlphabet.H) == 1
        assert scoring11.scoring.getScore(AminoAcidAlphabet.H, AminoAcidAlphabet.G) == -4
        assert scoring11.scoreThreshold == 0

        def scoring22 = AlignmentScoringProvider.loadScoring("2", false, tempFileName)

        assert scoring22.scoring.getScore(AminoAcidAlphabet.H, AminoAcidAlphabet.H) == 2
        assert scoring22.scoring.getScore(AminoAcidAlphabet.H, AminoAcidAlphabet.G) == -4
        assert scoring22.scoreThreshold == 1
    }
}
