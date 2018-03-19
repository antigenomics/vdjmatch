package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.impl.ScoringProvider
import com.milaboratory.core.alignment.Aligner
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.alignment.BLASTMatrix
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class AlignmentScoringTest {
    @Test
    void test() {
        def seq1 = new AminoAcidSequence("CASTLAAAPSATNDKLFF"),
            seq2 = new AminoAcidSequence("CASSLAPGATNEKLFF")
        Alignment nwAln = Aligner.alignGlobalLinear(
                LinearGapAlignmentScoring.getAminoAcidBLASTScoring(BLASTMatrix.BLOSUM62, -12),
                seq1, seq2),
                  nwAln11 = Aligner.alignGlobalLinear(
                          LinearGapAlignmentScoring.getAminoAcidBLASTScoring(BLASTMatrix.BLOSUM62, -12),
                          seq1, seq1),
                  nwAln22 = Aligner.alignGlobalLinear(
                          LinearGapAlignmentScoring.getAminoAcidBLASTScoring(BLASTMatrix.BLOSUM62, -12),
                          seq2, seq2)

        // Check we have same mutations
        def filterE = new SequenceFilter("sc", "CASTLAAAPSATNDKLFF",
                new SearchScope(3, 2, 0, 5, true),
                SM2AlignmentScoring.DEFAULT_BLOSUM62)
        def hit = ExampleSequenceColumn.SC.search(filterE).values().first()
        assert hit.mutations == nwAln.absoluteMutations

        // Test score calculation between sequence column and NW aligner
        assert hit.alignmentScore == (nwAln.score - Math.max(nwAln11.score, nwAln22.score))
    }

    @Test
    void vdjamLoadTest() {
        ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt")
    }

    // todo: residueWiseMax = true test

    @Test
    void vdjamScoringTest() {
        def filterE = new SequenceFilter("sc", "CASSDWGASSYEQYF",
                new SearchScope(3, 2, 4, true),
                ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt"))
        assert Math.round(100 * ExampleSequenceColumn.SC.search(filterE).values().first().alignmentScore) == Math.round(100 * -3.0665207f)

        filterE = new SequenceFilter("sc", "CLVGEGDNYQLIW",
                new SearchScope(3, 0, 3, true),
                ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt"))
        assert Math.round(100 * ExampleSequenceColumn.SC.search(filterE).values().first().alignmentScore) == Math.round(100 * -7.1182404f)

        // this also tests exhaustive mode - two hits, best is insertion first

        // todo: not same hit reversed?
        filterE = new SequenceFilter("sc", "CLVGETNAGKSTF",
                new SearchScope(2, 1, 3, true),
                ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt"))
        def score1 = ExampleSequenceColumn.SC.search(filterE).values().find {
            it.mutations.size() > 0 // no exact match
        }.alignmentScore

        assert Math.round(100 * score1) == Math.round(100 * -6.5490856f)

        filterE = new SequenceFilter("sc", "CAVGAGTNAGKSTF",
                new SearchScope(2, 1, 3, true),
                ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt"))
        def score2 = ExampleSequenceColumn.SC.search(filterE).values().find {
            it.mutations.size() > 0 // no exact match
        }.alignmentScore

        assert Math.round(100 * score1) == Math.round(100 * score2)
    }
}
