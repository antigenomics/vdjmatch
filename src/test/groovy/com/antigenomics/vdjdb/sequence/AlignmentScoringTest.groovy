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
        def hit = SequenceColumnTest.SC.search(filterE).values().first()
        assert hit.mutations == nwAln.absoluteMutations

        // Test score calculation between sequence column and NW aligner
        assert hit.alignmentScore == (nwAln.score - Math.max(nwAln11.score, nwAln22.score))
    }

    @Test
    void vdjamLoadTest() {
        ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt")
    }

    @Test
    void vdjamScoringTest() {
        // todo: case CAVGAGTNAGKSTF CLVGETNAGKSTF

        def filterE = new SequenceFilter("sc", "CASSDWGASSYEQYF",
                new SearchScope(3, 2, 4, true),
                ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt"))
        assert SequenceColumnTest.SC.search(filterE).values().first().alignmentScore == -3.0665207f

        filterE = new SequenceFilter("sc", "CLVGEGDNYQLIW",
                new SearchScope(3, 0, 3, true),
                ScoringProvider.loadVdjamScoring(0, false, "vdjam_legacy.txt"))
        assert SequenceColumnTest.SC.search(filterE).values().first().alignmentScore == -7.1182404f
    }
}
