package com.antigenomics.vdjdb.sequence

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
                LinearGapAlignmentScoring.getAminoAcidBLASTScoring(BLASTMatrix.BLOSUM62),
                seq1, seq2)

        //println nwAln
        //println nwAln.score

        def filterE = new SequenceFilter("sc", "CASTLAAAPSATNDKLFF",
                new SearchScope(3, 2, 0, 5, true),
                SubstitutionMatrixAlignmentScoring.DEFAULT_BLOSUM62)

        assert SequenceColumnTest.SC.search(filterE).values().first().mutations == nwAln.absoluteMutations
    }
}
