package com.antigenomics.vdjdb.sequence;

import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;

public interface AlignmentScoring {
    float computeScore(AminoAcidSequence query,
                       Mutations<AminoAcidSequence> mutations);

    ScoringType getScoringType();
}
