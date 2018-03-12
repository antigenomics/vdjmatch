package com.antigenomics.vdjdb.sequence;

import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class DummyAlignmentScoring implements AlignmentScoring {
    public static final DummyAlignmentScoring INSTANCE = new DummyAlignmentScoring();

    private DummyAlignmentScoring() {

    }

    @Override
    public float computeScore(AminoAcidSequence query,
                              Mutations<AminoAcidSequence> mutations) {
        return -mutations.size();
    }
}
