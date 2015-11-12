package com.antigenomics.vdjdb.sequence

import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class SequenceSearchResult {
    final Alignment alignment
    final double penalty

    SequenceSearchResult(Alignment<AminoAcidSequence> alignment, double penalty) {
        this.alignment = alignment
        this.penalty = penalty
    }
}
