package com.antigenomics.vdjdb.core2.sequence

import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.TreeSearchParameters

class SequenceSearchParameters {
    final AminoAcidSequence query
    final TreeSearchParameters treeSearchParameters
    final int depth

    SequenceSearchParameters(AminoAcidSequence query,
                             TreeSearchParameters treeSearchParameters,
                             int depth) {
        this.query = query
        this.treeSearchParameters = treeSearchParameters
        this.depth = depth
    }
}
