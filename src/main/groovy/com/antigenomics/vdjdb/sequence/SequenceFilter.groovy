package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.db.Filter
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.TreeSearchParameters

class SequenceFilter implements Filter {
    final String columnId
    final AminoAcidSequence query
    final TreeSearchParameters treeSearchParameters
    final int depth

    SequenceFilter(String columnId, AminoAcidSequence query,
                   TreeSearchParameters treeSearchParameters,
                   int depth) {
        this.columnId = columnId
        this.query = query
        this.treeSearchParameters = treeSearchParameters
        this.depth = depth
    }

    @Override
    boolean isSequenceFilter() {
        true
    }
}
