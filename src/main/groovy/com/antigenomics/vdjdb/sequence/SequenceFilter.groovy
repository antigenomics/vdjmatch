package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.db.Filter
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.TreeSearchParameters

class SequenceFilter implements Filter {
    final String columnId
    final AminoAcidSequence query
    final TreeSearchParameters treeSearchParameters
    final int depth

    SequenceFilter(String columnId, String query) {
        this(columnId, query, new TreeSearchParameters(2, 1, 1, 2))
    }

    SequenceFilter(String columnId, String query,
                   TreeSearchParameters treeSearchParameters) {
        this(columnId, Util.convert(query), treeSearchParameters, -1)
    }

    SequenceFilter(String columnId, String query,
                   TreeSearchParameters treeSearchParameters,
                   int depth) {
        this(columnId, Util.convert(query), treeSearchParameters, depth)
    }

    SequenceFilter(String columnId, AminoAcidSequence query,
                   TreeSearchParameters treeSearchParameters,
                   int depth) {
        if (query == null)
            throw new RuntimeException("Bad sequence filter query")
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
