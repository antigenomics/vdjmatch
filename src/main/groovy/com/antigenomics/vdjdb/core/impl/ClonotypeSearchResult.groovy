package com.antigenomics.vdjdb.core.impl

import com.antigenomics.vdjdb.core.db.Row
import com.antigenomics.vdjdb.core.sequence.SequenceSearchResult

class ClonotypeSearchResult implements Comparable<ClonotypeSearchResult> {
    final SequenceSearchResult result
    final Row row

    ClonotypeSearchResult(SequenceSearchResult result, Row row) {
        this.result = result
        this.row = row
    }

    @Override
    int compareTo(ClonotypeSearchResult o) {
        -result.penalty.compareTo(o.result.penalty)
    }
}
