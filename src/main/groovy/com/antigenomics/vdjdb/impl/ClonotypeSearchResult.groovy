package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.db.Row
import com.antigenomics.vdjdb.sequence.SequenceSearchResult

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
