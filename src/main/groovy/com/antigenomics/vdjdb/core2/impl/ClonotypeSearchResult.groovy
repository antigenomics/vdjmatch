package com.antigenomics.vdjdb.core2.impl

import com.antigenomics.vdjdb.core2.db.Row
import com.antigenomics.vdjdb.core2.sequence.SequenceSearchResult

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
