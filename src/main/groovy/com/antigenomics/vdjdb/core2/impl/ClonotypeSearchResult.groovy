package com.antigenomics.vdjdb.core2.impl

import com.antigenomics.vdjdb.core2.db.Row
import com.antigenomics.vdjdb.core2.sequence.SequenceSearchResult
import com.antigenomics.vdjtools.sample.Clonotype

class ClonotypeSearchResult implements Comparable<ClonotypeSearchResult> {
    final Clonotype clonotype
    final SequenceSearchResult result
    final Row row

    ClonotypeSearchResult(Clonotype clonotype, SequenceSearchResult result, Row row) {
        this.clonotype = clonotype
        this.result = result
        this.row = row
    }

    @Override
    int compareTo(ClonotypeSearchResult o) {
        -result.penalty.compareTo(o.result.penalty)
    }
}
