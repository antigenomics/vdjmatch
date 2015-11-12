package com.antigenomics.vdjdb.core2.db

import com.antigenomics.vdjdb.core2.sequence.SequenceSearchResult

class DatabaseSearchResult {
    final Row row
    final SequenceSearchResult[] sequenceSearchResults

    DatabaseSearchResult(Row row, SequenceSearchResult[] sequenceSearchResults) {
        this.row = row
        this.sequenceSearchResults = sequenceSearchResults
    }
}
