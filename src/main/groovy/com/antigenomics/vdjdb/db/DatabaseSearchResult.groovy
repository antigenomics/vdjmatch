package com.antigenomics.vdjdb.db

import com.antigenomics.vdjdb.sequence.SequenceSearchResult

class DatabaseSearchResult {
    final Row row
    final SequenceSearchResult[] sequenceSearchResults

    DatabaseSearchResult(Row row, SequenceSearchResult[] sequenceSearchResults) {
        this.row = row
        this.sequenceSearchResults = sequenceSearchResults
    }
}
