package com.antigenomics.vdjdb.core.db

import com.antigenomics.vdjdb.core.sequence.SequenceSearchResult

class DatabaseSearchResult {
    final Row row
    final SequenceSearchResult[] sequenceSearchResults

    DatabaseSearchResult(Row row, SequenceSearchResult[] sequenceSearchResults) {
        this.row = row
        this.sequenceSearchResults = sequenceSearchResults
    }
}
