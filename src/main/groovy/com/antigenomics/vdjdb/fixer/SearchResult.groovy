package com.antigenomics.vdjdb.fixer

class SearchResult {
    final int startInSegment, startInCdr3, matchSize

    SearchResult(int startInSegment, int startInCdr3, int matchSize) {
        this.startInSegment = startInSegment
        this.startInCdr3 = startInCdr3
        this.matchSize = matchSize
    }
}
