package com.antigenomics.vdjdb.fixer

class FixerResult {
    final String cdr3
    final FixType vFixType, jFixType

    FixerResult(String cdr3, FixType vFixType, FixType jFixType) {
        this.cdr3 = cdr3
        this.vFixType = vFixType
        this.jFixType = jFixType
    }
}
