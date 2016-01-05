package com.antigenomics.vdjdb.fixer

class OneSideFixerResult {
    final String cdr3
    final FixType fixType

    OneSideFixerResult(String cdr3, FixType fixType) {
        this.cdr3 = cdr3
        this.fixType = fixType
    }
}
