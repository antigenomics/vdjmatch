package com.antigenomics.vdjdb.fixer

class FixerResult {
    final String cdr3, vId, jId
    final boolean vCanonical, jCanonical
    final FixType vFixType, jFixType

    FixerResult(String cdr3, String vId, FixType vFixType, String jId, FixType jFixType) {
        this.cdr3 = cdr3
        this.vId = vId
        this.vFixType = vFixType
        this.jFixType = jFixType
        this.jId = jId
        this.vCanonical = cdr3.startsWith("C")
        this.jCanonical = cdr3.endsWith("F") || cdr3.endsWith("W")
    }

    static
    final String HEADER = "fixed.cdr3\tclosest.v.id\tclosest.j.id\tv.canonical\tj.canonical\tv.fix.type\tj.fix.type"

    @Override
    String toString() {
        [cdr3, vId, jId, vCanonical, jCanonical, vFixType, jFixType].join("\t")
    }
}
