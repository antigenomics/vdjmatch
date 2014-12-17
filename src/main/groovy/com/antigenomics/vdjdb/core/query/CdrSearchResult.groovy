package com.antigenomics.vdjdb.core.query

import com.antigenomics.vdjdb.core.db.CdrEntrySet
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class CdrSearchResult implements Comparable<CdrSearchResult> {
    private final AminoAcidSequence query
    private final Alignment<AminoAcidSequence> alignment
    private final CdrEntrySet cdrEntrySet

    CdrSearchResult(AminoAcidSequence query,
                    Alignment<AminoAcidSequence> alignment,
                    CdrEntrySet cdrEntrySet) {
        this.query = query
        this.alignment = alignment
        this.cdrEntrySet = cdrEntrySet
    }

    @Override
    public int compareTo(CdrSearchResult o) {
        alignment.score.compareTo(o.alignment.score)
    }

    public Alignment<AminoAcidSequence> getAlignment() {
        alignment
    }

    public AminoAcidSequence getQuery() {
        query
    }

    public CdrEntrySet getCdrEntrySet() {
        cdrEntrySet
    }

    @Override
    public String toString() {
        "Alignment:\n" + alignment.getAlignmentHelper().toString() + "\nDatabase entries:\n" +
                "v\tj\t" + cdrEntrySet.parent.ANNOTATION_HEADER + "\n" +
                cdrEntrySet.collect { it.toString() }.join("\n")
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        CdrSearchResult that = (CdrSearchResult) o

        cdrEntrySet == that.cdrEntrySet && query == that.query
    }

    @Override
    int hashCode() {
        31 * query.hashCode() + cdrEntrySet.hashCode()
    }
}
