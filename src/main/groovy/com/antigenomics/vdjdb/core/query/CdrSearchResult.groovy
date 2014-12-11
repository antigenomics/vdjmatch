package com.antigenomics.vdjdb.core.query

import com.antigenomics.vdjdb.core.db.CdrEntrySet
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class CdrSearchResult implements Comparable<CdrSearchResult> {
    private final Alignment<AminoAcidSequence> alignment
    private final CdrEntrySet cdrEntrySet

    CdrSearchResult(Alignment<AminoAcidSequence> alignment, CdrEntrySet cdrEntrySet) {
        this.alignment = alignment
        this.cdrEntrySet = cdrEntrySet
    }

    @Override
    public int compareTo(CdrSearchResult o) {
        alignment.score.compareTo(o.alignment.score)
    }

    public Alignment<AminoAcidSequence> getAlignment() {
        return alignment
    }

    public CdrEntrySet getCdrEntrySet() {
        return cdrEntrySet
    }

    @Override
    public String toString() {
        "Alignment:\n" + alignment.getAlignmentHelper().toString() + "\nDatabase hits:\n" +
                "v\tj\t" + cdrEntrySet.parent.ANNOTATION_HEADER + "\n" +
                cdrEntrySet.collect { it.toString() }.join("\n")
    }
}
