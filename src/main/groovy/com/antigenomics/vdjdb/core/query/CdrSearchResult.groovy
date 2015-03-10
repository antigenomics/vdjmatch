package com.antigenomics.vdjdb.core.query

import com.antigenomics.vdjdb.core.db.CdrEntrySet
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

/**
 * A match between query CDR3 amino acid sequence and a subject sequence from the database.
 * Note that the subject sequence can represent multiple database records.
 * {@see com.antigenomics.vdjdb.core.db.CdrEntrySet}
 */
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

    /**
     * Gets the alignment of query sequence ({@link com.milaboratory.core.alignment.Alignment#getSequence1()}) and
     * a subject sequence from the database. 
     * @return local alignment.
     */
    Alignment<AminoAcidSequence> getAlignment() {
        alignment
    }

    /**
     * Gets the query CDR3 amino acid sequence. 
     * @return CDR3 amino acid sequence.
     */
    AminoAcidSequence getQuery() {
        query
    }

    /**
     * Gets a database entry set that correspond to a subject sequence.
     * @return database entries.
     */
    CdrEntrySet getCdrEntrySet() {
        cdrEntrySet
    }

    @Override
    String toString() {
        "Alignment:\n" + alignment.getAlignmentHelper().toString() + "\nDatabase entries:\n" +
                "v\tj\t" + cdrEntrySet.parent.ANNOTATION_HEADER + "\n" +
                cdrEntrySet.collect { it.toString() }.join("\n")
    }

    @Override
    int compareTo(CdrSearchResult o) {
        alignment.score.compareTo(o.alignment.score)
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
