package com.antigenomics.vdjdb.web

/**
 * Created by mikesh on 7/23/17.
 */
class EpitopeSuggestion {
    String sequence // suggested sequence
    int substitutions, // subst diff from query
        indels, // insertion + deletion dist from query
        length, // sequence length
        count // count of all entries in VDJdb with this epitope

    EpitopeSuggestion(String sequence, int substitutions, int indels, int length, int count) {
        this.sequence = sequence
        this.substitutions = substitutions
        this.indels = indels
        this.length = length
        this.count = count
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        EpitopeSuggestion that = (EpitopeSuggestion) o

        if (count != that.count) return false
        if (indels != that.indels) return false
        if (length != that.length) return false
        if (substitutions != that.substitutions) return false
        if (sequence != that.sequence) return false

        return true
    }

    int hashCode() {
        int result
        result = sequence.hashCode()
        result = 31 * result + substitutions
        result = 31 * result + indels
        result = 31 * result + length
        result = 31 * result + count
        return result
    }
}
