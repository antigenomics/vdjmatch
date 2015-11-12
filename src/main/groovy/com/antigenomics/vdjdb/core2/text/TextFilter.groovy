package com.antigenomics.vdjdb.core2.text

import com.antigenomics.vdjdb.core2.db.Entry

abstract class TextFilter {
    final String value
    final boolean negative

    TextFilter(String value, boolean negative) {
        this.value = value
        this.negative = negative
    }

    protected abstract boolean passInner(Entry entry)

    public boolean pass(Entry entry) {
        negative ^ passInner(entry)
    }
}
