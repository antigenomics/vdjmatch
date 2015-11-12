package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.db.Entry
import com.antigenomics.vdjdb.db.Filter

abstract class TextFilter implements Filter {
    final String columnId, value
    final boolean negative

    TextFilter(String columnId, String value, boolean negative) {
        this.columnId = columnId
        this.value = value
        this.negative = negative
    }

    protected abstract boolean passInner(Entry entry)

    boolean pass(Entry entry) {
        negative ^ passInner(entry)
    }

    @Override
    boolean isSequenceFilter() {
        false
    }
}
