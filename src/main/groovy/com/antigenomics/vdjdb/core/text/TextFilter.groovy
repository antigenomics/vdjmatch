package com.antigenomics.vdjdb.core.text

import com.antigenomics.vdjdb.core.db.ColumnType
import com.antigenomics.vdjdb.core.db.Entry
import com.antigenomics.vdjdb.core.db.Filter

abstract class TextFilter implements Filter {
    final ColumnType columnType = ColumnType.Text
    final String columnId, value
    final boolean negative

    TextFilter(String columnId, String value, boolean negative) {
        this.columnId = columnId
        this.value = value
        this.negative = negative
    }

    protected abstract boolean passInner(Entry entry)

    public boolean pass(Entry entry) {
        negative ^ passInner(entry)
    }
}
