package com.antigenomics.vdjdb.core.text

import com.antigenomics.vdjdb.core.db.Entry

class ExactTextFilter extends TextFilter {
    ExactTextFilter(String columnId, String value, boolean negative) {
        super(columnId, value, negative)
    }

    @Override
    protected boolean passInner(Entry entry) {
        entry.value == value
    }
}
