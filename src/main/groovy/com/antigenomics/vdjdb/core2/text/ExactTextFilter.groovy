package com.antigenomics.vdjdb.core2.text

import com.antigenomics.vdjdb.core2.db.Entry

class ExactTextFilter extends TextFilter {
    ExactTextFilter(String value, boolean negative) {
        super(value, negative)
    }

    @Override
    protected boolean passInner(Entry entry) {
        entry.value == value
    }
}
