package com.antigenomics.vdjdb.core2.text

import com.antigenomics.vdjdb.core2.db.Entry

class SubstringTextFilter extends TextFilter {
    SubstringTextFilter(String columnId, String value, boolean negative) {
        super(columnId, value, negative)
    }

    @Override
    protected boolean passInner(Entry entry) {
        entry.value.contains(value)
    }
}
