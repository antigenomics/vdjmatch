package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.db.Entry

/**
 * Created by mikesh on 7/23/17.
 */
class LengthFilter extends TextFilter {
    final int length

    LengthFilter(String columnId, int length) {
        super(columnId, "", false)
        this.length = length
    }

    @Override
    protected boolean passInner(Entry entry) {
        return entry.value.size() == length
    }
}
