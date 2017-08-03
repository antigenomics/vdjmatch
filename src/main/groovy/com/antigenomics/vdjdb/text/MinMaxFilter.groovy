package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.db.Entry

/**
 * Created by mikesh on 7/23/17.
 */
class MinMaxFilter extends TextFilter {
    final int min
    final int max

    MinMaxFilter(String columnId, int min, int max) {
        super(columnId, "", false)
        this.min = min
        this.max = max
    }

    @Override
    protected boolean passInner(Entry entry) {
        return (entry.value.size() >= min) && (entry.value.size() <= max)
    }
}
