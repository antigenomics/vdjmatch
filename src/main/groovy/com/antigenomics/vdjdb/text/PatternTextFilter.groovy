package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.db.Entry

import java.util.regex.Pattern

class PatternTextFilter extends TextFilter {
    final Pattern pattern

    PatternTextFilter(String columnId, String value, boolean negative) {
        super(columnId, value, negative)
        pattern = Pattern.compile(value)
    }

    @Override
    protected boolean passInner(Entry entry) {
        entry.value =~ pattern
    }
}
