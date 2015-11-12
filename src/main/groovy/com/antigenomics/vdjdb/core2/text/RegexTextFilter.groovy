package com.antigenomics.vdjdb.core2.text

import com.antigenomics.vdjdb.core2.db.Entry

import java.util.regex.Pattern

class RegexTextFilter extends TextFilter {
    final Pattern pattern

    RegexTextFilter(String columnId, String value, boolean negative) {
        super(columnId, value, negative)
        pattern = Pattern.compile(value)
    }

    @Override
    protected boolean passInner(Entry entry) {
        entry.value =~ pattern
    }
}
