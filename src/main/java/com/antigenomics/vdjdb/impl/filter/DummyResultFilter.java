package com.antigenomics.vdjdb.impl.filter;

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult;

import java.util.List;

public final class DummyResultFilter implements ResultFilter {
    public static DummyResultFilter INSTANCE = new DummyResultFilter();

    @Override
    public List<ClonotypeSearchResult> filter(List<ClonotypeSearchResult> results) {
        return results;
    }
}
