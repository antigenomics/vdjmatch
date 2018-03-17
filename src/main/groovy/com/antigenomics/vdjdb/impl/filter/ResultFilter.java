package com.antigenomics.vdjdb.impl.filter;

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult;

import java.util.List;

public interface ResultFilter {
    List<ClonotypeSearchResult> filter(List<ClonotypeSearchResult> results);
}
