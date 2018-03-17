package com.antigenomics.vdjdb.impl.filter;

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult;

import java.util.List;
import java.util.stream.Collectors;

public class TopNResultFilter extends ScoreThresholdResultFilter {
    private final int n;

    public TopNResultFilter(float threshold, int n) {
        super(threshold);
        this.n = n;
    }

    @Override
    public List<ClonotypeSearchResult> filter(List<ClonotypeSearchResult> results) {
        List<ClonotypeSearchResult> filteredResults = super.filter(results);

        filteredResults.sort((o1, o2) -> -Float.compare(o1.getWeightedScore(), o2.getWeightedScore()));

        return filteredResults.stream().limit(n).collect(Collectors.toList());
    }
}
