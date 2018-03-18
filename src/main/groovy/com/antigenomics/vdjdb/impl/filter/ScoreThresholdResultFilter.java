package com.antigenomics.vdjdb.impl.filter;

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult;

import java.util.ArrayList;
import java.util.List;

public class ScoreThresholdResultFilter implements ResultFilter {
    protected final float threshold;

    public ScoreThresholdResultFilter(float threshold) {
        this.threshold = threshold;
    }

    @Override
    public List<ClonotypeSearchResult> filter(List<ClonotypeSearchResult> results) {
        List<ClonotypeSearchResult> filteredResults = new ArrayList<>();

        for (ClonotypeSearchResult result : results) {
            if (result.getScore() >= threshold) {
                filteredResults.add(result);
            }
        }

        return filteredResults;
    }

    public float getThreshold() {
        return threshold;
    }
}
