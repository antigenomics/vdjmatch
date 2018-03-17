package com.antigenomics.vdjdb.impl.filter;

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult;

import java.util.ArrayList;
import java.util.List;

public class MaxScoreResultFilter extends ScoreThresholdResultFilter {
    public MaxScoreResultFilter() {
        super(Float.NEGATIVE_INFINITY);
    }

    public MaxScoreResultFilter(float threshold) {
        super(threshold);
    }

    @Override
    public List<ClonotypeSearchResult> filter(List<ClonotypeSearchResult> results) {
        float maxScore = Float.NEGATIVE_INFINITY;
        for (ClonotypeSearchResult result : results) {
            maxScore = maxScore > result.getWeightedScore() ? maxScore : result.getWeightedScore();
        }

        List<ClonotypeSearchResult> filteredResults = new ArrayList<>();

        if (maxScore >= threshold) {
            for (ClonotypeSearchResult result : results) {
                if (result.getWeightedScore() == maxScore) {
                    filteredResults.add(result);
                }
            }
        }

        return filteredResults;
    }
}
