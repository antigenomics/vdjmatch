package com.antigenomics.vdjdb.impl.segment;

public interface SegmentScoring {
    SegmentScores computeScores(String v1, String v2,
                                String j1, String j2);
}
