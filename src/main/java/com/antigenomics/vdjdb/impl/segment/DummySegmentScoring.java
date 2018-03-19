package com.antigenomics.vdjdb.impl.segment;

public class DummySegmentScoring implements SegmentScoring {
    public static final DummySegmentScoring INSTANCE = new DummySegmentScoring();

    private DummySegmentScoring() {

    }

    @Override
    public SegmentScores computeScores(String v1, String v2, String j1, String j2) {
        return new SegmentScores(0, 0, 0, 0);
    }
}
