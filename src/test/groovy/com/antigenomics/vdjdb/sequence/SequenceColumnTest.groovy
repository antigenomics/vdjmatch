package com.antigenomics.vdjdb.sequence

import org.junit.Test

import static com.antigenomics.vdjdb.sequence.ExampleSequenceColumn.*

class SequenceColumnTest {
    @Test
    void exactSearchTest() {
        def filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF")
        assert SC.search(filter).size() == 1

        filter = new SequenceFilter("sc", "AASSLAPGATNEKLFF")
        assert SC.search(filter).size() == 0
    }

    @Test
    void scopeSearchTest1() {
        def filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(2, 1, 1, 2, false))

        assert SC.search(filter).size() == 7

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(2, 2, 1, 2, false))

        assert SC.search(filter).size() == 8

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(1, 1, 0, 1, false))

        assert SC.search(filter).size() == 3

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(1, 1, 1, 1, false))

        assert SC.search(filter).size() == 4
    }

    @Test
    void scopeSearchTest2() {
        def filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(1, 1, 1, 2))

        assert SC.search(filter).size() == 7

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(1, 2, 2))

        assert SC.search(filter).size() == 8


        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new SearchScope(1, 1, 2))

        assert SC.search(filter).size() == 5
    }

    @Test
    void exhaustiveTest() {
        // still better than working 3+ years in a row without a vacation...

        def filter = new SequenceFilter("sc", "CAWAAAGAAF",
                new SearchScope(1, 1, 1, 2, false),
                SM1AlignmentScoring.DEFAULT_BLOSUM62)

        def filterE = new SequenceFilter("sc", "CAWAAAGAAF",
                new SearchScope(1, 1, 1, 2, true),
                SM1AlignmentScoring.DEFAULT_BLOSUM62)

        assert SC.search(filter).values().first().alignmentScore < SC.search(filterE).values().first().alignmentScore
    }
}
