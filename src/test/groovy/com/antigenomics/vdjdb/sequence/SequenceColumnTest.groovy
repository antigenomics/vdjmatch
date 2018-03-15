package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.db.Database
import com.milaboratory.core.tree.TreeSearchParameters
import org.junit.Test

class SequenceColumnTest {
    static final SequenceColumn SC = new SequenceColumn("sc")
    static final Database DUMMY_DB = new Database([SC])
    static {
        DUMMY_DB.addEntries(
                [
                        // -
                        ["CASSLAPGAATNEKLFF"], // 1ins
                        ["CASSLAPGATNEKLFF"],
                        ["CASSLAPGAANEKLFF"],  // 1mm
                        ["CASSLAPGTNEKLFF"],   // 1del
                        ["CASSLAPGNNEKLFF"],   // 1del+1mm
                        ["CASSLPGATNAEKLFF"],  // 1del+1ins
                        ["ASSLAPGATNAEKLFF"],  // 1del+1ins
                        ["CASSLAPTNEKLFF"],    // 2del
                        // -
                        ["CAGAAAWAAF"],        // exhaustive test
                        ["CAAAAAAAF"],         // exhaustive test
                        // -
                        ["CASSDWGSYEQYF"],     // vdjam test1
                        ["CLVGDLTNYQLIW"],     // vdjam test2
                        ["CAVGAGTNAGKSTF"]     // vdjam test3
                ]
        )
    }

    @Test
    void exactSearchTest() {
        def filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF")
        assert SC.search(filter).size() == 1

        filter = new SequenceFilter("sc", "AASSLAPGATNEKLFF")
        assert SC.search(filter).size() == 0
    }

    @Test
    // FIXME
    void scopeSearchTest1() {
        def filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new TreeSearchParameters(1, 1, 1, 2),
                2)

        assert SC.search(filter).size() == 7

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new TreeSearchParameters(1, 2, 1, 2),
                2)

        assert SC.search(filter).size() == 8

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new TreeSearchParameters(1, 1, 0, 1))

        assert SC.search(filter).size() == 3

        filter = new SequenceFilter("sc", "CASSLAPGATNEKLFF",
                new TreeSearchParameters(1, 1, 1, 1))

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
