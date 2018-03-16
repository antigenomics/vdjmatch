package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.db.Database

class ExampleSequenceColumn {
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
                        ["CAVGAGTNAGKSTF"],    // vdjam test3-1
                        ["CLVGETNAGKSTF"],     // vdjam test3-2
                        // -
                        ["CASSLGANTIYF"]       // aggregate scoring test
                ]
        )
    }
}
