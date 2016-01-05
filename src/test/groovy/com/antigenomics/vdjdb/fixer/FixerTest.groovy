package com.antigenomics.vdjdb.fixer

import com.antigenomics.vdjdb.Util
import org.junit.Test

class FixerTest {
    static final Set<String> knownBadCases = new HashSet<>([".", "TRBV17", "TRBV13-1"])

    @Test
    void test() {
        def fixer = new Cdr3Fixer()

        boolean first = true
        Util.resourceAsStream("vdjdb_legacy.txt").splitEachLine('\t') { splitLine ->
            if (first) {
                first = false
                return
            }

            def v = splitLine[2], j = splitLine[3], species = splitLine[5]

            assert knownBadCases.contains(v) || fixer.getSegmentSeq(species, v) != null
            assert knownBadCases.contains(j) || fixer.getSegmentSeq(species, j) != null
        }
    }
}
