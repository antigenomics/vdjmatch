package com.antigenomics.vdjdb.fixer

import com.antigenomics.vdjdb.Util
import org.junit.Test

class FixerTest {
    static final Set<String> knownBadCases = new HashSet<>([".", "TRBV17", "TRBV13-1"])

    @Test
    void testConsistency() {
        def fixer = new Cdr3Fixer()

        fixer.segmentsByIdBySpecies.values().each { segmentsById ->
            segmentsById.each {
                if (it.key.contains("V")) {
                    assert it.value.startsWith("C")
                } else if (it.key.contains("J")) {
                    assert it.value.endsWith("F") || it.value.endsWith("W")
                }
            }
        }
    }

    @Test
    void testPresence() {
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

    @Test
    void testJ() {
        def cdr3 = "CASSPQTGTGGYGYTFG", id = "TRBJ1-2", species = "human"

        def fixer = new Cdr3Fixer()

        def segmentSeq = fixer.getSegmentSeq(species, id)

        def result = Cdr3Fixer.fix(cdr3.reverse(), segmentSeq.reverse())
        assert result.cdr3.reverse().endsWith("YTF")
        assert result.fixType == FixType.FixTrim

        cdr3 = "CASSPQTGTGGYGY"

        result = Cdr3Fixer.fix(cdr3.reverse(), segmentSeq.reverse())
        assert result.cdr3.reverse().endsWith("YTF")
        assert result.fixType == FixType.FixAdd
    }

    @Test
    void testV() {
        def cdr3 = "AACASRYRDDSYNEQFF", id = "TRBV7-9", species = "human"

        def fixer = new Cdr3Fixer()

        def segmentSeq = fixer.getSegmentSeq(species, id)

        def result = Cdr3Fixer.fix(cdr3, segmentSeq)
        assert result.cdr3.startsWith("CAS")
        assert result.fixType == FixType.FixTrim

        cdr3 = "ASSPQTGTGGYGY"

        result = Cdr3Fixer.fix(cdr3, segmentSeq)
        assert result.cdr3.startsWith("CAS")
        assert result.fixType == FixType.FixAdd
    }
}
