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

            assert knownBadCases.contains(v) || fixer.getClosestId(species, v) != null
            assert knownBadCases.contains(j) || fixer.getClosestId(species, j) != null
        }
    }

    @Test
    void testFixing() {
        def fixer = new Cdr3Fixer()

        boolean first = true
        def fixAttempted = 0, good = 0, fixed = 0, total = 0

        println FixerResult.HEADER

        Util.resourceAsStream("vdjdb_legacy.txt").splitEachLine('\t') { splitLine ->
            if (first) {
                first = false
                return
            }

            def cdr3 = splitLine[1], v = splitLine[2], j = splitLine[3], species = splitLine[5]

            def fixResult = fixer.fix(cdr3, v, j, species)

            assert fixResult.vCanonical || fixResult.vFixType == FixType.FailedBadSegment || fixResult.vFixType == FixType.FailedNoAlignment
            assert fixResult.jCanonical || fixResult.jFixType == FixType.FailedBadSegment || fixResult.jFixType == FixType.FailedNoAlignment

            [fixResult.vFixType, fixResult.jFixType].each {
                total++
                if (it.fixAttempted) {
                    fixAttempted++
                    if (it.good)
                        fixed++
                }
                if (it.good) {
                    good++
                } else {
                    println fixResult
                }
            }
        }

        println "Total=" + total
        println "Fix attempted=" + fixAttempted
        println "Fixed=" + fixed
        println "Good=" + good
    }

    @Test
    void testJ() {
        def cdr3 = "CASSPQTGTGGYGYTFG", id = "TRBJ1-2", species = "human"

        def fixer = new Cdr3Fixer()

        def segmentSeq = fixer.getSegmentSeq(species, fixer.getClosestId(species, id))

        def result = fixer.fix(cdr3.reverse(), segmentSeq.reverse())
        assert result.cdr3.reverse().endsWith("YTF")
        assert result.fixType == FixType.FixTrim

        cdr3 = "CASSPQTGTGGYGY"

        result = fixer.fix(cdr3.reverse(), segmentSeq.reverse())
        assert result.cdr3.reverse().endsWith("YTF")
        assert result.fixType == FixType.FixAdd
    }

    @Test
    void testV() {
        def cdr3 = "AACASRYRDDSYNEQFF", id = "TRBV7-9", species = "human"

        def fixer = new Cdr3Fixer()

        def segmentSeq = fixer.getSegmentSeq(species, fixer.getClosestId(species, id))

        def result = fixer.fix(cdr3, segmentSeq)
        assert result.cdr3.startsWith("CAS")
        assert result.fixType == FixType.FixTrim

        cdr3 = "ASSPQTGTGGYGY"

        result = fixer.fix(cdr3, segmentSeq)
        assert result.cdr3.startsWith("CAS")
        assert result.fixType == FixType.FixAdd
    }

    @Test
    void testVJ1() {
        def cdr3 = "CSVPFEGGTLETQH", vId = "TRBV29", jId = "TRBJ2-5", species = "human"

        def fixer = new Cdr3Fixer()

        def result = fixer.fix(cdr3, vId, jId, species)
        assert result.jFixType == FixType.FixReplace
        assert result.jCanonical
    }
}
