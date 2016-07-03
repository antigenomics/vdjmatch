package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.Util
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidAlphabet

import static com.milaboratory.core.sequence.AminoAcidSequence.ALPHABET

class AlignmentScoringProvider {
    static AlignmentScoring loadScoring(String scoringId = "default") {
        scoringId = scoringId.toLowerCase()
        def lines = Util.resourceAsStream("scoring.txt").readLines()

        def header = lines[0].toLowerCase().split("\t")

        def headerIndices = ["id", "parameter", "from", "to", "value"].collectEntries { colName ->
            [(colName): header.findIndexOf { colName.equalsIgnoreCase(it) }]
        }

        assert headerIndices.values().every { it >= 0 }

        int idCol = headerIndices["id"],
            parameterCol = headerIndices["parameter"],
            fromCol = headerIndices["from"],
            toCol = headerIndices["to"],
            valueCol = headerIndices["value"]

        lines = lines[1..-1].collect { it.split("\t") }.findAll {
            it[idCol].toLowerCase() == scoringId
        }

        if (lines.empty) {
            throw new RuntimeException("Failed to load scoring '$scoringId'")
        }

        def scoringMatrix = [], positionWeights = [], gapScore = Integer.MAX_VALUE, threshold = Integer.MAX_VALUE

        int n = ALPHABET.size()

        lines.each { splitLine ->
            switch (splitLine[parameterCol].toLowerCase()) {
                case "substitution":
                    scoringMatrix[ALPHABET.codeFromSymbol(splitLine[fromCol].charAt(0)) * n +
                            ALPHABET.codeFromSymbol(splitLine[toCol].charAt(0))] = splitLine[valueCol].toInteger()
                    break
                case "gap":
                    gapScore = splitLine[valueCol].toInteger()
                    break
                case "threshold":
                    threshold = splitLine[valueCol].toDouble()
                    break
                case "position_weight":
                    positionWeights[splitLine[fromCol].toInteger()] = splitLine[valueCol].toDouble()
                    break
            }
        }

        // Add zeros for incomplete/stop codons

        (0..<n).each { aa1 ->
            (0..<n).each { aa2 ->
                def i = (byte) aa1, j = (byte) aa2
                if (i == AminoAcidAlphabet.IncompleteCodon || i == AminoAcidAlphabet.Stop ||
                        j == AminoAcidAlphabet.IncompleteCodon || j == AminoAcidAlphabet.Stop) {
                    scoringMatrix[i * n + j] = 0
                }
            }
        }

        // Check that we've loaded everything

        assert scoringMatrix.every { it != null } && scoringMatrix.size() == n * n
        assert positionWeights.every { it != null } &&
                positionWeights.size() % 2 != 0
        assert gapScore != Integer.MAX_VALUE
        assert threshold != Integer.MAX_VALUE

        new AlignmentScoring(new LinearGapAlignmentScoring(ALPHABET, scoringMatrix as int[],
                gapScore), positionWeights as double[], threshold)
    }
}
