/*
 * Copyright 2016 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjdb.scoring

import com.antigenomics.vdjdb.Util
import com.milaboratory.core.alignment.LinearGapAlignmentScoring

import static com.milaboratory.core.sequence.AminoAcidSequence.ALPHABET

class AlignmentScoringProvider {
    static final int N = ALPHABET.size()
    static final String ID_COL = "id", PARAMETER_COL = "parameter",
                        FIRST_ARG_COL = "from", SECOND_ARG_COL = "to", VALUE_COL = "value"
    static final List<String> HEADER_COLS = [ID_COL, PARAMETER_COL, FIRST_ARG_COL, SECOND_ARG_COL, VALUE_COL]

    static saveScoring(Map<String, VdjdbAlignmentScoring> scoringsById, String fileName) {
        new File(fileName).withPrintWriter { pw ->
            pw.println(HEADER_COLS.join("\t"))

            scoringsById.each {
                def index = it.key, info = it.value

                (0..<N).each { aa1 ->
                    (0..<N).each { aa2 ->
                        def i = (byte) aa1, j = (byte) aa2
                        if (i != ALPHABET.IncompleteCodon && i != ALPHABET.Stop &&
                                j != ALPHABET.IncompleteCodon && j != ALPHABET.Stop) {
                            int score = info.scoring.getScore(i, j)

                            pw.println(index + "\tsubstitution\t" +
                                    ALPHABET.symbolFromCode(i) + "\t" +
                                    ALPHABET.symbolFromCode(j) + "\t" +
                                    score)
                        }
                    }
                }

                pw.println(index + "\tgap\tNA\tNA\t" + info.scoring.gapPenalty)
                pw.println(index + "\tthreshold\tNA\tNA\t" + info.scoreThreshold)
                info.positionWeights.eachWithIndex { double value, int i ->
                    pw.println(index + "\tposition_weight\t" + i + "\tNA\t" + value)
                }
            }
        }
    }

    static VdjdbAlignmentScoring loadScoring(String scoringId,
                                             boolean fromResource = true,
                                             String fileName = "solutions.txt") {
        def lines = (fromResource ? Util.resourceAsStream(fileName) : new File(fileName)).readLines()

        def header = lines[0].toLowerCase().split("\t")

        def headerIndices = HEADER_COLS.collectEntries { colName ->
            [(colName): header.findIndexOf { colName.equalsIgnoreCase(it) }]
        }

        if ( headerIndices.values().any { it < 0 }) {
            throw new RuntimeException("Failed to parse scoring file $fileName, critical columns are missing")
        }

        int idCol = headerIndices[ID_COL],
            parameterCol = headerIndices[PARAMETER_COL],
            fromCol = headerIndices[FIRST_ARG_COL],
            toCol = headerIndices[SECOND_ARG_COL],
            valueCol = headerIndices[VALUE_COL]

        scoringId = scoringId.toLowerCase()
        lines = lines[1..-1].collect { it.split("\t") }.findAll {
            it[idCol].toLowerCase() == scoringId
        }

        if (lines.empty) {
            throw new RuntimeException("Scoring with id '$scoringId' is missing in $fileName")
        }

        def scoringMatrix = [], positionWeights = [], gapScore = Integer.MAX_VALUE, threshold = Float.MAX_VALUE

        lines.each { splitLine ->
            switch (splitLine[parameterCol].toLowerCase()) {
                case "substitution":
                    scoringMatrix[ALPHABET.codeFromSymbol(splitLine[fromCol].charAt(0)) * N +
                            ALPHABET.codeFromSymbol(splitLine[toCol].charAt(0))] = splitLine[valueCol].toInteger()
                    break
                case "gap":
                    gapScore = splitLine[valueCol].toInteger()
                    break
                case "threshold":
                    threshold = splitLine[valueCol].toFloat()
                    break
                case "position_weight":
                    positionWeights[splitLine[fromCol].toInteger()] = splitLine[valueCol].toFloat()
                    break
            }
        }

        // Add zeros for incomplete/stop codons

        (0..<N).each { aa1 ->
            (0..<N).each { aa2 ->
                def i = (byte) aa1, j = (byte) aa2
                if (i == ALPHABET.IncompleteCodon || i == ALPHABET.Stop ||
                        j == ALPHABET.IncompleteCodon || j == ALPHABET.Stop) {
                    scoringMatrix[i * N + j] = 0
                }
            }
        }

        // Check that we've loaded everything

        assert scoringMatrix.every { it != null } && scoringMatrix.size() == N * N
        assert positionWeights.every { it != null } &&
                positionWeights.size() % 2 != 0
        assert gapScore < 0
        assert threshold != Float.MAX_VALUE

        new VdjdbAlignmentScoring(new LinearGapAlignmentScoring(ALPHABET, scoringMatrix as int[],
                gapScore), positionWeights as float[], threshold)
    }
}
