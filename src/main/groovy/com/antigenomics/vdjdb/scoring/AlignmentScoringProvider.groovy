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
import com.milaboratory.core.sequence.AminoAcidSequence
import sun.reflect.generics.reflectiveObjects.NotImplementedException

import static com.milaboratory.core.sequence.AminoAcidSequence.ALPHABET

class AlignmentScoringProvider {
    static final String DUMMY_SCORING_ID = "dummy",
            GAP_CHARACTER = "_"

    static AlignmentScoring loadScoring(String prefix,
                                        boolean fromResource = true) {
        if (prefix.equalsIgnoreCase(DUMMY_SCORING_ID)) {
            return DummyAlignmentScoring.INSTANCE
        }

        int n = ALPHABET.size()

        def substitutionMatrix = new float[n][n], gapPenalties = new float[n]

        (0..<n).each { gapPenalties[it] = -1.0f }

        def scoringFileName = prefix + ".scoring.txt"

        try {
            (fromResource ? Util.resourceAsStream(scoringFileName) :
                    new File(scoringFileName)).splitEachLine("\t") { splitLine ->
                def from = splitLine[0], to = splitLine[1], score = (float) splitLine[2].toDouble()

                if (from == GAP_CHARACTER) {
                    gapPenalties[ALPHABET.symbolToCode(to.charAt(0))] = score
                } else if (to == GAP_CHARACTER) {
                    gapPenalties[ALPHABET.symbolToCode(from.charAt(0))] = score
                } else {
                    def fromCode = ALPHABET.symbolToCode(from.charAt(0)),
                        toCode = ALPHABET.symbolToCode(to.charAt(0))
                    substitutionMatrix[fromCode][toCode] = score
                    substitutionMatrix[toCode][fromCode] = score
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error loading scoring $scoringFileName: " +
                    e.toString())
        }

        new VdjdbAlignmentScoring(substitutionMatrix, gapPenalties)
    }
}
