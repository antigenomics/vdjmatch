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

package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.impl.model.CloglogAggregateScoring
import com.antigenomics.vdjdb.impl.segment.PrecomputedSegmentScoring
import com.antigenomics.vdjdb.sequence.SubstitutionMatrixAlignmentScoring

import static com.milaboratory.core.sequence.AminoAcidSequence.ALPHABET

class ScoringProvider {
    // todo: also provide logit scoring
    static ScoringBundle loadScoringBundle(String species, String gene,
                                           String[] fileNames = ["score_coef.txt", "segm_score.txt", "vdjam.txt"],
                                           boolean fromResource = true) {
        def Cc1 = Float.NaN, Cc2 = Float.NaN, Cv = Float.NaN, Cj = Float.NaN, Cg = Float.NaN

        try {
            boolean firstLine = true
            (fromResource ? Util.resourceAsStream(fileNames[0]) :
                    new File(fileNames[0])).splitEachLine("\t") { splitLine ->
                if (firstLine) {
                    firstLine = false
                    // todo: use header to fetch column indices
                } else {
                    if (splitLine[0].equalsIgnoreCase(species) && splitLine[1].equalsIgnoreCase(gene)) {
                        (Cc1, Cc2, Cv, Cj, Cg) = splitLine[2..6].collect { it.toDouble() }
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error loading scoring ${fileNames[0]}: " +
                    e.toString())
        }

        new ScoringBundle(loadVdjamScoring(Cg, fileNames[2], fromResource),
                loadSegmentScoring(species, gene, fileNames[1], fromResource),
                new CloglogAggregateScoring(Cc1, Cc2, Cv, Cj))
    }

    static PrecomputedSegmentScoring loadSegmentScoring(String species, String gene,
                                                        String fileName = "segm_score.txt",
                                                        boolean fromResource = true) {
        Map<String, Map<String, float[]>> vScores = new HashMap<>()
        Map<String, Map<String, Float>> jScores = new HashMap<>()

        try {
            boolean firstLine = true
            (fromResource ? Util.resourceAsStream(fileName) :
                    new File(fileName)).splitEachLine("\t") { splitLine ->
                if (firstLine) {
                    firstLine = false
                    // todo: use header to fetch column indices
                } else {
                    if (splitLine[0].equalsIgnoreCase(species) && splitLine[1].equalsIgnoreCase(gene)) {
                        if (splitLine[2].equalsIgnoreCase("J")) {
                            def j1 = splitLine[3], j2 = splitLine[4],
                                score = (float) splitLine[7].toDouble()

                            def inner = jScores.getOrDefault(j1, new HashMap<String, Float>())
                            inner.put(j2, score)
                            jScores.put(j1, inner)
                        } else {
                            def v1 = splitLine[3], v2 = splitLine[4],
                                score1 = (float) splitLine[5].toDouble(),
                                score2 = (float) splitLine[6].toDouble(),
                                score3 = (float) splitLine[7].toDouble()

                            def inner = vScores.getOrDefault(v1, new HashMap<String, float[]>())
                            inner.put(v2, [score1, score2, score3] as float[])
                            vScores.put(v1, inner)
                        }
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error loading scoring $fileName: " +
                    e.toString())
        }

        return new PrecomputedSegmentScoring(vScores, jScores)
    }

    static SubstitutionMatrixAlignmentScoring loadVdjamScoring(float gapFactor = 0,
                                                               String fileName = "vdjam.txt",
                                                               boolean fromResource = true) {
        def substitutionMatrix = new float[ALPHABET.size()][ALPHABET.size()]

        try {
            boolean firstLine = true
            (fromResource ? Util.resourceAsStream(fileName) :
                    new File(fileName)).splitEachLine("\t") { splitLine ->
                if (firstLine) {
                    firstLine = false
                    // todo: use header to fetch column indices
                } else {
                    def from = splitLine[0], to = splitLine[1], score = (float) splitLine[2].toDouble()

                    def fromCode = ALPHABET.symbolToCode(from.charAt(0)),
                        toCode = ALPHABET.symbolToCode(to.charAt(0))
                    substitutionMatrix[fromCode][toCode] = score
                    substitutionMatrix[toCode][fromCode] = score
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error loading scoring $fileName: " +
                    e.toString())
        }

        return new SubstitutionMatrixAlignmentScoring(substitutionMatrix, gapFactor)
    }
}