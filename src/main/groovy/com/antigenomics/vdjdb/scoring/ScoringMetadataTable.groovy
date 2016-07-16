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

class ScoringMetadataTable {
    final List<ScoringMetadata> rows = new ArrayList<>()

    ScoringMetadataTable(String fileName = "roc.txt", boolean fromResource = true) {
        def lines = (fromResource ? Util.resourceAsStream(fileName) : new FileInputStream(fileName)).readLines()
        def iter = lines.iterator()
        def header = iter.next().split("\t")
        def idCol = header.findIndexOf { it.equalsIgnoreCase("id") },
            precisionCol = header.findIndexOf { it.equalsIgnoreCase("precision") },
            recallCol = header.findIndexOf { it.equalsIgnoreCase("recall") }

        if ([idCol, precisionCol, recallCol].any { it < 0 }) {
            throw new RuntimeException("Unable to parse scoring metadata from $fileName, critical columns are missing.")
        }

        while (iter.hasNext()) {
            def splitLine = iter.next().split("\t")
            // abs here: MOEA framework minimizes so precision/recall can be negative of corresponding values
            rows.add(new ScoringMetadata(splitLine[idCol],
                    Math.abs(splitLine[precisionCol].toFloat()),
                    Math.abs(splitLine[recallCol].toFloat())))
        }

        if (rows.empty) {
            throw new RuntimeException("No scoring metadata was loaded from $fileName")
        }
    }

    ScoringMetadata getOptimal() {
        rows.max { it.FScore }
    }

    ScoringMetadata getByPrecision(float precision) {
        if (precision < 0 || precision > 1) {
            throw new IllegalArgumentException()
        }

        rows.min { Math.abs(it.precision - precision) }
    }

    ScoringMetadata getByRecall(float recall) {
        if (recall < 0 || recall > 1) {
            throw new IllegalArgumentException()
        }

        rows.min { Math.abs(it.recall - recall) }
    }
}
