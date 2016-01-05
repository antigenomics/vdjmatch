/*
 * Copyright 2015 Mikhail Shugay
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


package com.antigenomics.vdjdb.fixer

import com.antigenomics.vdjdb.Util

class Cdr3Fixer {
    static final Map<String, String> speciesAliases = ['HomoSapiens': 'human',
                                                       'MusMusculus': 'mouse']

    final Map<String, Map<String, String>> segmentsByIdBySpecies = new HashMap<>()

    Cdr3Fixer() {
        speciesAliases.values().each {
            segmentsByIdBySpecies.put(it, new HashMap<String, String>())
        }
        Util.resourceAsStream("segments.txt").splitEachLine('\t') { splitLine ->
            if (!splitLine[0].startsWith("#")) {
                def species = speciesAliases[splitLine[0]],
                    type = splitLine[2].toLowerCase(), id = splitLine[3],
                    refPoint = splitLine[4].toInteger(), seq = splitLine[5]
                if (species != null && (type.startsWith("v") || type.startsWith("j"))) {

                    seq = type.startsWith("v") ? seq.substring(refPoint - 3) : seq.substring(0, refPoint + 4)

                    segmentsByIdBySpecies[species].put(id, seq)
                }
            }
        }
    }

    String getSegmentSeq(String species, String id) {
        def segmentsById = segmentsByIdBySpecies[species]

        if (segmentsById == null)
            return null

        segmentsById[id] ?: segmentsById[(id + "*01").toString()] ?: segmentsById[(id + "-1*01").toString()]
    }
}
