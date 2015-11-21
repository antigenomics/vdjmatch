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


package com.antigenomics.vdjdb

import com.milaboratory.core.sequence.AminoAcidSequence

class Util {
    static InputStream resourceAsStream(String resourceName) {
        Util.class.classLoader.getResourceAsStream(resourceName)
    }

    static AminoAcidSequence convert(String aaSeq) {
        aaSeq = aaSeq.trim()
        if (aaSeq.length() > 0 && aaSeq =~ /^[FLSYCWPHQRIMTNKVADEG]+$/)
            return new AminoAcidSequence(aaSeq)

        System.err.println("Error converting '$aaSeq' to amino acid sequences, entry will be skipped from search")
        null
    }

    static String simplifySegmentName(String segmentName) {
        segmentName = segmentName.split(",")[0] // take best match
        segmentName.split("\\*")[0] // trim allele if present
    }
}
