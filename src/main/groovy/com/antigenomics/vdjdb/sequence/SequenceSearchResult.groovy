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

package com.antigenomics.vdjdb.sequence

import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

/**
 * Sequence alignment result 
 */
class SequenceSearchResult {
    /**
     * Amino acid sequence alignment, computed using sequence tree
     */
    final Alignment alignment
    /**
     * Score as reported by the scoring algorithm
     */
    final double score

    /**
     * Creates a new amino acid sequence alignment result
     * @param alignment amino acid sequence alignment, computed using sequence tree
     * @param score score as reported by the scoring algorithm
     */
    SequenceSearchResult(Alignment<AminoAcidSequence> alignment, double score) {
        this.alignment = alignment
        this.score = score
    }
}
