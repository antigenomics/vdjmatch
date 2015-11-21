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
     * Amino acid sequence alignment
     */
    final Alignment alignment
    /**
     * Penalty as computed by the alignment algorithm 
     */
    final double penalty

    /**
     * Creates a new amino acid sequence alignment result
     * @param alignment amino acid sequence alignment
     * @param penalty penalty reported by the aligner
     */
    SequenceSearchResult(Alignment<AminoAcidSequence> alignment, double penalty) {
        this.alignment = alignment
        this.penalty = penalty
    }
}
