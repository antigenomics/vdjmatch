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

import com.milaboratory.core.tree.TreeSearchParameters

class SequenceSearcherPreset {
    final AlignmentScoring scoring
    final TreeSearchParameters parameters
    final boolean exhaustive, // if false stop with first variant from tree search, otherwise search all possible re-alignments
                  greedy // if true stop if new variant has more mismatches than the previous one

    final static SequenceSearcherPreset EXACT = new SequenceSearcherPreset()

    SequenceSearcherPreset(TreeSearchParameters parameters = new TreeSearchParameters(0, 0, 0, 0),
                           AlignmentScoring scoring = DummyAlignmentScoring.INSTANCE,
                           boolean exhaustive = false,
                           boolean greedy = true) {
        this.scoring = scoring
        this.parameters = parameters
        this.exhaustive = exhaustive
        this.greedy = greedy
    }
}
