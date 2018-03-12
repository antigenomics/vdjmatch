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

package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.impl.ScoringProvider
import com.milaboratory.core.tree.TreeSearchParameters

/**
 * Stores parameters for SequenceFilter specifying 'search scope' (number of allowed mismatches) and alignment
 * scoring method.
 */
class SearchScope {
    final int maxIndels
    final TreeSearchParameters parameters
    final boolean exhaustive, greedy

    final static SearchScope EXACT = new SearchScope(0, 0, 0, 0)

    /**
     * Create sequence search filter parameter set. Specifies maximal allowed number of substitutions, indels
     * (i.e. sum of insertions and deletions, or difference in query and target sequence lengths) and
     * Levenstein distance from query sequence. Note that resulting search is 'symmetric' as
     * number of insertions in query is the number of deletions in target. Used in 'Annotate' tab of VDJdb web app.
     *
     * @param maxSubstitutions maximal number of substitutions
     * @param maxIndels maximal number of insertion and deletion sum
     * @param maxLevensteinDistance maximal levenstein distance
     * @param scoring alignment scoring scheme
     * @param exhaustive if set to false stop with first variant from tree search, otherwise search all possible re-alignments
     * @param greedy if set to true stop if new variant has more mismatches than the previous one
     */
    SearchScope(int maxSubstitutions, int maxIndels, int maxLevensteinDistance,
                boolean exhaustive = false,
                boolean greedy = true) {
        this.parameters = new TreeSearchParameters(maxSubstitutions, maxIndels, maxIndels, maxLevensteinDistance)
        this.maxIndels = maxIndels
        this.exhaustive = exhaustive
        this.greedy = greedy
    }

    /**
     * Create sequence search filter parameter set. Specifies maximal allowed number of substitutions, insertions,
     * deletions and Levenstein distance from query sequence. Note that resulting search is not 'symmetric' as
     * number of insertions and deletions may not match. Used in 'Browse' tab of VDJdb web app.
     *
     * @param maxSubstitutions maximal number of substitutions
     * @param maxInsertions maximal number of insertions
     * @param maxDeletions maximal number of deletions
     * @param maxLevensteinDistance maximal levenstein distance
     * @param scoring alignment scoring scheme
     * @param exhaustive if set to false stop with first variant from tree search, otherwise search all possible re-alignments
     * @param greedy if set to true stop if new variant has more mismatches than the previous one
     */
    SearchScope(int maxSubstitutions, int maxInsertions, int maxDeletions, int maxLevensteinDistance,
                boolean exhaustive = false,
                boolean greedy = true) {
        this.parameters = new TreeSearchParameters(maxSubstitutions, maxInsertions, maxDeletions, maxLevensteinDistance)
        this.maxIndels = parameters.maxInsertions + parameters.maxDeletions
        this.exhaustive = exhaustive
        this.greedy = greedy
    }
}
