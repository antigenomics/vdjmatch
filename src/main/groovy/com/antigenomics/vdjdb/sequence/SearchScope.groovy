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

import com.milaboratory.core.tree.TreeSearchParameters

/**
 * Stores parameters for SequenceFilter specifying 'search scope' (number of allowed mismatches) and alignment
 * scoring method.
 */
class SearchScope {
    /**
     * Maximum number of allowed substitutions
     */
    final int maxSubstitutions
    /**
     * Maximum number of allowed deletions
     */
    final int maxDeletions
    /**
     * Maximum number of allowed insertions
     */
    final int maxInsertions
    /**
     * Maximum number of allowed insertions + deletions
     */
    final int maxIndels
    /**
     * Maximum number of allowed mutations (substitutions + indels, i.e. edit distance)
     */
    final int maxTotal
    /**
     * If set to true will search all possible re-alignments between two pairs of two sequences
     * and select the one with highest score. If set to false will take the first alignment
     */
    final boolean exhaustive
    /**
     * Do not consider alignment variants with more number of mismatches when doing exhaustive search within the scope
     */
    final boolean greedy

    final static SearchScope EXACT = new SearchScope(0, 0, 0, 0,
            false, false)

    /**
     * Create sequence search filter parameter set. Specifies maximal allowed number of substitutions, indels
     * (i.e. sum of insertions and deletions, or difference in query and target sequence lengths) and
     * edit distance from query sequence. Note that resulting search is 'symmetric' as
     * number of insertions in query is the number of deletions in target. Used in 'Annotate' tab of VDJdb web app.
     *
     * @param maxSubstitutions maximal number of substitutions
     * @param maxIndels maximal number of insertions/deletions
     * @param maxTotal maximal levenstein distance
     * @param exhaustive if set to false stop with first alignment from tree search, otherwise search all possible re-alignments
     * @param greedy if set to false, will also consider alignments with more mismatches than previous hits (only applicable within search scope and when exhaustive search is on)
     */
    SearchScope(int maxSubstitutions, int maxIndels, int maxTotal,
                boolean exhaustive = true, boolean greedy = true) {
        this.maxSubstitutions = maxSubstitutions
        this.maxTotal = maxTotal
        this.maxDeletions = maxIndels
        this.maxInsertions = maxIndels
        this.maxIndels = maxIndels
        this.exhaustive = exhaustive
        this.greedy = greedy
    }

    /**
     * Create sequence search filter parameter set. Specifies maximal allowed number of substitutions, insertions,
     * deletions and edit distance from query sequence. Note that resulting search is not 'symmetric' as
     * number of insertions and deletions may not match. Used in 'Browse' tab of VDJdb web app.
     *
     * @param maxSubstitutions maximal number of substitutions
     * @param maxDeletions maximal number of deletions
     * @param maxInsertions maximal number of insertions
     * @param maxTotal maximal levenstein distance
     * @param exhaustive if set to false stop with first alignment from tree search, otherwise search all possible re-alignments
     * @param greedy if set to false, will also consider alignments with more mismatches than previous hits (only applicable within search scope and when exhaustive search is on)
     */
    SearchScope(int maxSubstitutions, int maxDeletions, int maxInsertions, int maxTotal,
                boolean exhaustive = true, boolean greedy = true) {
        this.maxSubstitutions = maxSubstitutions
        this.maxTotal = maxTotal
        this.maxDeletions = maxDeletions
        this.maxInsertions = maxInsertions
        this.maxIndels = maxInsertions + maxDeletions
        this.exhaustive = exhaustive
        this.greedy = greedy
    }

    /**
     * Generate parameters for sequence tree map search
     * @return sequence tree map search parameters
     */
    TreeSearchParameters getTreeSearchParameters() {
        new TreeSearchParameters(
                maxSubstitutions,
                maxDeletions,
                maxInsertions,
                maxTotal,
                !exhaustive)
    }
}
