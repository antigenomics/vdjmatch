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

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.db.Filter
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.TreeSearchParameters

/**
 * An amino acid sequence search rule. Holds sequence search parameters (max number of mismatches and alignment
 * scoring method) and query sequence
 */
class SequenceFilter implements Filter {
    /**
     * Sequence column id 
     */
    final String columnId

    /**
     * Query sequence 
     */
    final AminoAcidSequence query

    /**
     * Sequence tree search parameters/thresholds
     */
    final TreeSearchParameters treeSearchParameters

    /**
     * Sequence alignment scoring
     */
    final AlignmentScoring alignmentScoring

    /**
     * Maximum number of allowed insertions + deletions
     */
    final int maxIndels

    /**
     * If set to true will search all possible re-alignments between two pairs of two sequences
     * and select the one with highest score. If set to false will take the first alignment
     */
    final boolean exhaustive

    /**
     * Do not consider variants with more number of mismatches when doing exhaustive search within scope
     */
    final boolean greedy

    /**
     * Creates a new amino acid sequence search rule
     * @param columnId sequence column id
     * @param query query sequence, will be converted to amino acid sequence
     * @param preset search scope parameters
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     */
    SequenceFilter(String columnId, String query,
                   SearchScope preset,
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE) {
        this(columnId, query, preset.parameters, preset.maxIndels, alignmentScoring, preset.exhaustive)
    }

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query query sequence, will be converted to amino acid sequence
     * @param treeSearchParameters alignment parameters
     * @param maxIndels maximum allowed sequence length difference
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     * @param exhaustive if set to false stop with first variant from tree search, otherwise search all possible re-alignments
     * @param greedy if set to false, will also consider variants with more mismatches than previous hits (only applicable within search scope and when exhaustive search is on)
     */
    SequenceFilter(String columnId, String query,
                   TreeSearchParameters treeSearchParameters = new TreeSearchParameters(0, 0, 0, 0, false),
                   int maxIndels = -1,
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                   boolean exhaustive = true,
                   boolean greedy = true) {
        this(columnId, Util.convert(query), treeSearchParameters, maxIndels, alignmentScoring, exhaustive, greedy)
    }

    /**
     * Creates a new amino acid sequence search rule
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param preset search scope parameters
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     */
    SequenceFilter(String columnId, AminoAcidSequence query,
                   SearchScope preset,
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE) {
        this(columnId, query, preset.parameters, preset.maxIndels, alignmentScoring, preset.exhaustive, preset.greedy)
    }

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param treeSearchParameters alignment parameters
     * @param maxIndels maximum allowed sequence length difference
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     * @param exhaustive if set to false stop with first variant from tree search, otherwise search all possible re-alignments
     * @param greedy if set to false, will also consider variants with more mismatches than previous hits (only applicable within search scope and when exhaustive search is on)
     */
    SequenceFilter(String columnId, AminoAcidSequence query,
                   TreeSearchParameters treeSearchParameters = new TreeSearchParameters(0, 0, 0, 0),
                   int maxIndels = -1,
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                   boolean exhaustive = true,
                   boolean greedy = true) {
        this.columnId = columnId
        this.query = query
        this.treeSearchParameters = treeSearchParameters
        this.maxIndels = maxIndels == -1 ? (treeSearchParameters.maxInsertions + treeSearchParameters.maxDeletions) : maxIndels
        this.alignmentScoring = alignmentScoring
        this.exhaustive = exhaustive
        this.greedy = greedy

        if (treeSearchParameters.greedy != !exhaustive)
            throw new RuntimeException("Conflicting tree search parameters (greedy) and exhaustive mode")

        if (query == null)
            throw new RuntimeException("Bad sequence filter query")
    }

    @Override
    boolean isSequenceFilter() {
        true
    }
}
