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
import com.antigenomics.vdjdb.scoring.AlignmentScoring
import com.antigenomics.vdjdb.scoring.DummyAlignmentScoring
import com.antigenomics.vdjdb.scoring.SequenceSearcherPreset
import com.antigenomics.vdjdb.scoring.VdjdbAlignmentScoring
import com.antigenomics.vdjdb.scoring.AlignmentScoringProvider
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.TreeSearchParameters

/**
 * An amino acid sequence search rule 
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
     * Alignment parameters 
     */
    final TreeSearchParameters treeSearchParameters

    /**
     * Alignment scoring
     */
    final AlignmentScoring alignmentScoring

    /**
     * If set to true will search all possible re-alignments between two pairs of two sequences
     * and select the one with highest score. If set to false will take the first alignment
     */
    final boolean exhaustive

    /**
     * If set to true stop if new variant has more mismatches than the previous one when doing exhaustive search, i.e.
     * will only consider and score cases with equal number of substitutions + indels.
     */
    final boolean greedy

    /**
     * Creates a new amino acid sequence search rule
     * @param columnId sequence column id
     * @param query query sequence, will be converted to amino acid sequence
     * @param preset alignment and scoring parameter preset
     */
    SequenceFilter(String columnId, String query,
                   SequenceSearcherPreset preset) {
        this(columnId, query, preset.parameters, preset.scoring, preset.exhaustive, preset.greedy)
    }

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query query sequence, will be converted to amino acid sequence
     * @param treeSearchParameters alignment parameters
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     * @param exhaustive if set to false stop with first variant from tree search, otherwise search all possible re-alignments
     * @param greedy if set to true stop if new variant has more mismatches than the previous one
     */
    SequenceFilter(String columnId, String query,
                   TreeSearchParameters treeSearchParameters = new TreeSearchParameters(0, 0, 0, 0),
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                   boolean exhaustive = false,
                   boolean greedy = true) {
        this(columnId, Util.convert(query), treeSearchParameters, alignmentScoring, exhaustive, greedy)
    }

    /**
     * Creates a new amino acid sequence search rule
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param preset alignment and scoring parameter preset
     */
    SequenceFilter(String columnId, AminoAcidSequence query,
                   SequenceSearcherPreset preset) {
        this(columnId, query, preset.parameters, preset.scoring, preset.exhaustive, preset.greedy)
    }

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param treeSearchParameters alignment parameters
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     * @param exhaustive if set to false stop with first variant from tree search, otherwise search all possible re-alignments
     * @param greedy if set to true stop if new variant has more mismatches than the previous one
     */
    SequenceFilter(String columnId, AminoAcidSequence query,
                   TreeSearchParameters treeSearchParameters = new TreeSearchParameters(0, 0, 0, 0),
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                   boolean exhaustive = false,
                   boolean greedy = true) {
        this.columnId = columnId
        this.query = query
        this.treeSearchParameters = treeSearchParameters
        this.alignmentScoring = alignmentScoring
        this.exhaustive = exhaustive
        this.greedy = greedy

        if (query == null)
            throw new RuntimeException("Bad sequence filter query")
    }

    @Override
    boolean isSequenceFilter() {
        true
    }
}
