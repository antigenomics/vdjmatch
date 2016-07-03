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
     * Search depth
     */
    final int depth
    /**
     * Alignment scoring
     */
    final AlignmentScoring alignmentScoring

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query query to be converted to amino acid sequence
     * @param treeSearchParameters alignment parameters
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     * @param depth search depth
     */
    SequenceFilter(String columnId, String query,
                   TreeSearchParameters treeSearchParameters = new TreeSearchParameters(5, 2, 2, 7),
                   int depth = -1,
                   AlignmentScoring alignmentScoring = AlignmentScoringProvider.loadScoring()) {
        this(columnId, Util.convert(query), treeSearchParameters, depth, alignmentScoring)
    }

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param treeSearchParameters alignment parameters
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     * @param depth search depth
     */
    SequenceFilter(String columnId, AminoAcidSequence query,
                   TreeSearchParameters treeSearchParameters = new TreeSearchParameters(5, 2, 2, 7),
                   int depth = -1,
                   AlignmentScoring alignmentScoring = AlignmentScoringProvider.loadScoring()) {
        if (query == null)
            throw new RuntimeException("Bad sequence filter query")

        this.columnId = columnId
        this.query = query
        this.treeSearchParameters = treeSearchParameters
        this.depth = depth
        this.alignmentScoring = alignmentScoring
    }

    @Override
    boolean isSequenceFilter() {
        true
    }
}
