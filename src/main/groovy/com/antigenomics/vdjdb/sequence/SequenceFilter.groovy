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
     * Sequence search parameters/thresholds
     */
    final SearchScope searchScope

    /**
     * Sequence alignment scoring
     */
    final AlignmentScoring alignmentScoring

    /**
     * Creates a new amino acid sequence search rule
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param searchScope sequence search parameters/thresholds
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     */
    SequenceFilter(String columnId, String query,
                   SearchScope searchScope = SearchScope.EXACT,
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE) {
        this(columnId, Util.convert(query), searchScope, alignmentScoring)
    }

    /**
     * Creates a new amino acid sequence search rule 
     * @param columnId sequence column id
     * @param query amino acid sequence query
     * @param searchScope sequence search parameters/thresholds
     * @param alignmentScoring alignment scoring containing substitution and gap scores, as well as a total score threshold
     */
    SequenceFilter(String columnId, AminoAcidSequence query,
                   SearchScope searchScope = SearchScope.EXACT,
                   AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE) {
        this.columnId = columnId
        this.query = query
        this.searchScope = searchScope
        this.alignmentScoring = alignmentScoring

        if (columnId == null || query == null)
            throw new RuntimeException("Bad sequence filter query")
    }

    @Override
    boolean isSequenceFilter() {
        true
    }
}
