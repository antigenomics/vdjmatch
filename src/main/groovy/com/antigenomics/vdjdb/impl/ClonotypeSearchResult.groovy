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

package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.db.Row
import com.antigenomics.vdjdb.db.SearchResult
import com.antigenomics.vdjdb.sequence.Hit
import groovy.transform.CompileStatic

/**
 * Clonotype search result 
 */
@CompileStatic
class ClonotypeSearchResult implements Comparable<ClonotypeSearchResult>, SearchResult {
    /**
     * CDR3 sequence alignment result 
     */
    final Hit hit

    /**
     * Database row that was found
     */
    final Row row

    /**
     * Clonotype id (order in sample), if specified, -1 otherwise
     */
    final int id

    /**
     * TCR similarity score
     */
    final float score

    /**
     * Weight/informativeness of a hit
     */
    final float weight

    /**
     * A product of score and weight
     */
    final float weightedScore

    /**
     * Creates a new clonotype search result
     * @param result CDR3 sequence alignment result
     * @param row database row
     * @param id clonotype id in sample
     * @param score full TCR similarity score
     * @param weight hit weight/informativeness
     */
    ClonotypeSearchResult(Hit hit, Row row, int id, float score, float weight) {
        this.hit = hit
        this.row = row
        this.id = id
        this.score = score
        this.weight = weight
        this.weightedScore = weight * score
    }

    @Override
    int compareTo(ClonotypeSearchResult o) {
        -Double.compare(weightedScore, o.weightedScore)
    }
}
