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

package com.antigenomics.vdjdb.stat

import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.impl.ClonotypeSearchResult
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import groovyx.gpars.GParsPool

/**
 * A database search summary implementation for clonotype sample 
 */
class ClonotypeSearchSearchSummary extends SearchSummary {
    /**
     * Query sample 
     */
    final Sample sample

    /**
     * Creates an empty clonotype sample search summary
     * @param database a clonotype database
     * @param columnNames column names to generate summary for
     * @param sample query sample
     */
    ClonotypeSearchSearchSummary(ClonotypeDatabase database, List<String> columnNames, Sample sample) {
        super(database, columnNames)
        this.sample = sample
    }

    /**
     * Appends clonotype database search results to the summary (in parallel) 
     * @param searchResult clonotype database search results
     */
    void append(Map<Clonotype, List<ClonotypeSearchResult>> searchResult) {
        GParsPool.withPool ExecUtil.THREADS, {
            searchResult.entrySet().eachParallel { Map.Entry<Clonotype, List<ClonotypeSearchResult>> entry ->
                append(entry.value*.row, entry.key.freq)
            }
        }
    }

    /**
     * Gets the number and frequency of clonotypes that were not found
     * @return missing clonotype counter
     */
    Counter getNotFound() {
        new CounterImpl(sample.diversity - foundCounter.uniqueCount,
                1.0 - foundCounter.weightedCount)
    }
}
