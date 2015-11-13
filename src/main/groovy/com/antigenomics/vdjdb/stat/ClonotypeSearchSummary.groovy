/*
 * Copyright 2013-{year} Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjdb.stat

import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.impl.ClonotypeSearchResult
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import groovyx.gpars.GParsPool

class ClonotypeSearchSummary extends SummaryStatistics {
    final Sample sample

    ClonotypeSearchSummary(Database database, List<String> columnNames, Sample sample) {
        super(database, columnNames)
        this.sample = sample
    }

    void append(Map<Clonotype, List<ClonotypeSearchResult>> searchResult) {
        GParsPool.withPool ExecUtil.THREADS, {
            searchResult.eachParallel {
                append(it.value*.row, it.key.freq)
            }
        }
    }

    Counter getNotFoundCounter() {
        new CounterImpl(sample.diversity - foundCounter.uniqueCount,
                1.0 - foundCounter.weightedCount)
    }
}
