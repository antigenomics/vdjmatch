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

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.misc.ExecUtil
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentHashMap
import java.util.function.Function

class ClonotypeSearchSummary {
    static final public List<String> FIELDS_STARBURST = ["mhc.class",
                                                  "mhc.a",
                                                  "mhc.b",
                                                  "antigen.species",
                                                  "antigen.gene",
                                                  "antigen.epitope"],
                              FIELDS_PLAIN_TEXT = ["mhc.class",
                                                   "antigen.species",
                                                   "antigen.gene",
                                                   "antigen.epitope"]

    final public Map<String, Map<String, ClonotypeCounter>> fieldCounters = new ConcurrentHashMap<>()

    final ClonotypeCounter totalCounter = new ClonotypeCounter(),
                           notFoundCounter

    ClonotypeSearchSummary(Map<Clonotype, List<ClonotypeSearchResult>> searchResults,
                           Sample sample,
                           List<String> columnNameList) {
        columnNameList.each {
            fieldCounters.put(it, new ConcurrentHashMap<String, ClonotypeCounter>())
        }

        GParsPool.withPool ExecUtil.THREADS, {
            searchResults.eachParallel { Map.Entry<Clonotype, List<ClonotypeSearchResult>> clonotypeResult ->
                clonotypeResult.value.each { result ->
                    columnNameList.each { columnId ->
                        def subMap = fieldCounters[columnId],
                            value = result.row[columnId].value

                        def counter = subMap.computeIfAbsent(value, countergen)
                        counter.update(clonotypeResult.key)
                    }
                }
                totalCounter.update(clonotypeResult.key)
            }
        }

        notFoundCounter = new ClonotypeCounter(sample.diversity - totalCounter.unique,
                sample.count - totalCounter.reads,
                sample.freq - totalCounter.frequency)
    }

    ClonotypeCounter getCounter(String columnId, String value) {
        if (!fieldCounters.containsKey(columnId)) {
            throw new RuntimeException("Column $columnId was not summarized.")
        }

        def counter = fieldCounters[columnId][value]

        counter ?: new ClonotypeCounter()
    }

    private static final Function<String, ClonotypeCounter> countergen = new Function<String, ClonotypeCounter>() {
        @Override
        ClonotypeCounter apply(String s) {
            new ClonotypeCounter()
        }
    }
}
