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
    static final List<String> FIELDS_AG = ["antigen.species", "antigen.gene"],
                              FIELDS_MHC_A = ["mhc.class", "mhc.a"],
                              FIELDS_MHC_B = ["mhc.class", "mhc.a", "mhc.b"],
                              FIELDS_PLAIN_TEXT = ["mhc.class", "mhc.a", "mhc.b",
                                                   "antigen.species", "antigen.gene"]

    final Map<String, Map<String, ClonotypeCounter>> columnSequenceCounters = new ConcurrentHashMap<>()

    final ClonotypeCounter totalCounter = new ClonotypeCounter(),
                           notFoundCounter

    private static final Function<String, Map<String, ClonotypeCounter>> mapgen =
            new Function<String, Map<String, ClonotypeCounter>>() {
                @Override
                Map<String, ClonotypeCounter> apply(String s) {
                    new ConcurrentHashMap<>()
                }
            }

    private static final Function<String, ClonotypeCounter> countergen = new Function<String, ClonotypeCounter>() {
        @Override
        ClonotypeCounter apply(String s) {
            new ClonotypeCounter()
        }
    }

    static ClonotypeSearchSummary createForStarburst(Map<Clonotype, List<ClonotypeSearchResult>> searchResults,
                                                     Sample sample) {
        def fields = listAllNameSequences(FIELDS_AG)
        fields.addAll(listAllNameSequences(FIELDS_MHC_A))
        fields.addAll(listAllNameSequences(FIELDS_MHC_B))
        new ClonotypeSearchSummary(searchResults, sample, fields)
    }

    static List<List<String>> listAllNameSequences(List<String> columnNames) {
        (0..<columnNames.size()).collect { columnNames[0..it] }
    }

    ClonotypeSearchSummary(Map<Clonotype, List<ClonotypeSearchResult>> searchResults,
                           Sample sample,
                           List<List<String>> columnNameListCombinations) {
        def counterTypes = columnNameListCombinations.collectEntries { [(it.join("\t")): it] }

        GParsPool.withPool ExecUtil.THREADS, {
            searchResults.eachParallel { Map.Entry<Clonotype, List<ClonotypeSearchResult>> clonotypeResult ->
                clonotypeResult.value.each { result ->
                    // todo can optimize by moving this upward, needs db though
                    counterTypes.each { Map.Entry<String, List<String>> counterType ->
                        def subMap = columnSequenceCounters.computeIfAbsent(counterType.key, mapgen)
                        // todo can optimize by pre-getting column indices
                        def combination = counterType.value.collect { result.row[it] }.join("\t")
                        def counter = subMap.computeIfAbsent(combination, countergen)
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

    ClonotypeCounter getCounter(List<String> columnNameSequence, List<String> columnValueSequence) {
        getCounters(columnNameSequence)[columnValueSequence.join("\t")]
    }

    Map<String, ClonotypeCounter> getCounters(List<String> columnNameSequence) {
        columnSequenceCounters[columnNameSequence.join("\t")]
    }

    Map<String, ClonotypeCounter> getCounters() {
        columnSequenceCounters.values().first()
    }
}
