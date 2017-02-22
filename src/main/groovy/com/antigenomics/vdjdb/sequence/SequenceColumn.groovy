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
import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.Entry
import com.antigenomics.vdjdb.db.Row
import com.antigenomics.vdjdb.db.SearchResult
import com.fasterxml.jackson.annotation.JsonIgnore
import com.milaboratory.core.Range
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import groovy.transform.CompileStatic

/**
 * A column containing amino acid sequences. The column is indexed automatically when new entries are added 
 */
@CompileStatic
class SequenceColumn extends Column {
    final SequenceTreeMap<AminoAcidSequence, List<Entry>> stm = new SequenceTreeMap(AminoAcidSequence.ALPHABET)

    /**
     * {@inheritDoc}
     */
    SequenceColumn(String name, Map<String, String> metadata = [:]) {
        super(name, metadata)
    }

    /**
     * Adds a new entry to the database. Entry is skipped if it is not an amino acid sequence. 
     * @param entry amino acid sequence entry
     */
    @Override
    void append(Entry entry) {
        if (entry.value.length() > 0) {
            def seq = Util.convert(entry.value)
            if (seq) {
                def entries = stm.get(seq)
                if (entries == null) {
                    stm.put(seq, entries = new ArrayList<Entry>())
                }
                entries.add(entry)
            }
        }
    }

    /**
     * Searches the column for an amino acid sequence
     * @param filter amino acid query
     * @return a map of rows that were found and corresponding sequence alignment results
     */
    Map<Row, Alignment> search(SequenceFilter filter) {
        def results = new HashMap<Row, Alignment>()

        def ni = stm.getNeighborhoodIterator(filter.query,
                filter.treeSearchParameters)

        def scoring = filter.alignmentScoring

        def baseScore = scoring.computeBaseScore(filter.query),
            refLength = filter.query.size()


        def resultBuffer = new HashMap<String, SearchResult>()
        List<Entry> entries

        while ((entries = ni.next()) != null) {
            def seq = entries.first().value
            def mutations = ni.currentMutations
            def previousResult = resultBuffer[seq]

            if (previousResult == null) {
                // compute new if we don't have any previous alignments with this sequence
                float score = scoring.computeScore(mutations, baseScore, refLength)
                def alignment = new Alignment(filter.query, mutations,
                        new Range(0, refLength), new Range(0, seq.size()),
                        score)

                resultBuffer.put(seq,
                        new SearchResult(alignment, entries))
            } else if (filter.exhaustive) {
                def previousMutationCount = previousResult.alignment.absoluteMutations.size()
                // exhaustive search - compare scores
                if (!filter.greedy || // non-greedy case: check even if we have more mutations
                        previousMutationCount >= mutations.size()) { // greedy: previous case has the same number of mutations or more
                    // Note: more mutations case N/A to current tree impl
                    float score = scoring.computeScore(mutations, baseScore, refLength)

                    if (score > previousResult.alignment.score) {
                        // replace if better score
                        def alignment = new Alignment(filter.query, mutations,
                                new Range(0, refLength), new Range(0, seq.size()),
                                score)

                        resultBuffer.put(seq,
                                new SearchResult(alignment, entries))
                    }
                }
            }
        }

        resultBuffer.values().each { result ->
            result.entries.each { entry ->
                results.put(entry.row, result.alignment)
            }
        }

        results
    }

    /**
     * Gets the set of all possible values in the column
     * @return a set of unique values in the column
     */
    @JsonIgnore
    @Override
    Set<String> getValues() {
        def values = new HashSet<String>()

        stm.values().each { List<Entry> entryList ->
            entryList.each { Entry entry ->
                values.add(entry.value)
            }
        }

        values
    }

    private class SearchResult {
        Alignment alignment
        List<Entry> entries

        SearchResult(Alignment alignment, List<Entry> entries) {
            this.alignment = alignment
            this.entries = entries
        }
    }
}
