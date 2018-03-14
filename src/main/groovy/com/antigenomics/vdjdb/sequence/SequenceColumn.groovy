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
import com.fasterxml.jackson.annotation.JsonIgnore
import com.milaboratory.core.mutations.Mutations
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
    Map<Row, Hit> search(SequenceFilter filter) {
        def ni = stm.getNeighborhoodIterator(filter.query,
                filter.treeSearchParameters) // 'search scope' iterator

        def scoring = filter.alignmentScoring // alignment and global scoring

        def matchBuffer = new HashMap<String, SequenceSearchResult>()
        List<Entry> entries

        while ((entries = ni.next()) != null) { // until no more alignments found within 'search scope'
            def matchSequence = entries.first().value
            def mutations = ni.currentMutations

            // need this workaround as it is not possible to implement (ins+dels) <= X scope with tree searcher
            // due to separate insertion and deletion counting
            if (mutations.countOfIndels() > filter.maxIndels)
                continue

            def previousResult = matchBuffer[matchSequence] // need buffer here as hits are not guaranteed to be ordered
            // by match sequence. using hashmap as we'll need to store
            // previous score

            if (previousResult == null) {
                // compute new if we don't have any previous alignments with this sequence
                float alignmentScore = scoring.computeScore(filter.query, mutations)
                matchBuffer.put(matchSequence, new SequenceSearchResult(alignmentScore, mutations, entries))
            } else if (filter.exhaustive) { // exhaustive search - compare scores
                float alignmentScore = scoring.computeScore(filter.query, mutations)

                if (alignmentScore > previousResult.alignmentScore) {
                    // replace if better score
                    matchBuffer.put(matchSequence, new SequenceSearchResult(alignmentScore, mutations, entries))
                }
            }
        }

        // Flatten results and wrap them into 'hits'

        def results = new HashMap<Row, Hit>()

        matchBuffer.each { matchKvp ->
            def searchResult = matchKvp.value
            searchResult.entries.each { entry -> // iterate through matched rows
                results.put(entry.row, new Hit(filter.query, // query sequence
                        searchResult.mutations, // query -> db match mutations
                        matchKvp.key.length(), // db match sequence length
                        searchResult.alignmentScore // alignment (e.g. CDR3 alignment) score
                ))
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

    private class SequenceSearchResult {
        final float alignmentScore
        final Mutations<AminoAcidSequence> mutations
        final List<Entry> entries

        SequenceSearchResult(float alignmentScore, Mutations<AminoAcidSequence> mutations, List<Entry> entries) {
            this.alignmentScore = alignmentScore
            this.mutations = mutations
            this.entries = entries
        }
    }
}
