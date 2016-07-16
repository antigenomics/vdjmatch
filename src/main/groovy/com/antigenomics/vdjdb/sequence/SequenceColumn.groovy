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

        // This hack excludes duplicates like
        //
        //   KLFF     KLFF     ...
        //   KLF-     KL-F
        def prevSeq = new HashSet<String>()

        while (true) {
            def entries = ni.next()

            if (entries) {
                def seq = entries.first().value
                if (!prevSeq.contains(seq)) {
                    def mutations = ni.currentMutations
                    float score = mutations.size() == 0 ? scoring.scoreThreshold : // always include exact match
                            scoring.computeScore(mutations, baseScore, refLength)

                    if (score >= scoring.scoreThreshold) {
                        def alignment = new Alignment(filter.query, mutations,
                                new Range(0, refLength), new Range(0, seq.size()),
                                score)
                        entries.each { entry ->
                            results.put(entry.row, alignment)
                        }
                    }
                    prevSeq.add(seq)
                }
            } else {
                break
            }
        }

        results
    }
}
