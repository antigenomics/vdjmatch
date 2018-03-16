package com.antigenomics.vdjdb.web

import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.milaboratory.core.mutations.MutationType
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.TreeSearchParameters

/**
 * Created by mikesh on 7/23/17.
 */
class EpitopeSuggestionGenerator {
    static final String EPITOPE_COLUMN_NAME = "antigen.epitope"

    static Map<String, List<EpitopeSuggestion>> generateSuggestions(VdjdbInstance vdjdbInstance,
                                                                    int maxMm = 3, int maxIndel = 2,
                                                                    int maxTotal = 5) {
        def epiColumn = vdjdbInstance.dbInstance[EPITOPE_COLUMN_NAME] as SequenceColumn

        epiColumn.values.collectEntries() { value ->
            def filter = new SequenceFilter(EPITOPE_COLUMN_NAME, new AminoAcidSequence(value),
                    new TreeSearchParameters(maxMm, maxIndel, maxIndel, maxTotal, false))

            def suggestionSet = new HashSet<EpitopeSuggestion>()
            epiColumn.search(filter).each { result ->
                def sequence = result.key[EPITOPE_COLUMN_NAME].value

                if (sequence != value) {
                    def muts = result.value.mutations

                    suggestionSet.add(new EpitopeSuggestion(sequence,
                            muts.countOf(MutationType.Substitution),
                            muts.countOfIndels(),
                            sequence.length(),
                            epiColumn.stm.get(new AminoAcidSequence(sequence)).size()))
                }
            }

            [(value): suggestionSet.sort { -it.substitutions - it.indels }]
        }
    }
}
