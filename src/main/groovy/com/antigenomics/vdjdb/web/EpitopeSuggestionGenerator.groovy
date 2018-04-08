package com.antigenomics.vdjdb.web

import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.milaboratory.core.mutations.MutationType
import com.milaboratory.core.sequence.AminoAcidSequence

/**
 * Created by mikesh on 7/23/17.
 */
class EpitopeSuggestionGenerator {
    static Map<String, List<EpitopeSuggestion>> generateSuggestions(VdjdbInstance vdjdbInstance,
                                                                    int maxMm = 3, int maxIndel = 2,
                                                                    int maxTotal = 5) {
        def epiColumn = vdjdbInstance.dbInstance[ClonotypeDatabase.EPITOPE_COL_DEFAULT] as SequenceColumn

        epiColumn.values.collectEntries() { value ->
            def filter = new SequenceFilter(ClonotypeDatabase.EPITOPE_COL_DEFAULT, new AminoAcidSequence(value),
                    new SearchScope(maxMm, maxIndel, maxTotal, false))

            def suggestionSet = new HashSet<EpitopeSuggestion>()
            epiColumn.search(filter).each { result ->
                def sequence = result.key[ClonotypeDatabase.EPITOPE_COL_DEFAULT].value

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
