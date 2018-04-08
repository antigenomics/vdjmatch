package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.HashSetStringGen
import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.db.Entry
import com.antigenomics.vdjdb.impl.ClonotypeDatabase

class EpitopeSizeFilterUtil {
    static Map<String, Integer> generateEpitopeCounts(VdjdbInstance vdjdbInstance,
                                                      String species = null,
                                                      String gene = null) {
        def epiCdr3Map = new HashMap<String, Set<String>>()

        def filters = []

        if (species) {
            filters << new ExactTextFilter(ClonotypeDatabase.SPECIES_COL_DEFAULT, species, false)
        }
        if (gene) {
            filters << new ExactTextFilter(ClonotypeDatabase.GENE_COL_DEFAULT, gene, false)
        }

        if (filters) {
            vdjdbInstance = vdjdbInstance.filter(filters)
        }

        vdjdbInstance.dbInstance.rows.each { row ->
            def epi = row[ClonotypeDatabase.EPITOPE_COL_DEFAULT].value,
                cdr3 = row[ClonotypeDatabase.CDR3_COL_DEFAULT].value

            epiCdr3Map.computeIfAbsent(epi, HashSetStringGen.INSTANCE).add(cdr3)
        }

        (epiCdr3Map.collectEntries { [(it.key): it.value.size()] } as Map<String, Integer>)
    }

    static TextFilter createEpitopeSizeFilter(VdjdbInstance vdjdbInstance,
                                              String species,
                                              String gene,
                                              int countThreshold) {
        def goodEpitopes = new HashSet<String>()

        generateEpitopeCounts(vdjdbInstance, species, gene).each {
            if (it.value >= countThreshold) {
                goodEpitopes.add(it.key)
            }
        }

        new Filter(goodEpitopes)
    }

    private static class Filter extends TextFilter {
        final Set<String> values

        Filter(Set<String> values) {
            super(ClonotypeDatabase.EPITOPE_COL_DEFAULT, null, false)
            this.values = values
        }

        @Override
        protected boolean passInner(Entry entry) {
            values.contains(entry.value)
        }
    }
}
