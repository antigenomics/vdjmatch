package com.antigenomics.vdjdb.core.query

import com.antigenomics.vdjdb.core.db.CdrDatabase
import com.antigenomics.vdjdb.core.db.CdrEntrySet
import com.milaboratory.core.Range
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.mutations.Mutations
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters

class CdrDatabaseSearcher {
    private final SequenceTreeMap<AminoAcidSequence, CdrEntrySet> stm
    private final TreeSearchParameters params
    private final int depth
    private final CdrDatabase database

    public CdrDatabaseSearcher(CdrDatabase database) {
        this(database, 2, 1, 1, 3, 10)
    }

    public CdrDatabaseSearcher(CdrDatabase database,
                               int maxSubstitutions, int maxInsertions, int maxDeletions,
                               int maxMutations, int depth) {
        this.stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET)
        this.database = database
        database.each {
            stm.put(new AminoAcidSequence(it.cdr3aa), it)
        }
        this.params = new TreeSearchParameters(maxSubstitutions, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    public List<CdrSearchResult> search(String cdr3aa) {
        search(new AminoAcidSequence(cdr3aa))
    }

    public List<CdrSearchResult> search(AminoAcidSequence cdr3aa) {
        def results = new ArrayList<CdrSearchResult>(depth)
        def ni = stm.getNeighborhoodIterator(cdr3aa, params)

        // This hack excludes duplicates like
        //
        //   KLFF     KLFF     ...
        //   KLF-     KL-F
        def prevCdr = new HashSet<String>()

        int i = 0
        while (true) {
            def entry = ni.next()

            if (entry) {
                def alignment = ni.getCurrentAlignment()
                def match = entry.cdr3aa
                if (!prevCdr.contains(match)) {
                    results.add(new CdrSearchResult(cdr3aa, alignment, entry))
                    prevCdr.add(match)
                }
            } else {
                break
            }

            if (depth > 0 && ++i == depth)
                break
        }

        results.sort()
    }

    public CdrSearchResult lucky(String cdr3aa) {
        def result = search(cdr3aa)
        result.size() > 0 ? result[0] : null
    }

    public CdrSearchResult exact(String cdr3aa) {
        def match = database[cdr3aa]
        if (!match)
            return null

        def query = new AminoAcidSequence(cdr3aa)

        def dummyRange = new Range(0, cdr3aa.length())
        def dummyAlignment = new Alignment(query,
                new Mutations(AminoAcidSequence.ALPHABET),
                dummyRange, dummyRange,
                cdr3aa.length())

        new CdrSearchResult(query, dummyAlignment, match)
    }
}
