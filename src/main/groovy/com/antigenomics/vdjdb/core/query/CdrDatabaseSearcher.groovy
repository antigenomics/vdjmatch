package com.antigenomics.vdjdb.core.query

import com.antigenomics.vdjdb.core.db.CdrDatabase
import com.antigenomics.vdjdb.core.db.CdrEntrySet
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters

class CdrDatabaseSearcher {
    private final SequenceTreeMap<AminoAcidSequence, CdrEntrySet> stm
    private final TreeSearchParameters params
    private final int depth

    public CdrDatabaseSearcher(CdrDatabase database) {
        this(database, 2, 1, 1, 3, 10)
    }

    public CdrDatabaseSearcher(CdrDatabase database,
                               int maxSubstitutions, int maxInsertions, int maxDeletions,
                               int maxMismatches, int depth) {
        this.stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET)
        database.each {
            stm.put(new AminoAcidSequence(it.cdr3aa), it)
        }
        this.params = new TreeSearchParameters(maxSubstitutions, maxInsertions, maxDeletions, maxMismatches)
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

        for (int i = 0; i < depth; i++) {
            def entry = ni.next()
            if (entry) {
                def alignment = ni.getCurrentAlignment()
                def match = entry.cdr3aa
                if (!prevCdr.contains(match)) {
                    results.add(new CdrSearchResult(alignment, entry))
                    prevCdr.add(match)
                }
            } else {
                break
            }
        }
        results.sort()
    }

    public CdrSearchResult lucky(String cdr3aa) {
        search(cdr3aa)[0]
    }
}
