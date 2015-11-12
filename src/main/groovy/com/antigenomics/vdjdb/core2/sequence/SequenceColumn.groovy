package com.antigenomics.vdjdb.core2.sequence

import com.antigenomics.vdjdb.core2.db.Column
import com.antigenomics.vdjdb.core2.db.ColumnType
import com.antigenomics.vdjdb.core2.db.Entry
import com.antigenomics.vdjdb.core2.db.Row
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import groovy.transform.CompileStatic

@CompileStatic
class SequenceColumn extends Column {
    final SequenceTreeMap<AminoAcidSequence, List<Entry>> stm = new SequenceTreeMap(AminoAcidSequence.ALPHABET)

    SequenceColumn(String name, List<String> metadata) {
        super(name, metadata, ColumnType.Sequence)
    }

    @Override
    void append(Entry entry) {
        if (entry.value.length() > 0) {
            def seq = convert(entry.value)
            def entries = stm.get(seq)
            if (entries == null) {
                stm.put(seq, entries = new ArrayList<Entry>())
            }
            entries.add(entry)
        }
    }

    Map<Row, SequenceSearchResult> search(SequenceSearchParameters parameters) {
        def results = new HashMap<Row, SequenceSearchResult>()
        def ni = stm.getNeighborhoodIterator(parameters.query,
                parameters.treeSearchParameters)

        // This hack excludes duplicates like
        //
        //   KLFF     KLFF     ...
        //   KLF-     KL-F
        def prevCdr = new HashSet<String>()

        int i = 0
        while (true) {
            def entries = ni.next()

            if (entries) {
                def seq = entries[0].value
                if (!prevCdr.contains(seq)) {
                    entries.each { entry ->
                        results.put(entry.row, new SequenceSearchResult(ni.currentAlignment, ni.penalty))
                    }
                    prevCdr.add(seq)
                }
            } else {
                break
            }

            if (parameters.depth > 0 && ++i == parameters.depth)
                break
        }

        results
    }

    private static AminoAcidSequence convert(String aaSeq) {
        if (aaSeq =~ /^[FLSYCWPHQRIMTNKVADEG]+$/)
            return new AminoAcidSequence(aaSeq)
        throw new RuntimeException("Illegal character in amino acid sequence string $aaSeq")
    }
}
