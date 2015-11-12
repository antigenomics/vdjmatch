package com.antigenomics.vdjdb.sequence

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.ColumnType
import com.antigenomics.vdjdb.db.Entry
import com.antigenomics.vdjdb.db.Row
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import groovy.transform.CompileStatic

@CompileStatic
class SequenceColumn extends Column {
    final SequenceTreeMap<AminoAcidSequence, List<Entry>> stm = new SequenceTreeMap(AminoAcidSequence.ALPHABET)

    SequenceColumn(String name, Map<String, String> metadata = [:]) {
        super(name, metadata, ColumnType.Sequence)
    }

    @Override
    void append(Entry entry) {
        if (entry.value.length() > 0) {
            def seq = Util.convert(entry.value)
            def entries = stm.get(seq)
            if (entries == null) {
                stm.put(seq, entries = new ArrayList<Entry>())
            }
            entries.add(entry)
        }
    }

    Map<Row, SequenceSearchResult> search(SequenceFilter filter) {
        def results = new HashMap<Row, SequenceSearchResult>()

        def ni = stm.getNeighborhoodIterator(filter.query,
                filter.treeSearchParameters)

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

            if (filter.depth > 0 && ++i == filter.depth)
                break
        }

        results
    }
}
