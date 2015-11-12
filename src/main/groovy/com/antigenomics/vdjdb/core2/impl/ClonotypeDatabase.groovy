package com.antigenomics.vdjdb.core2.impl

import com.antigenomics.vdjdb.core2.db.Column
import com.antigenomics.vdjdb.core2.db.ColumnType
import com.antigenomics.vdjdb.core2.db.Database
import com.antigenomics.vdjdb.core2.sequence.SequenceSearchParameters
import com.antigenomics.vdjdb.core2.text.ExactTextFilter
import com.antigenomics.vdjdb.core2.text.TextFilter
import com.antigenomics.vdjtools.sample.Clonotype
import com.milaboratory.core.tree.TreeSearchParameters
import groovy.transform.CompileStatic

@CompileStatic
class ClonotypeDatabase extends Database {
    static final String CDR_COL = "cdr3", V_COL = "v", J_COL = "j"
    final TreeSearchParameters treeSearchParameters
    final int depth

    ClonotypeDatabase(List<Column> columns, List<List<String>> entries, Map<String, TextFilter> filters,
                      int maxMismatches, int maxInsertions, int maxDeletions, int maxMutations, int depth) {
        super(columns, entries, filters)
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    ClonotypeDatabase(InputStream source, InputStream metadata, Map<String, TextFilter> filters,
                      int maxMismatches, int maxInsertions, int maxDeletions, int maxMutations, int depth) {
        super(source, metadata, filters)
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    @Override
    protected boolean checkColumns() {
        columns.any { it.name == CDR_COL && it.columnType == ColumnType.Sequence } &&
                columns.any { it.name == V_COL && it.columnType == ColumnType.Text } &&
                columns.any { it.name == J_COL && it.columnType == ColumnType.Text }
    }

    List<ClonotypeSearchResult> search(Clonotype clonotype,
                                       boolean matchV, boolean matchJ) {
        def filterMap = new HashMap<String, TextFilter>()

        if (matchV) {
            filterMap.put(V_COL, new ExactTextFilter(clonotype.v, false))
        }
        if (matchJ) {
            filterMap.put(J_COL, new ExactTextFilter(clonotype.j, false))
        }

        def results = search(filterMap,
                [(CDR_COL): new SequenceSearchParameters(clonotype.cdr3aaBinary, treeSearchParameters, depth)])

        results.collect {
            new ClonotypeSearchResult(clonotype, it.value.values().first(), it.key)
        }.sort()
    }
}
