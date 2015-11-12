package com.antigenomics.vdjdb.core2.impl

import com.antigenomics.vdjdb.core2.db.Column
import com.antigenomics.vdjdb.core2.db.ColumnType
import com.antigenomics.vdjdb.core2.db.Database
import com.antigenomics.vdjdb.core2.sequence.SequenceFilter
import com.antigenomics.vdjdb.core2.text.ExactTextFilter
import com.antigenomics.vdjdb.core2.text.TextFilter
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import com.milaboratory.core.tree.TreeSearchParameters
import groovy.transform.CompileStatic
import groovyx.gpars.GParsPool

@CompileStatic
class ClonotypeDatabase extends Database {
    static final String CDR_COL = "cdr3", V_COL = "v", J_COL = "j"
    final TreeSearchParameters treeSearchParameters
    final int depth

    ClonotypeDatabase(List<Column> columns,
                      int maxMismatches, int maxInsertions, int maxDeletions, int maxMutations, int depth) {
        super(columns)
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    ClonotypeDatabase(InputStream metadata,
                      int maxMismatches, int maxInsertions, int maxDeletions, int maxMutations, int depth) {
        super(metadata)
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
        def filters = new ArrayList<TextFilter>()

        if (matchV) {
            filters.add(new ExactTextFilter(V_COL, clonotype.v, false))
        }
        if (matchJ) {
            filters.add(new ExactTextFilter(J_COL, clonotype.j, false))
        }

        def results = search(filters,
                [new SequenceFilter(CDR_COL, clonotype.cdr3aaBinary, treeSearchParameters, depth)])

        results.collect {
            new ClonotypeSearchResult(clonotype, it.sequenceSearchResults[0], it.row)
        }.sort()
    }

    List<ClonotypeSearchResult> search(Sample sample, boolean matchV, boolean matchJ) {
        List<List<ClonotypeSearchResult>> results
        GParsPool.withPool ExecUtil.THREADS, {
            results = sample.collectParallel { Clonotype clonotype ->
                search(clonotype, matchV, matchJ)
            }
        }
        results.flatten()
    }
}
