package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.ColumnType
import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.TextFilter
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import com.milaboratory.core.tree.TreeSearchParameters
import groovy.transform.CompileStatic
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentHashMap

class ClonotypeDatabase extends Database {
    static final String CDR_COL = "cdr3", V_COL = "v.segm", J_COL = "j.segm"
    final TreeSearchParameters treeSearchParameters
    final int depth

    ClonotypeDatabase(List<Column> columns) {
        this(columns, 2, 1, 1, 2, -1)
    }

    ClonotypeDatabase(List<Column> columns,
                      int maxMismatches, int maxInsertions, int maxDeletions, int maxMutations, int depth) {
        super(columns)
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    ClonotypeDatabase(InputStream metadata) {
        this(metadata, 2, 1, 1, 2, -1)
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

    @CompileStatic
    List<ClonotypeSearchResult> search(Clonotype clonotype,
                                       boolean matchV, boolean matchJ) {
        def filters = new ArrayList<TextFilter>()

        if (matchV) {
            filters.add(new ExactTextFilter(V_COL, Util.simplifySegmentName(clonotype.v), false))
        }
        if (matchJ) {
            filters.add(new ExactTextFilter(J_COL, Util.simplifySegmentName(clonotype.j), false))
        }

        def results = search(filters,
                [new SequenceFilter(CDR_COL, clonotype.cdr3aaBinary, treeSearchParameters, depth)])

        results.collect {
            new ClonotypeSearchResult(it.sequenceSearchResults[0], it.row)
        }.sort()
    }

    Map<Clonotype, List<ClonotypeSearchResult>> search(Sample sample) {
        search(sample, false, false)
    }

    Map<Clonotype, List<ClonotypeSearchResult>> search(Sample sample, boolean matchV, boolean matchJ) {
        def results = new ConcurrentHashMap<Clonotype, List<ClonotypeSearchResult>>()
        GParsPool.withPool ExecUtil.THREADS, {
            sample.eachParallel { Clonotype clonotype ->
                results.put(clonotype, search(clonotype, matchV, matchJ))
            }
        }
        results
    }
}
