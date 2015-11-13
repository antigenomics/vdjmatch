package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.TextColumn
import com.antigenomics.vdjdb.text.TextFilter
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import com.milaboratory.core.tree.TreeSearchParameters
import groovy.transform.CompileStatic
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentHashMap

class ClonotypeDatabase extends Database {
    final static String CDR3_COL_DEFAULT = "cdr3", V_COL_DEFAULT = "v.segm", J_COL_DEFAULT = "j.segm"

    final String cdr3ColName, vColName, jColName
    final TreeSearchParameters treeSearchParameters
    final int depth
    final boolean matchV, matchJ

    ClonotypeDatabase(List<Column> columns, boolean matchV = false, boolean matchJ = false,
                      int maxMismatches = 2, int maxInsertions = 1, int maxDeletions = 1, int maxMutations = 2, int depth = -1,
                      String cdr3ColName = CDR3_COL_DEFAULT, String vColName = V_COL_DEFAULT, String jColName = J_COL_DEFAULT) {
        super(columns)
        this.cdr3ColName = cdr3ColName
        this.vColName = vColName
        this.jColName = jColName
        this.matchV = matchV
        this.matchJ = matchJ
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    ClonotypeDatabase(InputStream metadata, boolean matchV = false, boolean matchJ = false,
                      int maxMismatches = 2, int maxInsertions = 1, int maxDeletions = 1, int maxMutations = 2, int depth = -1,
                      String cdr3ColName = CDR3_COL_DEFAULT, String vColName = V_COL_DEFAULT, String jColName = J_COL_DEFAULT) {
        super(metadata)

        this.cdr3ColName = cdr3ColName
        this.vColName = vColName
        this.jColName = jColName
        this.matchV = matchV
        this.matchJ = matchJ
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
    }

    @Override
    protected boolean checkColumns() {
        columns.any { it.name == cdr3ColName && it instanceof SequenceColumn } &&
                columns.any { it.name == vColName && it instanceof TextColumn } &&
                columns.any { it.name == jColName && it instanceof TextColumn }
    }

    @CompileStatic
    List<ClonotypeSearchResult> search(Clonotype clonotype,
                                       boolean matchV, boolean matchJ) {
        def filters = new ArrayList<TextFilter>()

        if (matchV) {
            filters.add(new ExactTextFilter(vColName, Util.simplifySegmentName(clonotype.v), false))
        }
        if (matchJ) {
            filters.add(new ExactTextFilter(jColName, Util.simplifySegmentName(clonotype.j), false))
        }

        def results = search(filters,
                [new SequenceFilter(cdr3ColName, clonotype.cdr3aaBinary, treeSearchParameters, depth)])

        results.collect {
            new ClonotypeSearchResult(it.sequenceSearchResults[0], it.row)
        }.sort()
    }

    Map<Clonotype, List<ClonotypeSearchResult>> search(Sample sample) {
        def results = new ConcurrentHashMap<Clonotype, List<ClonotypeSearchResult>>()
        GParsPool.withPool ExecUtil.THREADS, {
            sample.eachParallel { Clonotype clonotype ->
                results.put(clonotype, search(clonotype, matchV, matchJ))
            }
        }
        results
    }
}
