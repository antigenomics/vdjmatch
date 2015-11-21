/*
 * Copyright 2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

/**
 * A database implementation holding clonotypes, that is unpaired antigen receptor chains 
 */
class ClonotypeDatabase extends Database {
    final static String CDR3_COL_DEFAULT = "cdr3", V_COL_DEFAULT = "v.segm", J_COL_DEFAULT = "j.segm",
                        SPECIES_COL_DEFAULE = "species", CHAIN_COL_DEFAULT = "chain"

    final String cdr3ColName, vColName, jColName, speciesColName, chainColName
    final TreeSearchParameters treeSearchParameters
    final int depth
    final boolean matchV, matchJ

    /**
     * Creates an empty clonotype database 
     * @param columns a list of database columns
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param maxMismatches max allowed mismatches in CDR3 region when searching
     * @param maxInsertions max allowed insertions in CDR3 region when searching
     * @param maxDeletions max allowed deletions in CDR3 region when searching
     * @param maxMutations max allowed mutations (mismatches or indels) in CDR3 region when searching
     * @param depth sequence tree scanning depth
     * @param cdr3ColName CDR3 containing column name
     * @param vColName Variable segment containing column name
     * @param jColName Joining segment containing column name
     * @param speciesColName species column name
     * @param chainColName receptor chain column name
     */
    ClonotypeDatabase(List<Column> columns, boolean matchV = false, boolean matchJ = false,
                      int maxMismatches = 2, int maxInsertions = 1, int maxDeletions = 1, int maxMutations = 2, int depth = -1,
                      String cdr3ColName = CDR3_COL_DEFAULT, String vColName = V_COL_DEFAULT, String jColName = J_COL_DEFAULT,
                      String speciesColName = SPECIES_COL_DEFAULE, String chainColName = CHAIN_COL_DEFAULT) {
        super(columns)
        this.cdr3ColName = cdr3ColName
        this.vColName = vColName
        this.jColName = jColName
        this.matchV = matchV
        this.matchJ = matchJ
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
        this.speciesColName = speciesColName
        this.chainColName = chainColName
    }

    /**
     * Creates an empty clonotype database using plain-text metadata file. 
     * The file should contain {@value #NAME_COL} column with column names and {@value #TYPE_COL} column with column types.
     * The {@value #SEQ_TYPE_METADATA_ENTRY} type specifies an amino acid sequence column, text column is created otherwise. 
     * @param metadata metadata file stream
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param maxMismatches max allowed mismatches in CDR3 region when searching
     * @param maxInsertions max allowed insertions in CDR3 region when searching
     * @param maxDeletions max allowed deletions in CDR3 region when searching
     * @param maxMutations max allowed mutations (mismatches or indels) in CDR3 region when searching
     * @param depth sequence tree scanning depth
     * @param cdr3ColName CDR3 containing column name
     * @param vColName Variable segment containing column name
     * @param jColName Joining segment containing column name
     * @param speciesColName species column name
     * @param chainColName receptor chain column name
     */
    ClonotypeDatabase(InputStream metadata, boolean matchV = false, boolean matchJ = false,
                      int maxMismatches = 2, int maxInsertions = 1, int maxDeletions = 1, int maxMutations = 2, int depth = -1,
                      String cdr3ColName = CDR3_COL_DEFAULT, String vColName = V_COL_DEFAULT, String jColName = J_COL_DEFAULT,
                      String speciesColName = SPECIES_COL_DEFAULE, String chainColName = CHAIN_COL_DEFAULT) {
        super(metadata)

        this.cdr3ColName = cdr3ColName
        this.vColName = vColName
        this.jColName = jColName
        this.matchV = matchV
        this.matchJ = matchJ
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, maxInsertions, maxDeletions, maxMutations)
        this.depth = depth
        this.speciesColName = speciesColName
        this.chainColName = chainColName
    }

    @Override
    protected boolean checkColumns() {
        columns.any { it.name == cdr3ColName && it instanceof SequenceColumn } &&
                columns.any { it.name == vColName && it instanceof TextColumn } &&
                columns.any { it.name == jColName && it instanceof TextColumn }
    }

    /**
     * Adds database entries from a given file to the database. First line should 
     * contain column names that should contain those specified during database creation, in any order.
     * Only records corresponding to specified species and chain are retained 
     * @param source a file with database table
     * @param species species name
     * @param chain receptor chain name 
     */
    void addEntries(InputStream source, String species, String chain) {
        addEntries(source, [new ExactTextFilter(speciesColName, species, false),
                            new ExactTextFilter(chainColName, chain, false)])
    }

    /**
     * Adds a matrix of strings (entries) to the database. 
     * Each row should have the number of strings equal to the number of columns in database. 
     * Only records corresponding to specified species and chain are retained
     * @param entries a matrix of strings
     * @param species species name
     * @param chain receptor chain name 
     */
    public void addEntries(List<List<String>> entries, String species, String chain) {
        addEntries(entries, [new ExactTextFilter(speciesColName, species, false),
                             new ExactTextFilter(chainColName, chain, false)])
    }

    /**
     * Searches a database for a given clonotype 
     * @param clonotype a clonotype
     * @return clonotype search result
     */
    @CompileStatic
    List<ClonotypeSearchResult> search(Clonotype clonotype) {
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

    /**
     * Searches a database for a given clonotype sample (in parallel)
     * @param clonotype a clonotype sample
     * @return a map containing search results for every clonotype that was found at least once
     */
    Map<Clonotype, List<ClonotypeSearchResult>> search(Sample sample) {
        def results = new ConcurrentHashMap<Clonotype, List<ClonotypeSearchResult>>()
        GParsPool.withPool ExecUtil.THREADS, {
            sample.eachParallel { Clonotype clonotype ->
                def result = search(clonotype)
                if (!result.empty) {
                    results.put(clonotype, result)
                }
            }
        }
        results
    }
}
