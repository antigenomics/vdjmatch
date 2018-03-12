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

import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.impl.model.AggregateScoring
import com.antigenomics.vdjdb.impl.model.DummyAggregateScoring
import com.antigenomics.vdjdb.impl.segment.DummySegmentScoring
import com.antigenomics.vdjdb.impl.segment.SegmentScoring
import com.antigenomics.vdjdb.sequence.AlignmentScoring
import com.antigenomics.vdjdb.sequence.DummyAlignmentScoring
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.SegmentFilter
import com.antigenomics.vdjdb.text.TextColumn
import com.antigenomics.vdjdb.text.TextFilter
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.milaboratory.core.sequence.AminoAcidSequence
import groovy.transform.CompileStatic
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentHashMap

/**
 * A database implementation holding clonotypes, that is unpaired antigen (T-cell) receptor gene sequences
 */
class ClonotypeDatabase extends Database {
    final static String CDR3_COL_DEFAULT = "cdr3", V_COL_DEFAULT = "v.segm", J_COL_DEFAULT = "j.segm",
                        SPECIES_COL_DEFAULT = "species", GENE_COL_DEFAULT = "gene"

    final String cdr3ColName, vColName, jColName, speciesColName, geneColName
    final int vColIdx, jColIdx, cdr3ColIdx // speed up fetching corresponding entries
    final SearchScope searchParameters
    final AggregateScoring aggregateScoring
    final SegmentScoring segmentScoring
    final AlignmentScoring alignmentScoring
    final boolean matchV, matchJ

    /**
     * Creates an empty clonotype database 
     * @param columns a list of database columns
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param searchParameters search parameter set, sets the search scope
     * @param alignmentScoring alignment scoring function
     * @param segmentScoring segment scoring function
     * @param aggregateScoring score aggregating function
     * @param cdr3ColName CDR3 containing column name
     * @param vColName Variable segment containing column name
     * @param jColName Joining segment containing column name
     * @param speciesColName species column name
     * @param geneColName receptor gene column name
     */
    // todo: filter & informativeness factory
    ClonotypeDatabase(List<Column> columns, boolean matchV = false, boolean matchJ = false,
                      SearchScope searchParameters = SearchScope.EXACT,
                      AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                      SegmentScoring segmentScoring = DummySegmentScoring.INSTANCE,
                      AggregateScoring aggregateScoring = DummyAggregateScoring.INSTANCE,
                      String cdr3ColName = CDR3_COL_DEFAULT, String vColName = V_COL_DEFAULT, String jColName = J_COL_DEFAULT,
                      String speciesColName = SPECIES_COL_DEFAULT, String geneColName = GENE_COL_DEFAULT) {
        super(columns)

        this.cdr3ColName = cdr3ColName
        this.vColName = vColName
        this.jColName = jColName
        this.matchV = matchV
        this.matchJ = matchJ
        this.searchParameters = searchParameters
        this.alignmentScoring = alignmentScoring
        this.segmentScoring = segmentScoring
        this.aggregateScoring = aggregateScoring
        this.speciesColName = speciesColName
        this.geneColName = geneColName
        this.vColIdx = columns.findIndexOf { it.name == vColName }
        this.jColIdx = columns.findIndexOf { it.name == jColName }
        this.cdr3ColIdx = columns.findIndexOf { it.name == cdr3ColName }
    }

    /**
     * Creates an empty clonotype database using plain-text metadata file. 
     * The file should contain {@value #NAME_COL} column with column names and {@value #TYPE_COL} column with column types.
     * The {@value #SEQ_TYPE_METADATA_ENTRY} type specifies an amino acid sequence column, text column is created otherwise. 
     * @param metadata metadata file stream
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param searchParameters search parameter set
     * @param alignmentScoring alignment scoring function
     * @param segmentScoring segment scoring function
     * @param aggregateScoring score aggregating function
     * @param cdr3ColName CDR3 containing column name
     * @param vColName Variable segment containing column name
     * @param jColName Joining segment containing column name
     * @param speciesColName species column name
     * @param geneColName receptor gene column name
     */
    ClonotypeDatabase(InputStream metadata, boolean matchV = false, boolean matchJ = false,
                      SearchScope searchParameters = SearchScope.EXACT,
                      AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                      SegmentScoring segmentScoring = DummySegmentScoring.INSTANCE,
                      AggregateScoring aggregateScoring = DummyAggregateScoring.INSTANCE,
                      String cdr3ColName = CDR3_COL_DEFAULT, String vColName = V_COL_DEFAULT, String jColName = J_COL_DEFAULT,
                      String speciesColName = SPECIES_COL_DEFAULT, String geneColName = GENE_COL_DEFAULT) {
        super(metadata)

        this.cdr3ColName = cdr3ColName
        this.vColName = vColName
        this.jColName = jColName
        this.matchV = matchV
        this.matchJ = matchJ
        this.searchParameters = searchParameters
        this.alignmentScoring = alignmentScoring
        this.segmentScoring = segmentScoring
        this.aggregateScoring = aggregateScoring
        this.speciesColName = speciesColName
        this.geneColName = geneColName
        this.vColIdx = columns.findIndexOf { it.name == vColName }
        this.jColIdx = columns.findIndexOf { it.name == jColName }
        this.cdr3ColIdx = columns.findIndexOf { it.name == cdr3ColName }
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
     * Only records corresponding to specified species and gene are retained
     * @param source a file with database table
     * @param species species name
     * @param gene receptor gene name
     */
    void addEntries(InputStream source, String species, String gene) {
        addEntries(source, [new ExactTextFilter(speciesColName, species, false),
                            new ExactTextFilter(geneColName, gene, false)])
    }

    /**
     * Adds a matrix of strings (entries) to the database. 
     * Each row should have the number of strings equal to the number of columns in database. 
     * Only records corresponding to specified species and gene are retained
     * @param entries a matrix of strings
     * @param species species name
     * @param gene receptor gene name
     */
    void addEntries(List<List<String>> entries, String species, String gene) {
        addEntries(entries, [new ExactTextFilter(speciesColName, species, false),
                             new ExactTextFilter(geneColName, gene, false)])
    }

    /**
     * Searches a database for a given clonotype 
     * @param clonotype a clonotype
     * @return clonotype search result
     */
    @CompileStatic
    List<ClonotypeSearchResult> search(Clonotype clonotype, int id = -1) {
        search(clonotype.v, clonotype.j, clonotype.cdr3aaBinary, id)
    }

    /**
     * Searches a database for a given clonotype 
     * @param v clonotype V segment name
     * @param j clonotype J segment name
     * @param cdr3aa clonotype CDR3 amino acid sequence 
     * @return clonotype search result
     */
    @CompileStatic
    List<ClonotypeSearchResult> search(String v, String j, String cdr3aa, int id = -1) {
        search(v, j, new AminoAcidSequence(cdr3aa), id)
    }

    /**
     * Searches a database for a given clonotype 
     * @param v clonotype V segment name
     * @param j clonotype J segment name
     * @param cdr3aa clonotype CDR3 amino acid sequence 
     * @return clonotype search result
     */
    @CompileStatic
    List<ClonotypeSearchResult> search(String v, String j, AminoAcidSequence cdr3aa, int id = -1) {
        def filters = new ArrayList<TextFilter>()

        if (matchV) {
            filters.add(new SegmentFilter(vColName, v))
        }
        if (matchJ) {
            filters.add(new SegmentFilter(jColName, j))
        }

        def results = search(filters,
                [new SequenceFilter(cdr3ColName, cdr3aa, searchParameters)])

        // todo: filter results
        results.collect {
            def hit = it.hits[0]

            def dbV = it.row[vColIdx].value,
                dbJ = it.row[jColIdx].value,
                dbCdr3 = it.row[cdr3ColIdx].value

            def segmentScores = segmentScoring.computeScores(v, dbV, j, dbJ),
                fullScore = aggregateScoring.computeFullScore(segmentScores.vScore,
                        segmentScores.cdr1Score,
                        segmentScores.cdr2Score,
                        hit.alignmentScore,
                        segmentScores.jScore
                )
            new ClonotypeSearchResult(hit, it.row, id, fullScore, 1)
        }.sort()
    }

    /**
     * Searches a database for a given clonotype sample (in parallel)
     * @param sample a sample
     * @return a map containing search results for every clonotype that was found at least once
     */
    Map<Clonotype, List<ClonotypeSearchResult>> search(Sample sample) {
        def results = new ConcurrentHashMap<Clonotype, List<ClonotypeSearchResult>>()
        GParsPool.withPool ExecUtil.THREADS, {
            sample.eachWithIndexParallel { Clonotype clonotype, int id ->
                if (clonotype.coding) {
                    def result = search(clonotype, id)
                    if (!result.empty) {
                        results.put(clonotype, result)
                    }
                }
            }
        }
        results
    }
}
