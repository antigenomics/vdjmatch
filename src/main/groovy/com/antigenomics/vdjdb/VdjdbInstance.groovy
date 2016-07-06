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


package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.db.ExpressionFilterBatch
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.scoring.SequenceSearcherPreset
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.LevelFilter
import com.antigenomics.vdjdb.text.TextColumn
import com.antigenomics.vdjdb.text.TextFilter

/**
 * A simple API to operate with current VDJdb database.
 * The database is loaded from local database copy. In case it doesn't exist,
 * the latest database release is downloaded from GitHub.
 *
 * This class contains methods to create, filter and wrap the database.
 *
 * @see <a href="https://github.com/antigenomics/vdjdb-db">vdjdb-db repository</a>
 */
class VdjdbInstance {
    static final String SCORE_COLUMN_DEFAULT = "vdjdb.score"

    final Database dbInstance

    VdjdbInstance(boolean useFatDb = true) {
        this(checkDbAndGetMetadata(useFatDb),
                new FileInputStream(Util.HOME_DIR + (useFatDb ? "/vdjdb.txt" : "/vdjdb.slim.txt")))
    }

    VdjdbInstance(Database dbInstance) {
        this.dbInstance = dbInstance
    }

    private static InputStream checkDbAndGetMetadata(boolean useFatDb) {
        Util.checkDatabase()
        new FileInputStream(Util.HOME_DIR + (useFatDb ? "/vdjdb.meta.txt" : "/vdjdb.slim.meta.txt"))
    }

    VdjdbInstance(InputStream metadata, InputStream entries) {
        dbInstance = new Database(metadata)
        dbInstance.addEntries(entries)
    }

    VdjdbInstance(File metadata, File entries) {
        this(new FileInputStream(metadata), new FileInputStream(entries))
    }

    /**
     * Gets the header of current VDJdb database.
     * @return a list of columns
     */
    List<Column> getHeader() {
        dbInstance.columns.collect {
            it instanceof SequenceColumn ? new SequenceColumn(it.name, it.metadata) :
                    new TextColumn(it.name, it.metadata)
        }
    }

    /**
     * Creates an instance of current VDJdb database and applies text and sequence filters if specified.
     * @param textFilters text filters to apply
     * @param sequenceFilters sequence filters to apply
     * @return a VdjdbInstance object
     */
    VdjdbInstance filter(List<TextFilter> textFilters = [],
                         List<SequenceFilter> sequenceFilters = []) {
        new VdjdbInstance(Database.create(dbInstance.search(textFilters, sequenceFilters), dbInstance))
    }

    /**
     * Creates an instance of current VDJdb database and applies text and sequence filters if specified.
     * Pre-filtering is performed using a runtime-evaluated logical expression, containing database column
     * names highlighted with '__', e.g. {@code __source__=~/(EBV|influenza)/} or
     * {@code __source__=="EBV" || __source__=="influenza"}.
     * @param source a file with database table
     * @param expression a logical expression that will be compiled to filter or (String)null
     * @return a VdjdbInstance object
     */
    VdjdbInstance filter(String expression) {
        new VdjdbInstance(Database.create(dbInstance.search(new ExpressionFilterBatch(dbInstance, expression), []), dbInstance))
    }

    /**
     * Converts database instance to a clonotype database, providing means for
     * browsing with {@link com.antigenomics.vdjtools.sample.Clonotype} and
     * {@link com.antigenomics.vdjtools.sample.Sample} vdjtools objects that facilitates
     * searching of specific clonotypes in database.
     *
     * Clonotype searcher parameters and filtering of database records based on VDJdb confidence score can be set here.
     *
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param searchParameters cdr3 matching parameter preset
     * @param species species name
     * @param gene receptor gene name
     * @param vdjdbRecordConfidenceThreshold VDJdb record confidence score threshold
     * @return a clonotype database object
     */
    ClonotypeDatabase asClonotypeDatabase(boolean matchV = false, boolean matchJ = false,
                                          SequenceSearcherPreset searchParameters = SequenceSearcherPreset.byName("dummy"),
                                          String species = null, String gene = null,
                                          int vdjdbRecordConfidenceThreshold = -1) {
        def cdb = new ClonotypeDatabase(header, matchV, matchJ, searchParameters)

        def filters = []

        if (species) {
            filters << new ExactTextFilter(cdb.speciesColName, species, false)
        }
        if (gene) {
            filters << new ExactTextFilter(cdb.geneColName, gene, false)
        }
        if (vdjdbRecordConfidenceThreshold > 0) {
            filters << new LevelFilter(SCORE_COLUMN_DEFAULT, vdjdbRecordConfidenceThreshold.toString(), false)
        }

        cdb.addEntries(dbInstance.rows.collect { row -> row.entries.collect { it.value } }, filters)

        cdb
    }
}