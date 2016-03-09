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
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.TextColumn
import com.antigenomics.vdjdb.text.TextFilter

/**
 * A simple API to operate with current VDJdb database.
 * The database is loaded from local database copy. In case it doesn't exist,
 * latest database release is downloaded using GitHub API.
 *
 * This class contains methods to create, filter and wrap the database.
 *
 * @see <a href="https://github.com/antigenomics/vdjdb-db">vdjdb-db repository</a>
 */
class VdjdbInstance {
    final Database dbInstance

    VdjdbInstance() {
        this(checkDbAndGetMetadata(), new FileInputStream(Util.HOME_DIR + "/vdjdb.txt"))
    }

    private static InputStream checkDbAndGetMetadata() {
        Util.checkDatabase()
        new FileInputStream(Util.HOME_DIR + "/vdjdb.meta.txt")
    }

    VdjdbInstance(InputStream metadata, InputStream entries) {
        dbInstance = new Database(metadata)
        dbInstance.addEntries(entries)
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
     * @return a database instance
     */
    Database filter(List<TextFilter> textFilters = [],
                    List<SequenceFilter> sequenceFilters = []) {
        Database.create(dbInstance.search(textFilters, sequenceFilters), dbInstance)
    }

    /**
     * Creates an instance of current VDJdb database and applies text and sequence filters if specified.
     * @param db database to filter
     * @param textFilters text filters to apply
     * @param sequenceFilters sequence filters to apply
     * @return a database instance
     */
    static Database filter(Database db, List<TextFilter> textFilters = [],
                           List<SequenceFilter> sequenceFilters = []) {
        Database.create(db.search(textFilters, sequenceFilters), db)
    }

    /**
     * Converts database instance to a clonotype database, providing means for 
     * browsing with {@link com.antigenomics.vdjtools.sample.Clonotype} and 
     * {@link com.antigenomics.vdjtools.sample.Sample} vdjtools objects that facilitates
     * searching of specific clonotypes in database.
     *
     * Clonotype searcher parameters can be set here. 
     *
     * @param db database to convert
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param maxMismatches max allowed mismatches in CDR3 region when searching
     * @param maxInsertions max allowed insertions in CDR3 region when searching
     * @param maxDeletions max allowed deletions in CDR3 region when searching
     * @param maxMutations max allowed mutations (mismatches or indels) in CDR3 region when searching
     * @param depth sequence tree scanning depth
     *
     * @return a clonotype database object
     */
    ClonotypeDatabase asClonotypeDatabase(boolean matchV = false, boolean matchJ = false,
                                          int maxMismatches = 2, int maxInsertions = 1,
                                          int maxDeletions = 1, int maxMutations = 2, int depth = -1) {
        def cdb = new ClonotypeDatabase(dbInstance.columns, matchV, matchJ,
                maxMismatches, maxInsertions,
                maxDeletions, maxMutations, depth)

        cdb.addEntries(dbInstance.rows.collect { row -> row.entries.collect { it.value } })

        cdb
    }

    /**
     * Converts database instance to a clonotype database, providing means for 
     * browsing with {@link com.antigenomics.vdjtools.sample.Clonotype} and 
     * {@link com.antigenomics.vdjtools.sample.Sample} vdjtools objects that facilitates
     * searching of specific clonotypes in database.
     *
     * Clonotype searcher parameters can be set here. 
     *
     * @param db database to convert
     * @param matchV should Variable segment matching be performed when searching
     * @param matchJ should Joining segment matching be performed when searching
     * @param maxMismatches max allowed mismatches in CDR3 region when searching
     * @param maxInsertions max allowed insertions in CDR3 region when searching
     * @param maxDeletions max allowed deletions in CDR3 region when searching
     * @param maxMutations max allowed mutations (mismatches or indels) in CDR3 region when searching
     * @param depth sequence tree scanning depth
     *
     * @return a clonotype database object
     */
    static ClonotypeDatabase asClonotypeDatabase(Database db, boolean matchV = false, boolean matchJ = false,
                                                 int maxMismatches = 2, int maxInsertions = 1,
                                                 int maxDeletions = 1, int maxMutations = 2, int depth = -1) {
        def cdb = new ClonotypeDatabase(db.columns, matchV, matchJ,
                maxMismatches, maxInsertions,
                maxDeletions, maxMutations, depth)


        cdb.addEntries(db.rows.collect { row -> row.entries.collect { it.value } })

        cdb
    }
}