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

package com.antigenomics.vdjdb.db

import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.sequence.SequenceSearchResult
import com.antigenomics.vdjdb.text.TextColumn
import com.antigenomics.vdjdb.text.TextFilter
import groovy.transform.CompileStatic

/**
 * Base class for database implementations 
 */
class Database {
    static final String NAME_COL = "name", TYPE_COL = "type", SEQ_TYPE_METADATA_ENTRY = "seq"

    final List<Row> rows = new ArrayList<>()
    final List<Column> columns = new ArrayList<>()

    protected final Map<String, Integer> columnId2Index = new HashMap<>()

    /**
     * Creates a new database and fills it with database search results
     * @param searchResults database search results
     * @param template database template, only used to create columns in new database if results are empty
     * @return a new database created from search results (deep copy)
     */
    static Database create(List<SearchResult> searchResults, Database template = null) {
        if (searchResults.size() == 0 && template == null)
            throw new RuntimeException("Cannot create database, empty search results and template database not specified")
        template = template ?: searchResults[0].row.parent

        def database = new Database(template.columns.collect {
            it instanceof SequenceColumn ? new SequenceColumn(it.name, it.metadata) :
                    new TextColumn(it.name, it.metadata)
        })

        database.addEntries(searchResults.collect { result -> result.row.entries.collect { it.value } })

        database
    }

    /**
     * Creates an empty database
     * @param columns a list of database columns
     */
    public Database(List<Column> columns) {
        columns.each {
            if (columnId2Index.containsKey(it.name)) {
                throw new RuntimeException("Error creating database. Column names should be unique ($it)")
            }
            columnId2Index.put(it.name, this.columns.size())
            this.columns.add(it)
        }

        checkColumns()
    }

    /**
     * Creates an empty database using plain-text metadata file. 
     * The file should contain {@value #NAME_COL} column with column names and {@value #TYPE_COL} column with column types.
     * The {@value #SEQ_TYPE_METADATA_ENTRY} type specifies an amino acid sequence column, text column is created otherwise. 
     * @param metadata metadata file stream
     */
    public Database(InputStream metadata) {
        boolean first = true
        def metadataField2Index = new HashMap<String, Integer>()

        metadata.splitEachLine("\t") { List<String> splitLine ->
            if (first) {
                for (int i = 0; i < splitLine.size(); i++) {
                    if (metadataField2Index.containsKey(splitLine[i])) {
                        throw new RuntimeException("Error creating database. Repetitive metadata fields are not allowed (${splitLine[i]})")
                    }
                    metadataField2Index.put(splitLine[i], i)
                }
                if (!metadataField2Index.containsKey(NAME_COL)) {
                    throw new RuntimeException("Error creating database. Name field '$NAME_COL' is missing in metadata")
                }
                if (!metadataField2Index.containsKey(TYPE_COL)) {
                    throw new RuntimeException("Error creating database. Column type field '$TYPE_COL' is missing in metadata")
                }
                first = false
            } else {
                String name = splitLine[metadataField2Index[NAME_COL]],
                       type = splitLine[metadataField2Index[TYPE_COL]]
                if (columnId2Index.containsKey(name)) {
                    throw new RuntimeException("Error creating database. Column names should be unique ($name)")
                }
                def columnMetadata = (Map<String, String>) metadataField2Index.collectEntries {
                    [(it.key): splitLine[it.value]]
                }
                columnMetadata.remove(NAME_COL)
                columnMetadata.remove(TYPE_COL)
                def column = type == SEQ_TYPE_METADATA_ENTRY ?
                        new SequenceColumn(name, columnMetadata) : new TextColumn(name, columnMetadata)
                columnId2Index.put(name, columns.size())
                columns.add(column)
            }
        }

        checkColumns()
    }

    /**
     * Override it if you need to check for presence of specific columns
     * @return true if all necessary columns are present, false otherwise
     */
    protected boolean checkColumns() {

    }

    /**
     * Adds database entries from a given file to the database. First line should 
     * contain column names that should contain those specified during database creation, in any order. 
     * @param source a file with database table
     * @param filters a list of text column filters to apply during loading stage
     */
    void addEntries(InputStream source, List<TextFilter> filters = []) {
        addEntries(source, new ColumnwiseFilterBatch(this, filters))
    }

    /**
     * Adds database entries from a given file to the database. First line should 
     * contain column names that should contain those specified during database creation, in any order.
     * Pre-filtering is performed using a runtime-evaluated logical expression, containing database column 
     * names highlighted with '__', e.g. {@code __source__=~/(EBV|influenza)/} or 
     * {@code __source__=="EBV" || __source__=="influenza"}.
     * @param source a file with database table
     * @param expression a logical expression that will be compiled to filter or (String)null
     */
    void addEntries(InputStream source, String expression) {
        addEntries(source, expression ? new ExpressionFilterBatch(this, expression) : DummyFilterBatch.INSTANCE)
    }

    @CompileStatic
    protected void addEntries(InputStream source, FilterBatch filters) {
        def first = true
        def indexConversion = new HashMap<Integer, Integer>()

        source.splitEachLine("\t") { List<String> splitLine ->
            if (first) {
                // Read in column ids from the table, make table <-> columns index map
                def columnSet = new HashSet<String>()
                for (int i = 0; i < splitLine.size(); i++) {
                    def index = columnId2Index[splitLine[i]]
                    if (index != null) {
                        columnSet.add(splitLine[i])
                        indexConversion.put(i, index)
                    }
                }
                if (!columnSet.containsAll(columnId2Index.keySet())) {
                    throw new RuntimeException("Error filling database. Some columns specified in database table (${columnSet}) " +
                            "are missing in metadata (${columnId2Index.keySet()})")
                }
                first = false
            } else {
                def row = new Row(this)

                // Fill row, convert indices
                indexConversion.each {
                    row.entries[it.value] = new Entry(columns[it.value], row, splitLine[it.key])
                }

                // Add to database if passes filters
                if (filters.pass(row)) {
                    row.entries.eachWithIndex { Entry it, Integer ind -> columns[ind].append(it) }
                    rows.add(row)
                }
            }
        }
    }

    /**
     * Adds a matrix of strings (entries) to the database. 
     * Each row should have the number of strings equal to the number of columns in database. 
     * @param entries a matrix of strings
     * @param filters a list of text column filters to apply during loading stage
     */
    public void addEntries(List<List<String>> entries, List<TextFilter> filters = []) {
        addEntries(entries, new ColumnwiseFilterBatch(this, filters))
    }

    /**
     * Adds a matrix of strings (entries) to the database. 
     * Each row should have the number of strings equal to the number of columns in database.
     * Pre-filtering is performed using a runtime-evaluated logical expression, containing database column 
     * names highlighted with '__', e.g. {@code __source__=~/(EBV|influenza)/} or 
     * {@code __source__=="EBV" || __source__=="influenza"}.
     * @param entries a matrix of strings
     * @param expression a logical expression that will be compiled to filter or (String)null
     */
    public void addEntries(List<List<String>> entries, String expression) {
        addEntries(entries, new ExpressionFilterBatch(this, expression))
    }

    @CompileStatic
    protected void addEntries(List<List<String>> entries, FilterBatch filters) {
        entries.each { List<String> splitLine ->
            if (splitLine.size() != columns.size()) {
                throw new RuntimeException("Error filling database. Row size and number of columns don't match: $splitLine")
            }

            def row = new Row(this)

            // Fill row
            for (int i = 0; i < splitLine.size(); i++) {
                row.entries[i] = new Entry(columns[i], row, splitLine[i])
            }

            // Add to database if passes filters
            if (filters.pass(row)) {
                row.entries.eachWithIndex { Entry it, Integer ind -> columns[ind].append(it) }
                rows.add(row)
            }
        }
    }

    @CompileStatic
    protected int[] getFilterColIds(List<? extends Filter> filters) {
        def filterColIds = new int[filters.size()]

        for (int i = 0; i < filters.size(); i++) {
            filterColIds[i] = getColumnIndexAndCheck(filters[i])
        }

        filterColIds
    }

    @CompileStatic
    private int getColumnIndexAndCheck(Filter filter) {
        def colIndex = columnId2Index[filter.columnId]
        if (colIndex == null) {
            throw new RuntimeException("Bad filter: column $filter.columnId doesn't exist")
        }
        if (filter.sequenceFilter && !(columns[colIndex] instanceof SequenceColumn)) {
            throw new RuntimeException("Bad filter: sequence filter can only be applied to sequence column")
        }
        colIndex
    }

    /**
     * Searches a given database using provided list of text column and sequence column filters.
     * Text filters are applied first, then sequence search is run if specified.
     * @param textFilters list of text column filters, could be empty.
     * @param sequenceFilters list of sequence column filters, could be empty.
     * @return list of database search results.
     */
    @CompileStatic
    List<DatabaseSearchResult> search(
            List<TextFilter> textFilters,
            List<SequenceFilter> sequenceFilters) {
        search(new ColumnwiseFilterBatch(this, textFilters), sequenceFilters)
    }

    /**
     * Searches a given database using provided list of text column and sequence column filters.
     * Batch filter is applied first, then sequence search is run if specified.
     * @param filterBatch a batch filter.
     * @param sequenceFilters list of sequence column filters, could be empty.
     * @return list of database search results.
     */
    @CompileStatic
    List<DatabaseSearchResult> search(
            FilterBatch filterBatch,
            List<SequenceFilter> sequenceFilters) {
        def sequenceFilterColIds = getFilterColIds(sequenceFilters)

        List<DatabaseSearchResult> results

        if (sequenceFilters.empty) {
            results = rows.findAll { filterBatch.pass(it) }
                    .collect { new DatabaseSearchResult(it, new SequenceSearchResult[0]) }
        } else {
            results = new ArrayList<>()

            def sequenceSearchResults = new Map<Row, SequenceSearchResult>[sequenceFilters.size()]

            for (int i = 0; i < sequenceFilters.size(); i++) {
                sequenceSearchResults[i] = ((SequenceColumn) columns[sequenceFilterColIds[i]]).search(sequenceFilters[i])
            }

            def minRowSet = sequenceSearchResults.min { it.keySet().size() }.keySet()

            OUTER:
            for (Row row : minRowSet) {
                if (filterBatch.pass(row)) {
                    def sequenceSearchResultsByRow = new SequenceSearchResult[sequenceFilters.size()]

                    for (int i = 0; i < sequenceFilters.size(); i++) {
                        def sequenceSearchResult = sequenceSearchResults[i][row]
                        if (!sequenceSearchResult)
                            continue OUTER
                        sequenceSearchResultsByRow[i] = sequenceSearchResult
                    }

                    results.add(new DatabaseSearchResult(row, sequenceSearchResultsByRow))
                }
            }
        }

        results
    }

    /**
     * Gets a column index by column identifier 
     * @param name column identifier
     * @return column index
     */
    int getColumnIndex(String name) {
        def index = columnId2Index[name]
        if (index == null)
            throw new RuntimeException("Column $name not found")
        index
    }

    /**
     * Checks if the database has a column with a given identifier 
     * @param name column identifier
     * @return true if specified column is in the database, false otherwise
     */
    boolean hasColumn(String name) {
        columnId2Index.containsKey(name)
    }

    /**
     * Gets an entire column by identifier 
     * @param name column identifier
     * @return database column
     */
    Column getAt(String name) {
        columns[getColumnIndex(name)]
    }

    /**
     * Gets an entire row by index 
     * @param index row index
     * @return database row
     */
    Row getAt(int index) {
        rows[index]
    }

    /**
     * Gets database header
     * @return tab-separated column identifier string
     */
    String getHeader() {
        columns.collect { it.name }.join("\t")
    }

    @Override
    String toString() {
        "columns: ${columns.collect { it.name }.join(",")}, rows: ${rows.size()}"
    }
}
