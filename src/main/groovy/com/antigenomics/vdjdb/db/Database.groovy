package com.antigenomics.vdjdb.db

import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.sequence.SequenceSearchResult
import com.antigenomics.vdjdb.text.TextColumn
import com.antigenomics.vdjdb.text.TextFilter
import groovy.transform.CompileStatic

class Database {
    static final String NAME_COL = "name", TYPE_COL = "type", SEQ_TYPE_METADATA_ENTRY = "seq"

    final List<Row> rows = new ArrayList<>()
    final List<Column> columns = new ArrayList<>()

    private final Map<String, Integer> columnId2Index = new HashMap<>()

    static Database create(List<DatabaseSearchResult> searchResults, Database template = null) {
        if(searchResults.size() == 0 && template == null)
            throw new RuntimeException("Cannot create database, empty search results and template database not specified")
        template = template ?: searchResults[0].row.parent
        
        def database = new Database(template.columns.collect {
            it instanceof SequenceColumn ? new SequenceColumn(it.name, it.metadata) :
                    new TextColumn(it.name, it.metadata)
        })

        database.addEntries(searchResults.collect { result -> result.row.entries.collect { it.value } })
        
        database
    }

    public Database(List<Column> columns) {
        columns.each {
            if (columnId2Index.containsKey(it.name)) {
                throw new RuntimeException("Column names should be unique")
            }
            columnId2Index.put(it.name, this.columns.size())
            this.columns.add(it)
        }

        checkColumns()
    }

    public Database(InputStream metadata) {
        boolean first = true
        def metadataField2Index = new HashMap<String, Integer>()

        metadata.splitEachLine("\t") { List<String> splitLine ->
            if (first) {
                for (int i = 0; i < splitLine.size(); i++) {
                    if (metadataField2Index.containsKey(splitLine[i])) {
                        throw new RuntimeException("Repetitive metadata fields are not allowed")
                    }
                    metadataField2Index.put(splitLine[i], i)
                }
                if (!metadataField2Index.containsKey(NAME_COL)) {
                    throw new RuntimeException("Name field is missing from metadata")
                }
                if (!metadataField2Index.containsKey(TYPE_COL)) {
                    throw new RuntimeException("Column type field is missing from metadata")
                }
                first = false
            } else {
                String name = splitLine[metadataField2Index[NAME_COL]],
                       type = splitLine[metadataField2Index[TYPE_COL]]
                if (columnId2Index.containsKey(name)) {
                    throw new RuntimeException("Column names should be unique")
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

    protected boolean checkColumns() {

    }

    @CompileStatic
    public void addEntries(InputStream source, List<TextFilter> filters = []) {
        def filterColIds = getFilterColIds(filters)

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
                    throw new RuntimeException("Some columns specified in database table (${columnSet}) " +
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
                if (pass(filters, filterColIds, row)) {
                    row.entries.eachWithIndex { Entry it, Integer ind -> columns[ind].append(it) }
                    rows.add(row)
                }
            }
        }
    }

    @CompileStatic
    public void addEntries(List<List<String>> entries, List<TextFilter> filters = []) {
        def filterColIds = getFilterColIds(filters)

        entries.each { List<String> splitLine ->
            if (splitLine.size() != columns.size()) {
                throw new RuntimeException("Row size and number of columns don't match")
            }

            def row = new Row(this)

            // Fill row
            for (int i = 0; i < splitLine.size(); i++) {
                row.entries[i] = new Entry(columns[i], row, splitLine[i])
            }

            // Add to database if passes filters
            if (pass(filters, filterColIds, row)) {
                row.entries.eachWithIndex { Entry it, Integer ind -> columns[ind].append(it) }
                rows.add(row)
            }
        }
    }

    @CompileStatic
    private int[] getFilterColIds(List<? extends Filter> filters) {
        def filterColIds = new int[filters.size()]

        for (int i = 0; i < filters.size(); i++) {
            filterColIds[i] = getColumnIndexAndCheck(filters[i])
        }

        filterColIds
    }

    @CompileStatic
    private static boolean pass(List<TextFilter> filters, int[] filterColIds, Row row) {
        for (int i = 0; i < filters.size(); i++) {
            if (!filters[i].pass((Entry) row[filterColIds[i]])) {
                return false
            }
        }
        true
    }

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

    @CompileStatic
    List<DatabaseSearchResult> search(
            List<TextFilter> textFilters,
            List<SequenceFilter> sequenceFilters) {
        def textFilterColIds = getFilterColIds(textFilters),
            sequenceFilterColIds = getFilterColIds(sequenceFilters)

        List<DatabaseSearchResult> results

        if (sequenceFilters.empty) {
            results = rows.findAll { pass(textFilters, textFilterColIds, it) }
                    .collect { new DatabaseSearchResult(it, new SequenceSearchResult[0]) }
        } else {
            results = new ArrayList<>()

            def sequenceSearchResults = new Map<Row, SequenceSearchResult>[sequenceFilters.size()]

            for (int i = 0; i < sequenceFilters.size(); i++) {
                sequenceSearchResults[i] = ((SequenceColumn) columns[sequenceFilterColIds[i]]).search(sequenceFilters[i])
            }

            def minRowSet = sequenceSearchResults.min { it.keySet().size() }.keySet()

            for (Row row : minRowSet) {
                if (pass(textFilters, textFilterColIds, row)) {
                    def sequenceSearchResultsByRow = new SequenceSearchResult[sequenceFilters.size()]

                    for (int i = 0; i < sequenceFilters.size(); i++) {
                        def sequenceSearchResult = sequenceSearchResults[i][row]
                        if (!sequenceSearchResult)
                            continue
                        sequenceSearchResultsByRow[i] = sequenceSearchResult
                    }

                    results.add(new DatabaseSearchResult(row, sequenceSearchResultsByRow))
                }
            }
        }

        results
    }

    int getColumnId(String name) {
        columnId2Index[name]
    }

    Column getAt(String name) {
        columns[getColumnId(name)]
    }

    Row getAt(int index) {
        rows[index]
    }
}
