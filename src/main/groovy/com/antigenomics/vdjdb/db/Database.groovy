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

    protected final Map<String, Integer> columnId2Index = new HashMap<>()

    static Database create(List<DatabaseSearchResult> searchResults, Database template = null) {
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

    protected boolean checkColumns() {

    }

    public void addEntries(InputStream source, List<TextFilter> filters = []) {
        addEntries(source, new ColumnwiseFilterBatch(this, filters))
    }

    public void addEntries(InputStream source, String expression) {
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

    public void addEntries(List<List<String>> entries, List<TextFilter> filters = []) {
        addEntries(entries, new ColumnwiseFilterBatch(this, filters))
    }

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

    @CompileStatic
    List<DatabaseSearchResult> search(
            List<TextFilter> textFilters,
            List<SequenceFilter> sequenceFilters) {
        def textFilterBatch = new ColumnwiseFilterBatch(this, textFilters)
        def sequenceFilterColIds = getFilterColIds(sequenceFilters)

        List<DatabaseSearchResult> results

        if (sequenceFilters.empty) {
            results = rows.findAll { textFilterBatch.pass(it) }
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
                if (textFilterBatch.pass(row)) {
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

    int getColumnIndex(String name) {
        def index = columnId2Index[name]
        if (index == null)
            throw new RuntimeException("Column $name not found")
        index
    }

    boolean hasColumn(String name) {
        columnId2Index.containsKey(name)
    }

    Column getAt(String name) {
        columns[getColumnIndex(name)]
    }

    Row getAt(int index) {
        rows[index]
    }

    String getHeader() {
        columns.collect { it.name }.join("\t")
    }

    @Override
    String toString() {
        "columns: ${columns.collect { it.name }.join(",")}, rows: ${rows.size()}"
    }
}
