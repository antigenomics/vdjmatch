package com.antigenomics.vdjdb.core2.db

import com.antigenomics.vdjdb.core2.sequence.SequenceColumn
import com.antigenomics.vdjdb.core2.sequence.SequenceFilter
import com.antigenomics.vdjdb.core2.sequence.SequenceSearchResult
import com.antigenomics.vdjdb.core2.text.TextColumn
import com.antigenomics.vdjdb.core2.text.TextFilter
import groovy.transform.CompileStatic

@CompileStatic
class Database {
    static final String NAME_COL = "name", TYPE_COL = "type"

    final List<Row> rows = new ArrayList<>()
    final List<Column> columns = new ArrayList<>()

    private final Map<String, Integer> columnId2Index = new HashMap<>()

    public Database(List<Column> columns) {
        columns.each {
            if (columnId2Index.containsKey(it.name)) {
                throw new RuntimeException("Column names should be unique")
            }
            columnId2Index.put(it.name, columns.size())
            columns.add(it)
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
                def column = ColumnType.getByName(type) == ColumnType.Sequence ?
                        new SequenceColumn(name, columnMetadata) : new TextColumn(name, columnMetadata)
                columnId2Index.put(name, columns.size())
                columns.add(column)
            }
        }

        checkColumns()
    }

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
                    throw new RuntimeException("Some columns specified in metadata are missing in database table")
                }
                first = false
            } else {
                def row = new RowImpl(rows.size())

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

    public void addEntries(List<List<String>> entries, List<TextFilter> filters = []) {
        def filterColIds = getFilterColIds(filters)

        entries.each { List<String> splitLine ->
            if (splitLine.size() != columns.size()) {
                throw new RuntimeException("Row size and number of columns don't match")
            }

            def row = new RowImpl(rows.size())

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

    private int[] getFilterColIds(List<Filter> filters) {
        def filterColIds = new int[filters.size()]

        for (int i = 0; i < filters.size(); i++) {
            filterColIds[i] = getColumnIndexAndCheck(filters[i])
        }

        filterColIds
    }

    private boolean pass(List<TextFilter> filters, int[] filterColIds, Row row) {
        for (int i = 0; i < filters.size(); i++) {
            if (!filters[i].pass(row[filterColIds[i]])) {
                return false
            }
        }
        true
    }

    protected boolean checkColumns() {

    }

    private int getColumnIndexAndCheck(Filter filter) {
        def colIndex = columnId2Index[filter.columnId]
        if (colIndex == null) {
            throw new RuntimeException("Bad filter: column $filter.columnId doesn't exist")
        }
        if (columns[colIndex].columnType != filter.columnType) {
            throw new RuntimeException("Bad filter: filter and column type don't match")
        }
        colIndex
    }

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

    private class RowImpl extends Row {
        final int index

        RowImpl(int index) {
            super(new Entry[columns.size()])
            this.index = index
        }

        @Override
        Entry getAt(String columnId) {
            entries[columnId2Index[columnId]]
        }

        @Override
        boolean equals(o) {
            if (this.is(o)) return true
            if (getClass() != o.class) return false

            RowImpl row = (RowImpl) o

            if (index != row.index) return false

            return true
        }

        @Override
        int hashCode() {
            return index
        }
    }
}
