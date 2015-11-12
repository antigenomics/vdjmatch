package com.antigenomics.vdjdb.core2.db

import com.antigenomics.vdjdb.core2.sequence.SequenceColumn
import com.antigenomics.vdjdb.core2.sequence.SequenceSearchParameters
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
    private final Map<String, Integer> metadataField2Index = new HashMap<>()

    public Database(List<Column> columns, List<List<String>> entries, Map<String, TextFilter> filters = [:]) {
        columns.each {
            if (columnId2Index.containsKey(it.name)) {
                throw new RuntimeException("Column names should be unique")
            }
            columnId2Index.put(it.name, columns.size())
            columns.add(it)
        }

        checkColumns()

        entries.each { List<String> splitLine ->
            if (splitLine.size() != columns.size()) {
                throw new RuntimeException("Row size and number of columns don't match")
            }

            def row = new RowImpl(rows.size())

            // Fill row
            splitLine.eachWithIndex { String it, Integer ind ->
                row.entries[ind] = new Entry(columns[ind], row, it)
            }

            // Add to database if passes filters
            if (filters.every { it.value.pass(row[it.key]) }) {
                row.entries.eachWithIndex { Entry it, Integer ind -> columns[ind].append(it) }
                rows.add(row)
            }
        }
    }

    public Database(InputStream source, InputStream metadata, Map<String, TextFilter> filters = [:]) {
        boolean first = true
        metadata.splitEachLine("\t") { List<String> splitLine ->
            if (first) {
                splitLine.eachWithIndex { String it, Integer ind ->
                    if (metadataField2Index.containsKey(it)) {
                        throw new RuntimeException("Repetitive metadata fields are not allowed")
                    }
                    metadataField2Index.put(it, ind)
                }
                if (!checkMetadata()) {
                    throw new RuntimeException("Critical metadata fields are missing")
                }
                first = false
            } else {
                def name = splitLine[metadataField2Index[NAME_COL]],
                    type = splitLine[metadataField2Index[TYPE_COL]]
                if (columnId2Index.containsKey(name)) {
                    throw new RuntimeException("Column names should be unique")
                }
                def column = ColumnType.getByName(type) == ColumnType.Sequence ?
                        new SequenceColumn(name, splitLine) : new TextColumn(name, splitLine)
                columnId2Index.put(name, columns.size())
                columns.add(column)
            }
        }

        checkColumns()

        first = true

        source.splitEachLine("\t") { List<String> splitLine ->
            def _columnId2Index = new HashMap<String, Integer>()
            if (first) {
                splitLine.eachWithIndex { String it, Integer ind ->
                    _columnId2Index.put(it, ind)
                }
                if (!_columnId2Index.keySet().containsAll(columnId2Index)) {
                    throw new RuntimeException("Some columns specified in metadata are missing in database table")
                }
                first = false
            } else {
                def row = new RowImpl(rows.size())

                // Fill row
                splitLine.eachWithIndex { String it, Integer ind ->
                    row.entries[ind] = new Entry(columns[ind], row, it)
                }

                // Add to database if passes filters
                if (filters.every { it.value.pass(row[it.key]) }) {
                    row.entries.eachWithIndex { Entry it, Integer ind -> columns[ind].append(it) }
                    rows.add(row)
                }
            }
        }
    }

    protected boolean checkMetadata() {
        metadataField2Index.containsKey(NAME_COL) && metadataField2Index.containsKey(TYPE_COL)
    }

    protected boolean checkColumns() {

    }

    Map<Row, Map<String, SequenceSearchResult>> search(
            Map<String, TextFilter> filters,
            Map<String, SequenceSearchParameters> sequenceSearchParameters) {
        Map<Row, Map<String, SequenceSearchResult>> foundRows

        if (sequenceSearchParameters.empty) {
            foundRows = rows.findAll { row -> filters.every { it.value.pass(row[it.key]) } }
                    .collectEntries { [(it): new HashMap<>()] }
        } else { //todo: optimize for one seq search
            def results = sequenceSearchParameters.collectEntries {
                if (!columnId2Index.containsKey(it.key)) {
                    throw new RuntimeException("Column $it.key doesn't exist")
                }
                def column = getColumn(it.key)
                if (column.columnType != ColumnType.Sequence) {
                    throw new RuntimeException("Sequenced-based search for non-sequence column")
                }
                [(it.key): ((SequenceColumn) column).search(it.value)]
            }

            foundRows = new HashMap<Row, Map<String, SequenceSearchResult>>()

            ((Map<String, Map<Row, SequenceSearchResult>>) results).values()
                    .min { Map<Row, SequenceSearchResult> it -> it.size() } // search rows from smaller map
                    .keySet().each { row ->

                // on-fly filter
                if (filters.every { it.value.pass(row[it.key]) }) {
                    // All results for a given row:
                    def rowResults = ((Map<String, Map<Row, SequenceSearchResult>>) results).collectEntries {
                        [(it.key): it.value[row]]
                    }

                    if (!rowResults.any { it.value == null }) {
                        foundRows.put(row, (Map<String, SequenceSearchResult>) rowResults)
                    }
                }
            }
        }

        foundRows
    }

    String getMetadata(String columnName, String metadataField) {
        getColumn(columnName).metadata[metadataField2Index[metadataField]]
    }

    Column getColumn(String name) {
        columns[columnId2Index[name]]
    }

    private class RowImpl extends Row {
        final int index

        RowImpl(int index) {
            super(new Entry[columns.size()])
            this.index = index
        }

        Entry getAt(int index) {
            entries[index]
        }

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
