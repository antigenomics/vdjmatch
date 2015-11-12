package com.antigenomics.vdjdb.core.db

abstract class Column {
    final String name
    final Map<String, String> metadata
    final ColumnType columnType

    Column(String name, Map<String, String> metadata, ColumnType columnType) {
        this.name = name
        this.metadata = metadata
        this.columnType = columnType
    }

    abstract void append(Entry entry)

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Column column = (Column) o

        if (name != column.name) return false

        return true
    }

    @Override
    int hashCode() {
        return name.hashCode()
    }
}