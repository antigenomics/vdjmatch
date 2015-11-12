package com.antigenomics.vdjdb.db

abstract class Column {
    final String name
    final Map<String, String> metadata

    Column(String name, Map<String, String> metadata) {
        this.name = name
        this.metadata = metadata
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