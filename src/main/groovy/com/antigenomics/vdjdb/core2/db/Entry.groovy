package com.antigenomics.vdjdb.core2.db

class Entry {
    final Column column
    final Row row
    final String value

    Entry(Column column, Row row, String value) {
        this.column = column
        this.row = row
        this.value = value
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Entry entry = (Entry) o

        if (row != entry.row) return false
        if (value != entry.value) return false

        return true
    }

    @Override
    int hashCode() {
        int result
        result = row.hashCode()
        result = 31 * result + value.hashCode()
        return result
    }
}