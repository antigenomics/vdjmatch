package com.antigenomics.vdjdb.db

import groovy.transform.CompileStatic

@CompileStatic
class Row {
    final Database parent
    final int index
    final Entry[] entries

    Row(Database parent) {
        this.parent = parent
        this.index = parent.rows.size()
        this.entries = new Entry[parent.columns.size()]
    }

    Entry getAt(int index) {
        entries[index]
    }

    Entry getAt(String name) {
        entries[parent.getColumnIndex(name)]
    }

    @Override
    boolean equals(o) {
        index == ((Row) o).index
    }

    @Override
    int hashCode() {
        index
    }

    @Override
    String toString() {
        entries.collect { Entry it -> it.value }.join("\t")
    }
}