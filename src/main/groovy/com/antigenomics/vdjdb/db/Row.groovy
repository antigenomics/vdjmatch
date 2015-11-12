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
        entries[parent.getColumnId(name)]
    }

    boolean equals(o) {
        index == ((Row) o).index
    }

    int hashCode() {
        index
    }
}