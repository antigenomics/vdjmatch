package com.antigenomics.vdjdb.core.db

abstract class Row {
    final Entry[] entries

    Row(Entry[] entries) {
        this.entries = entries
    }

    Entry getAt(int index) {
        entries[index]
    }

    abstract Entry getAt(String columnId)
}