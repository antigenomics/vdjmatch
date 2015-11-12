package com.antigenomics.vdjdb.core.db

enum ColumnType {
    Sequence("seq"), Text("txt")

    final String shortName

    ColumnType(String shortName) {
        this.shortName = shortName
    }

    static ColumnType getByName(String name) {
        name = name.toUpperCase()
        values().find { it.name() == name }
    }
}
