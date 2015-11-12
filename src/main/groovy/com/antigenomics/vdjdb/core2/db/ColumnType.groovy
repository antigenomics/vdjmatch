package com.antigenomics.vdjdb.core2.db

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
