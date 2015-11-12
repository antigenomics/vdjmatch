package com.antigenomics.vdjdb.db

enum ColumnType {
    Sequence("seq"), Text("txt")

    final String shortName

    ColumnType(String shortName) {
        this.shortName = shortName
    }

    static ColumnType getByName(String name) {
        name = name.toLowerCase()
        values().find { it.name() == name }
    }
}
