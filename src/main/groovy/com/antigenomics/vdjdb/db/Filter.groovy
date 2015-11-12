package com.antigenomics.vdjdb.db

interface Filter {
    String getColumnId()

    ColumnType getColumnType()
}
