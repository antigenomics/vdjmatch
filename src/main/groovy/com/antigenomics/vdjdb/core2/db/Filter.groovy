package com.antigenomics.vdjdb.core2.db

interface Filter {
    String getColumnId()

    ColumnType getColumnType()
}
