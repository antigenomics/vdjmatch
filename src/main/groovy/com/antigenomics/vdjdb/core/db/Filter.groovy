package com.antigenomics.vdjdb.core.db

interface Filter {
    String getColumnId()

    ColumnType getColumnType()
}
