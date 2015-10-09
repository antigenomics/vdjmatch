package com.antigenomics.vdjdb.filters;

/**
 * Created by bvdmitri on 24.09.15.
 */

public interface Filter {
    public String getStatement();
    public String getFieldName();
    public String getFieldValue();
}

