package com.antigenomics.vdjdb.filters;

/**
 * Created by bvdmitri on 24.09.15.
 */
public class MatchFilter implements Filter {
    private String fieldName;
    private String fieldValue;
    private boolean match = true;

    public MatchFilter(String fieldName, String fieldValue) {
        this.fieldName = fieldName;
        this.fieldValue = fieldValue;
    }

    public MatchFilter(String fieldName, String fieldValue, boolean match) {
        this.fieldName = fieldName;
        this.fieldValue = fieldValue;
        this.match = match;
    }

    public void setMatch(boolean match) {
        this.match = match;
    }

    public String getStatement() {
        String matchStatement = match ? "='" : "!='";
        return " " + fieldName + matchStatement + fieldValue + "'";
    }

    @Override
    public String getFieldName() {
        return fieldName;
    }

    @Override
    public String getFieldValue() {
        return fieldValue;
    }
}
