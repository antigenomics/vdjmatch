package com.antigenomics.vdjdb.filters;

/**
 * Created by bvdmitri on 24.09.15.
 */

public class FuzzyFilter implements Filter{
    private String fieldName;
    private String fieldValue;
    private int distance;
    private int insertionCost = 1;
    private int deletionCost = 1;
    private int substitutionCost = 1;

    public FuzzyFilter(String fieldName, String fieldValue, int distance) {
        this.fieldName = fieldName;
        this.fieldValue = fieldValue;
        this.distance = distance;
    }

    public FuzzyFilter(String fieldName, String fieldValue, int distance, int insertionCost, int deletionCost, int substitutionCost) {
        this.fieldName = fieldName;
        this.fieldValue = fieldValue;
        this.distance = distance;
        this.insertionCost = insertionCost < 0 ? 1000 : insertionCost;
        this.deletionCost = deletionCost < 0 ? 1000 : deletionCost;
        this.substitutionCost = substitutionCost < 0 ? 1000 : substitutionCost;
    }

    /**
     * Levenshtein function of PostgreSQL database
     * levenshtein(text source, text target, int ins_cost, int del_cost, int sub_cost) returns int
     * The cost parameters specify how much to charge for a character insertion, deletion, or substitution, respectively.
     * @return statement for sql query;
     */
    public String getStatement() {
        return " levenshtein(" + fieldName + ",'" + fieldValue + "'," + insertionCost + "," + deletionCost + "," + substitutionCost + ") <= " + distance + " ";
    }
}
