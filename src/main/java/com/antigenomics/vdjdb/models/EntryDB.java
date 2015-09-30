package com.antigenomics.vdjdb.models;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Created by bvdmitri on 23.09.15.
 */

public class EntryDB {
    public static final String DB_TABLE_NAME = "cdr_entry_db";
    private static final int V_FIELD_NUM = 2;
    private static final int J_FIELD_NUM = 3;
    private static final int CHAIN_FIELD_NUM = 4;
    private static final int MHC_FIELD_NUM = 5;
    private static final int NATURE_FIELD_NUM = 6;
    private static final int DISEASE_FIELD_NUM = 7;
    private static final int ORIGIN_FIELD_NUM = 8;
    private static final int ANTIGEN_SEQ_FIELD_NUM = 9;
    private static final int ANTIGEN_NAME_FIELD_NUM = 10;
    private static final int METHOD_FIELD_NUM = 11;
    private static final int GENBANK_FIELD_NUM = 12;
    private static final int REFERENCE_FIELD_NUM = 13;
    public static final int PARENT_ID_FIELD_NUM = 14;

    public static enum Fields {
        V("v", V_FIELD_NUM),
        J("j", J_FIELD_NUM),
        CHAIN("chain", CHAIN_FIELD_NUM),
        MHC("mhc", MHC_FIELD_NUM),
        NATURE("nature", NATURE_FIELD_NUM),
        DISEASE("disease", DISEASE_FIELD_NUM),
        ORIGIN("origin", ORIGIN_FIELD_NUM),
        ANTIGEN_SEQ("antigen_seq", ANTIGEN_SEQ_FIELD_NUM),
        ANTIGEN_NAME("antigen_name", ANTIGEN_NAME_FIELD_NUM),
        METHOD("method", METHOD_FIELD_NUM),
        GENBANK("genbank", GENBANK_FIELD_NUM),
        REFERENCE("reference", REFERENCE_FIELD_NUM),
        PARENT_ID("parent_id", PARENT_ID_FIELD_NUM);

        private String fieldName;
        private int columnNum;

        private Fields(String fieldName, int columnNum) {
            this.fieldName = fieldName;
            this.columnNum = columnNum;
        }

        public int getColumnNum() {
            return columnNum;
        }

        public String getFieldName() {
            return fieldName;
        }
    }

    public final String v;
    public final String j;
    public final String chain;
    public final String mhc;
    public final String nature;
    public final String disease;
    public final String origin;
    public final String antigen_seq;
    public final String antigen_name;
    public final String method;
    public final String genbank;
    public final String reference;
    public final Long parentId;

    private CdrEntrySetDB parent;

    public EntryDB(ResultSet resultSet, CdrEntrySetDB parent) throws SQLException {
        this.v = resultSet.getString(V_FIELD_NUM);
        this.j = resultSet.getString(J_FIELD_NUM);
        this.chain = resultSet.getString(CHAIN_FIELD_NUM);
        this.mhc = resultSet.getString(MHC_FIELD_NUM);
        this.nature = resultSet.getString(NATURE_FIELD_NUM);
        this.disease = resultSet.getString(DISEASE_FIELD_NUM);
        this.origin = resultSet.getString(ORIGIN_FIELD_NUM);
        this.antigen_seq = resultSet.getString(ANTIGEN_SEQ_FIELD_NUM);
        this.antigen_name = resultSet.getString(ANTIGEN_NAME_FIELD_NUM);
        this.method = resultSet.getString(METHOD_FIELD_NUM);
        this.genbank = resultSet.getString(GENBANK_FIELD_NUM);
        this.reference = resultSet.getString(REFERENCE_FIELD_NUM);
        this.parentId = resultSet.getLong(PARENT_ID_FIELD_NUM);
        this.parent = parent;
    }

    public String getV() {
        return v;
    }

    public String getJ() {
        return j;
    }

    public String getChain() {
        return chain;
    }

    public String getMhc() {
        return mhc;
    }

    public String getNature() {
        return nature;
    }

    public String getDisease() {
        return disease;
    }

    public String getOrigin() {
        return origin;
    }

    public String getAntigen_seq() {
        return antigen_seq;
    }

    public String getAntigen_name() {
        return antigen_name;
    }

    public String getMethod() {
        return method;
    }

    public String getGenbank() {
        return genbank;
    }

    public String getReference() {
        return reference;
    }

    public CdrEntrySetDB getParent() {
        return parent;
    }

    public String getCdr3() {
        return parent.getCdr3();
    }

    public void setParent(CdrEntrySetDB parent) {
        this.parent = parent;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        EntryDB entryDB = (EntryDB) o;

        if (v != null ? !v.equals(entryDB.v) : entryDB.v != null) return false;
        if (j != null ? !j.equals(entryDB.j) : entryDB.j != null) return false;
        if (chain != null ? !chain.equals(entryDB.chain) : entryDB.chain != null) return false;
        if (mhc != null ? !mhc.equals(entryDB.mhc) : entryDB.mhc != null) return false;
        if (nature != null ? !nature.equals(entryDB.nature) : entryDB.nature != null) return false;
        if (disease != null ? !disease.equals(entryDB.disease) : entryDB.disease != null) return false;
        if (origin != null ? !origin.equals(entryDB.origin) : entryDB.origin != null) return false;
        if (antigen_seq != null ? !antigen_seq.equals(entryDB.antigen_seq) : entryDB.antigen_seq != null) return false;
        if (antigen_name != null ? !antigen_name.equals(entryDB.antigen_name) : entryDB.antigen_name != null)
            return false;
        if (method != null ? !method.equals(entryDB.method) : entryDB.method != null) return false;
        if (genbank != null ? !genbank.equals(entryDB.genbank) : entryDB.genbank != null) return false;
        if (reference != null ? !reference.equals(entryDB.reference) : entryDB.reference != null) return false;
        if (parentId != null ? !parentId.equals(entryDB.parentId) : entryDB.parentId != null) return false;
        return !(parent != null ? !parent.equals(entryDB.parent) : entryDB.parent != null);

    }

    @Override
    public int hashCode() {
        int result = v != null ? v.hashCode() : 0;
        result = 31 * result + (j != null ? j.hashCode() : 0);
        result = 31 * result + (chain != null ? chain.hashCode() : 0);
        result = 31 * result + (mhc != null ? mhc.hashCode() : 0);
        result = 31 * result + (nature != null ? nature.hashCode() : 0);
        result = 31 * result + (disease != null ? disease.hashCode() : 0);
        result = 31 * result + (origin != null ? origin.hashCode() : 0);
        result = 31 * result + (antigen_seq != null ? antigen_seq.hashCode() : 0);
        result = 31 * result + (antigen_name != null ? antigen_name.hashCode() : 0);
        result = 31 * result + (method != null ? method.hashCode() : 0);
        result = 31 * result + (genbank != null ? genbank.hashCode() : 0);
        result = 31 * result + (reference != null ? reference.hashCode() : 0);
        result = 31 * result + (parentId != null ? parentId.hashCode() : 0);
        result = 31 * result + (parent != null ? parent.hashCode() : 0);
        return result;
    }

    public void printEntry() {
        String parentString = parent != null ? parent.getCdr3() : " ";
        System.out.printf("%10s %10s %10s %20s %10s %10s %10s %20s %20s %20s %20s %30s %s", v, j, chain, mhc, nature, disease, origin, antigen_seq, antigen_name, method, genbank, parentString, reference);
        System.out.println();
    }

    public static void printHeader() {
        System.out.printf("%10s %10s %10s %20s %10s %10s %10s %20s %20s %20s %20s %30s %s", "v", "j", "chain", "mhc", "nature", "disease", "origin", "antigen_seq", "antigen_name", "method", "genbank", "cdr3", "reference");
        System.out.println();
    }

    @Override
    public String toString() {
        String parentString = parent != null ? parent.getCdr3() : " ";
        return v + "\t" + j + "\t" + chain + "\t" + mhc + "\t" + nature + "\t" + disease + "\t" + origin + "\t" + antigen_seq + "\t" + antigen_name + "\t" + method + "\t" + genbank + "\t" + reference + "\t" + parentString + "\n";
    }
}
