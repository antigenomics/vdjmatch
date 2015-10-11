package com.antigenomics.vdjdb.models;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by bvdmitri on 23.09.15.
 */

public class CdrEntrySetDB {
    public static final String DB_TABLE_NAME = "cdr_entry_set_db";
    private static final int ID_FIELD_NUM = 1;
    private static final int CDR_FIELD_NUM = 2;

    public static enum Fields {
        ID("id", "ID",ID_FIELD_NUM),
        CDR3("cdr3", "CDR3",CDR_FIELD_NUM);

        private String fieldName;
        private String name;
        private int columnNum;

        private Fields(String fieldName, String name, int columnNum) {
            this.fieldName = fieldName;
            this.name = name;
            this.columnNum = columnNum;
        }

        public String getFieldName() {
            return fieldName;
        }

        public int getColumnNum() {
            return columnNum;
        }

        public String getName() {
            return name;
        }
    }

    private Long id;
    private String cdr3;
    private List<EntryDB> cdrEntries;

    public CdrEntrySetDB(String cdr3, Long id) {
        this.id = id;
        this.cdr3 = cdr3;
        this.cdrEntries = new ArrayList<EntryDB>();
    }

    public CdrEntrySetDB(ResultSet resultSet) throws SQLException {
        this(resultSet.getString(CDR_FIELD_NUM), resultSet.getLong(ID_FIELD_NUM));
    }

    public Long getId() {
        return id;
    }

    public String getCdr3() {
        return cdr3;
    }

    public void addEntry(EntryDB entryDB) {
        cdrEntries.add(entryDB);
    }

    public List<EntryDB> getCdrEntries() {
        return cdrEntries;
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(cdr3).append("\n");
        for (EntryDB cdrEntry : cdrEntries) {
            stringBuilder.append(cdrEntry.toString());
        }
        return stringBuilder.toString();
    }
}
