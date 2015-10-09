package com.antigenomics.vdjdb;

import com.antigenomics.vdjdb.core.db.CdrEntry;
import com.antigenomics.vdjdb.filters.Filter;
import com.antigenomics.vdjdb.models.CdrEntrySetDB;
import com.antigenomics.vdjdb.models.EntryDB;

import java.sql.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by bvdmitri on 23.09.15.
 */

public class DatabaseSearcher {
    private final static String url = "jdbc:postgresql://127.0.0.1:5432/";
    private final String dbName;
    private final String user;
    private final String password;

    private Connection connection = null;
    private Statement statement = null;

    public DatabaseSearcher(String dbName, String user, String password) {
        this.dbName = dbName;
        this.user = user;
        this.password = password;
    }

    public List<CdrEntrySetDB> search(List<Filter> filters) throws SQLException {
        List<Filter> setFilters = new ArrayList<>();
        List<Filter> entryFilters = new ArrayList<>();
        for (Filter filter : filters) {
            String fieldName = filter.getFieldName();
            String fieldValue = filter.getFieldValue();
            if (Objects.equals(fieldName, "cdr3")) {
                setFilters.add(filter);
            } else {
                try {
                    EntryDB.Fields field = EntryDB.Fields.valueOf(fieldName.toUpperCase());
                    if (field == null) {
                        System.out.println("Skipping pattern filter [" + fieldName + "," + fieldValue + "]: invalid fieldName");
                    } else {
                        entryFilters.add(filter);
                    }
                } catch (Exception ignored) {
                    System.out.println("Skipping pattern filter [" + fieldName + "," + fieldValue + "]: invalid fieldName");
                }
            }
        }
        if (setFilters.size() > 0 && entryFilters.size() > 0)
            return findSetAndEntries(setFilters, entryFilters);
        else if (entryFilters.size() > 0)
            return findEntries(filters);
        else
            return findSet(filters);
    }

    public List<CdrEntrySetDB> search(Filter ... filters) throws SQLException {
        return search(Arrays.asList(filters));
    }

    public List<CdrEntrySetDB> findEntries(List<Filter> filters) throws SQLException {
        ResultSet entryResult = null;
        ResultSet setResult = null;
        List<CdrEntrySetDB> cdrEntrySetDBList = new ArrayList<>();
        List<EntryDB> list = new ArrayList<EntryDB>();
        StringBuilder query = new StringBuilder("select * from " + EntryDB.DB_TABLE_NAME);
        if (filters.size() > 0) query.append(" where");
        for (int i = 0; i < filters.size(); i++) {
            query.append(filters.get(i).getStatement());
            if (i != filters.size() - 1) query.append(" and ");
        }
        entryResult = statement.executeQuery(query.toString());
        while (entryResult.next()) {
            list.add(new EntryDB(entryResult, null));
        }
        for (EntryDB entryDB : list) {
            setResult = statement.executeQuery("select * from " + CdrEntrySetDB.DB_TABLE_NAME + " where id = '" + entryDB.parentId + "'");
            if (setResult.next()) {
                CdrEntrySetDB parent = new CdrEntrySetDB(setResult);
                parent.addEntry(entryDB);
                entryDB.setParent(parent);
                cdrEntrySetDBList.add(parent);
            }
        }
        entryResult.close();
        if (setResult != null) {
            setResult.close();
        }
        return cdrEntrySetDBList;
    }

    public List<CdrEntrySetDB> findEntries(Filter ... filters) throws SQLException {
        return findEntries(Arrays.asList(filters));
    }

    public List<CdrEntrySetDB> findSet(List<Filter> filters) throws SQLException {
        ResultSet cdrEntrySetResult = null;
        ResultSet entryResult = null;
        List<CdrEntrySetDB> list = new ArrayList<CdrEntrySetDB>();
        StringBuilder query = new StringBuilder("select * from " + CdrEntrySetDB.DB_TABLE_NAME);
        if (filters.size() > 0) query.append(" where");
        for (int i = 0; i < filters.size(); i++) {
            query.append(filters.get(i).getStatement());
            if (i != filters.size() - 1) query.append(" and ");
        }
        cdrEntrySetResult = statement.executeQuery(query.toString());
        while (cdrEntrySetResult.next()) {
            list.add(new CdrEntrySetDB(cdrEntrySetResult));
        }
        for (CdrEntrySetDB entrySetDB : list) {
            entryResult = statement.executeQuery("select * from " + EntryDB.DB_TABLE_NAME + " where parent_id = '" + entrySetDB.getId() + "'");
            while (entryResult.next()) {
                entrySetDB.addEntry(new EntryDB(entryResult, entrySetDB));
            }
        }
        cdrEntrySetResult.close();
        if (entryResult != null) {
            entryResult.close();
        }
        return list;
    }

    public List<CdrEntrySetDB> findSet(Filter ... filters) throws SQLException {
        return findSet(Arrays.asList(filters));
    }

    public List<CdrEntrySetDB> findSetAndEntries(List<Filter> setFilters, List<Filter> entryFilters) throws SQLException {
        ResultSet cdrEntrySetResult = null;
        ResultSet entryResult = null;
        List<CdrEntrySetDB> list = new ArrayList<CdrEntrySetDB>();
        StringBuilder query = new StringBuilder("select * from " + CdrEntrySetDB.DB_TABLE_NAME);
        if (setFilters != null) {
        if (setFilters.size() > 0) query.append(" where");
            for (int i = 0; i < setFilters.size(); i++) {
                query.append(setFilters.get(i).getStatement());
                if (i != setFilters.size() - 1) query.append(" and ");
            }
        }
        cdrEntrySetResult = statement.executeQuery(query.toString());
        while (cdrEntrySetResult.next()) {
            list.add(new CdrEntrySetDB(cdrEntrySetResult));
        }
        List<CdrEntrySetDB> filteredList = new ArrayList<CdrEntrySetDB>();
        for (CdrEntrySetDB entrySetDB : list) {
            StringBuilder queryBuilder = new StringBuilder("select * from ");
            queryBuilder.append(EntryDB.DB_TABLE_NAME)
                    .append(" where parent_id = '")
                    .append(entrySetDB.getId())
                    .append("' ");
            if (entryFilters != null) {
                for (Filter entryFilter : entryFilters) {
                    queryBuilder.append(" and ").append(entryFilter.getStatement());
                }
            }
            entryResult = statement.executeQuery(queryBuilder.toString());
            while (entryResult.next()) {
                entrySetDB.addEntry(new EntryDB(entryResult, entrySetDB));
            }
            if (entrySetDB.getCdrEntries().size() > 0) filteredList.add(entrySetDB);
        }
        cdrEntrySetResult.close();
        if (entryResult != null) {
            entryResult.close();
        }
        return filteredList;
    }

    public void open() throws SQLException {
        connection = DriverManager.getConnection(url + dbName, user, password);
        connection.setReadOnly(true);
        statement = connection.createStatement();
    }

    public void close() throws SQLException {
        if (statement != null) {
            statement.close();
        }
        if (connection != null) {
            connection.close();
        }
    }
}
