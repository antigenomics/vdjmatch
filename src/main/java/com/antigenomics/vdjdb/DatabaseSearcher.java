package com.antigenomics.vdjdb;

import com.antigenomics.vdjdb.filters.Filter;
import com.antigenomics.vdjdb.models.CdrEntrySetDB;
import com.antigenomics.vdjdb.models.EntryDB;

import java.sql.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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

    public List<EntryDB> findEntries(Filter... filters) {
        return findEntries(Arrays.asList(filters));
    }

    public List<EntryDB> findEntries(List<Filter> filters) {
        ResultSet entryResult = null;
        ResultSet setResult = null;
        List<EntryDB> list = new ArrayList<EntryDB>();
        try {
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
                    entryDB.setParent(new CdrEntrySetDB(setResult));
                }
            }
            return list;
        } catch (SQLException ex) {
            Logger lgr = Logger.getLogger(DatabaseSearcher.class.getName());
            lgr.log(Level.SEVERE, ex.getMessage(), ex);
            return null;
        } finally {
            try {
                if (entryResult != null) {
                    entryResult.close();
                }
                if (setResult != null) {
                    setResult.close();
                }
            } catch (SQLException ex) {
                Logger lgr = Logger.getLogger(DatabaseSearcher.class.getName());
                lgr.log(Level.WARNING, ex.getMessage(), ex);
            }
        }
    }

    public List<CdrEntrySetDB> findSet(Filter... filters) {
        return findSet(Arrays.asList(filters));
    }

    public List<CdrEntrySetDB> findSet(List<Filter> filters) {
        ResultSet cdrEntrySetResult = null;
        ResultSet entryResult = null;
        List<CdrEntrySetDB> list = new ArrayList<CdrEntrySetDB>();
        try {
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
            return list;
        } catch (SQLException ex) {
            Logger lgr = Logger.getLogger(DatabaseSearcher.class.getName());
            lgr.log(Level.SEVERE, ex.getMessage(), ex);
            return null;
        } finally {
            try {
                if (cdrEntrySetResult != null) {
                    cdrEntrySetResult.close();
                }
                if (entryResult != null) {
                    entryResult.close();
                }
            } catch (SQLException ex) {
                Logger lgr = Logger.getLogger(DatabaseSearcher.class.getName());
                lgr.log(Level.WARNING, ex.getMessage(), ex);
            }
        }
    }

    public List<CdrEntrySetDB> findSetAndEntries(List<Filter> setFilters, List<Filter> entryFilters) {
        ResultSet cdrEntrySetResult = null;
        ResultSet entryResult = null;
        List<CdrEntrySetDB> list = new ArrayList<CdrEntrySetDB>();
        try {
            StringBuilder query = new StringBuilder("select * from " + CdrEntrySetDB.DB_TABLE_NAME);
            if (setFilters.size() > 0) query.append(" where");
            for (int i = 0; i < setFilters.size(); i++) {
                query.append(setFilters.get(i).getStatement());
                if (i != setFilters.size() - 1) query.append(" and ");
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
                for (Filter entryFilter : entryFilters) {
                    queryBuilder.append(" and ").append(entryFilter.getStatement());
                }
                entryResult = statement.executeQuery(queryBuilder.toString());
                while (entryResult.next()) {
                    entrySetDB.addEntry(new EntryDB(entryResult, entrySetDB));
                }
                if (entrySetDB.getCdrEntries().size() > 0) filteredList.add(entrySetDB);
            }
            return filteredList;
        } catch (SQLException ex) {
            Logger lgr = Logger.getLogger(DatabaseSearcher.class.getName());
            lgr.log(Level.SEVERE, ex.getMessage(), ex);
            return null;
        } finally {
            try {
                if (cdrEntrySetResult != null) {
                    cdrEntrySetResult.close();
                }
                if (entryResult != null) {
                    entryResult.close();
                }
            } catch (SQLException ex) {
                Logger lgr = Logger.getLogger(DatabaseSearcher.class.getName());
                lgr.log(Level.WARNING, ex.getMessage(), ex);
            }
        }
    }

    public void executeSql(String sql) throws SQLException {
        statement.executeQuery(sql);
    }

    public void open() {
        try {
            connection = DriverManager.getConnection(url + dbName, user, password);
            connection.setReadOnly(true);
            statement = connection.createStatement();
        } catch (SQLException ex) {
            Logger logger = Logger.getLogger(DatabaseSearcher.class.getName());
            logger.log(Level.SEVERE, ex.getMessage(), ex);
        }
    }

    public void close() {
        try {
            if (statement != null) {
                statement.close();
            }
            if (connection != null) {
                connection.close();
            }
        } catch (SQLException ex) {
            Logger logger = Logger.getLogger(DatabaseSearcher.class.getName());
            logger.log(Level.WARNING, ex.getMessage(), ex);
        }
    }
}
