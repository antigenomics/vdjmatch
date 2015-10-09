package com.antigenomics.vdjdb.browsing

import com.antigenomics.vdjdb.DatabaseSearcher
import com.antigenomics.vdjdb.filters.Filter
import com.antigenomics.vdjdb.filters.FuzzyFilter
import com.antigenomics.vdjdb.filters.MatchFilter
import com.antigenomics.vdjdb.filters.PatternFilter
import com.antigenomics.vdjdb.models.CdrEntrySetDB
import com.antigenomics.vdjdb.models.EntryDB
import groovyjarjarcommonscli.Option

import java.util.logging.Level
import java.util.logging.Logger

/**
 * Created by bvdmitri on 30.09.15.
 */

def printHeader() {
    EntryDB.printHeader()
}

def printEntries(List<EntryDB> entries) {
    entries.each {
        entry ->
            entry.printEntry()
    }
}

def cli = new CliBuilder(usage: "SearchDatabase [options] ")
cli.h("display help message")
cli.d(longOpt: "database", argName: "string", args: 1, "PostgreSQL database name", required: true)
cli.u(longOpt: "user", argName: "string", args: 1, "PostgreSQL user name", required: true)
cli.p(longOpt: "password", argName: "string", args: 1, "PostgreSQL user password", required: true)
cli.fm(longOpt: "filterMatch", args: Option.UNLIMITED_VALUES, valueSeparator: ',', required: false, "Match Filters")
cli.fp(longOpt: "filterPattern", args: Option.UNLIMITED_VALUES, valueSeparator: ',', required: false, "Pattern Filters")
cli.ff(longOpt: "filterFuzzy", args: Option.UNLIMITED_VALUES, valueSeparator: ',', required: false, "Fuzzy Filters")
cli.sf(longOpt: "showFields", "Show available fields")
cli.e(longOpt: "errors", "Show detailed information about errors")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(-1)
}

if (opt.h) {
    cli.usage()
    System.exit(-1)
}

if (opt.sf) {
    println "Fields: "
    EntryDB.Fields.values().each {
        field ->
            if (!field.fieldName.contains("id"))
                print field.fieldName + " "
    }
    CdrEntrySetDB.Fields.values().each {
        field ->
            if (!field.fieldName.contains("id"))
                print field.fieldName + " "
    }
    println ""
    System.exit(0)
}

def dbName = (String) (opt.d ?: null),
    userName = (String) (opt.u ?: null),
    password = (String) (opt.p ?: null),
    filtersMatch = opt.fms ?: null,
    filtersPattern = opt.fps ?: null,
    filtersFuzzy = opt.ffs ?: null

if (dbName == null || userName == null || password == null) {
    println "Invalid parameters for -d (PostgreSQL database name) -u (PostgreSQL user name) and -p (PostgreSQL use password) options"
    System.exit(-1)
}

List<Filter> entryFilterList = new ArrayList<>()
List<Filter> setFilterList = new ArrayList<>();

filtersMatch.each {
    filter ->
        String f = (String) filter
        String[] rules = f.split("=")
        try {
            String fieldName = rules[0]
            String valueField = rules[1]
            Boolean match = true
            if (rules.length == 3) {
                match = Boolean.parseBoolean(rules[2])
            }
            MatchFilter matchFilter = new MatchFilter(fieldName, valueField, match)
            if (fieldName == "cdr3") {
                setFilterList.add(matchFilter)
            } else {
                try {
                    def field = EntryDB.Fields.valueOf(fieldName.toUpperCase())
                    if (field == null) {
                        println "Skipping match filter [" + fieldName + "," + valueField + "]: invalid fieldName"
                    } else {
                        entryFilterList.add(matchFilter)
                    }
                } catch (Exception ignored) {
                    println "Skipping match filter [" + fieldName + "," + valueField + "]: invalid fieldName"
                }
            }
        } catch (Exception ignored) {
            println "Skipping match filter " + rules + ": must be [<fieldName:string>, <value:string>]"
        }
}

filtersPattern.each {
    filter ->
        String f = (String) filter
        String[] rules = f.split("=")
        try {
            String fieldName = rules[0]
            String valueField = rules[1]
            Boolean match = true
            if (rules.length == 3) {
                match = Boolean.parseBoolean(rules[2])
            }
            PatternFilter matchFilter = new PatternFilter(fieldName, valueField, match)
            if (fieldName == "cdr3") {
                setFilterList.add(matchFilter)
            } else {
                try {
                    def field = EntryDB.Fields.valueOf(fieldName.toUpperCase())
                    if (field == null) {
                        println "Skipping pattern filter [" + fieldName + "," + valueField + "]: invalid fieldName"
                    } else {
                        entryFilterList.add(matchFilter)
                    }
                } catch (Exception ignored) {
                    println "Skipping pattern filter [" + fieldName + "," + valueField + "]: invalid fieldName"
                }
            }
        } catch (Exception ignored) {
            println "Skipping pattern filter " + rules + ": must be [<fieldName:string>, <value:string>]"
        }
}

filtersFuzzy.each {
    filter ->
        String f = (String) filter
        String[] rules = f.split("=")
        try {
            String fieldName = rules[0]
            String valueField = rules[1]
            int distance = Integer.parseInt(rules[2])
            FuzzyFilter fuzzyFilter = new FuzzyFilter(fieldName, valueField, distance)
            if (fieldName == "cdr3") {
                setFilterList.add(fuzzyFilter)
            } else {
                try {
                    def field = EntryDB.Fields.valueOf(fieldName.toUpperCase())
                    if (field == null) {
                        println "Skipping fuzzy filter [" + fieldName + "," + valueField + "]: invalid fieldName"
                    } else {
                        entryFilterList.add(fuzzyFilter)
                    }
                } catch (Exception ignored) {
                    println "Skipping fuzzy filter [" + fieldName + "," + valueField + "]: invalid fieldName"
                }
            }
        } catch (Exception ignored) {
            println "Skipping fuzzy filter " + rules + ": must be [<fieldName:string>, <value:string>, <distance:int>]"
        }
}

DatabaseSearcher databaseSearcher;

try {
    databaseSearcher = new DatabaseSearcher(dbName, userName, password)
    databaseSearcher.open()
} catch (Exception ex) {
    if (opt.e) {
        Logger logger = Logger.getLogger(DatabaseSearcher.class.getName());
        logger.log(Level.WARNING, ex.getMessage(), ex);
    } else {
        println "An error has occurred while accessing PostgreSQL database. Use '-e' option if you want to see detailed information about this error"
    }
    System.exit(-1);
}

try {
    printHeader()
    if (entryFilterList.size() > 0 && setFilterList.size() > 0) {
        List<CdrEntrySetDB> entries = databaseSearcher.findSetAndEntries(setFilterList, entryFilterList)
        for (CdrEntrySetDB cdrEntrySetDB : entries) {
            printEntries(cdrEntrySetDB.cdrEntries)
        }
    } else if (entryFilterList.size() > 0) {
        List<EntryDB> entries = databaseSearcher.findEntries(entryFilterList)
        printEntries(entries)
    } else if (setFilterList.size() > 0) {
        List<CdrEntrySetDB> entries = databaseSearcher.findSet(setFilterList)
        for (CdrEntrySetDB cdrEntrySetDB : entries) {
            printEntries(cdrEntrySetDB.cdrEntries)
        }
    } else {
        List<CdrEntrySetDB> entries = databaseSearcher.findSet()
        for (CdrEntrySetDB cdrEntrySetDB : entries) {
            printEntries(cdrEntrySetDB.cdrEntries)
        }
    }
    databaseSearcher.close()
} catch (Exception ex) {
    if (opt.e) {
        Logger logger = Logger.getLogger(DatabaseSearcher.class.getName());
        logger.log(Level.WARNING, ex.getMessage(), ex);
    } else {
        println "An error has occurred. Use '-e' option if you want to see detailed information about errors"
    }
    System.exit(-1);
}