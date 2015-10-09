package com.antigenomics.vdjdb;

import com.antigenomics.vdjdb.filters.MatchFilter;
import com.antigenomics.vdjdb.filters.PatternFilter;
import com.antigenomics.vdjdb.models.CdrEntrySetDB;
import com.antigenomics.vdjdb.models.EntryDB;
import org.junit.Assert;

import java.util.List;
import java.util.regex.Pattern;

/**
 * Created by bvdmitri on 24.09.15.
 */

public class Tests {

    @org.junit.Test
    public void patternFilterTest() {
        System.out.println("Pattern filter test");
        try {
            DatabaseSearcher database = new DatabaseSearcher("mydb", "postgres", "postgres");
            database.open();
            String vPattern = "TRBV29-_";
            String jPattern = "TRBJ1-_";
            PatternFilter VFilter = new PatternFilter(EntryDB.Fields.V.getFieldName(), vPattern);
            PatternFilter JFilter = new PatternFilter(EntryDB.Fields.J.getFieldName(), jPattern);
            List<CdrEntrySetDB> setentries = database.findEntries(VFilter, JFilter);
            for (CdrEntrySetDB setentry : setentries) {
                List<EntryDB> entries = setentry.getCdrEntries();
                for (EntryDB entry : entries) {
                    Assert.assertTrue(Pattern.matches("^TRBV29-[a-zA-Z0-9]$", entry.getV()));
                    Assert.assertTrue(Pattern.matches("^TRBJ1-[a-zA-Z0-9]$", entry.getJ()));
                }
            }
            VFilter.setMatch(false);
            JFilter.setMatch(false);
            List<CdrEntrySetDB> otherEntries = database.findEntries(VFilter, JFilter);
            for (CdrEntrySetDB otherEntry : otherEntries) {
                List<EntryDB> entries = otherEntry.getCdrEntries();
                for (EntryDB entry : entries) {
                    Assert.assertFalse(Pattern.matches("^TRBV29-[a-zA-Z0-9]$", entry.getV()));
                    Assert.assertFalse(Pattern.matches("^TRBJ1-[a-zA-Z0-9]$", entry.getJ()));
                }
            }
            database.close();
            System.out.println("Passed");
        } catch (Exception e) {
            System.out.println();
            System.out.println("Failed");
        }
    }

    @org.junit.Test
    public void matchFilterTest() {
        try {
        System.out.println("Match filter test");
        DatabaseSearcher database = new DatabaseSearcher("mydb", "postgres", "postgres");
        database.open();
        String cdrString = "CSVGTGGTNEKLF";
        MatchFilter matchFilter = new MatchFilter(CdrEntrySetDB.Fields.CDR3.getFieldName(), cdrString, true);
        List<CdrEntrySetDB> setList = database.findSet(matchFilter);
        for (CdrEntrySetDB entrySetDB : setList) {
            Assert.assertEquals(entrySetDB.getCdr3(), cdrString);
        }
        matchFilter.setMatch(false);
        setList = database.findSet(matchFilter);
        for (CdrEntrySetDB entrySetDB : setList) {
            Assert.assertNotEquals(entrySetDB.getCdr3(), cdrString);
        }
        database.close();
        } catch (Exception e) {
            System.out.println();
            System.out.println("Failed");
        }
    }
}
