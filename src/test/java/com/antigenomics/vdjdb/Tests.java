package com.antigenomics.vdjdb;

import com.antigenomics.vdjdb.filters.MatchFilter;
import com.antigenomics.vdjdb.filters.PatternFilter;
import com.antigenomics.vdjdb.models.CdrEntrySetDB;
import com.antigenomics.vdjdb.models.EntryDB;
import com.antigenomics.vdjdb.util.PGUtil;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;
import java.util.regex.Pattern;

/**
 * Created by bvdmitri on 24.09.15.
 */

public class Tests {

    @org.junit.Test
    public void patternFilterTest() {
        System.out.println("Pattern filter test");
        Database database = new Database("mydb", "postgres", "postgres");
        database.open();
        String vPattern = "TRBV29-_";
        String jPattern = "TRBJ1-_";
        PatternFilter VFilter = new PatternFilter(EntryDB.Fields.V.getFieldName(), vPattern);
        PatternFilter JFilter = new PatternFilter(EntryDB.Fields.J.getFieldName(), jPattern);
        List<EntryDB> entries = database.findEntries(VFilter, JFilter);
        for (EntryDB entry : entries) {
            Assert.assertTrue(Pattern.matches("^TRBV29-[a-zA-Z0-9]$", entry.getV()));
            Assert.assertTrue(Pattern.matches("^TRBJ1-[a-zA-Z0-9]$", entry.getJ()));
        }
        VFilter.setMatch(false);
        JFilter.setMatch(false);
        List<EntryDB> otherEntries = database.findEntries(VFilter, JFilter);
        for (EntryDB otherEntry : otherEntries) {
            Assert.assertFalse(Pattern.matches("^TRBV29-[a-zA-Z0-9]$", otherEntry.getV()));
            Assert.assertFalse(Pattern.matches("^TRBJ1-[a-zA-Z0-9]$", otherEntry.getJ()));
        }
        database.close();
        System.out.println("Passed");
    }

    @org.junit.Test
    public void matchFilterTest() {
        System.out.println("Match filter test");
        Database database = new Database("mydb", "postgres", "postgres");
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
        System.out.println("Passed");
    }

    @Test
    public void initializationTest() {
        System.out.println("Installing database test");
        PGUtil.installPGDataBase("postgres", "postgres");
        System.out.println("Passed");
    }
}
