package com.antigenomics.vdjdb.db

import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.TextColumn
import org.junit.Test

import static com.antigenomics.vdjdb.Util.resourceAsStream

class DbTest {
    @Test
    public void createTest1() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        assert database.columns.size() == resourceAsStream("vdjdb_legacy.meta").readLines().size() - 1
    }

    @Test
    public void createTest2() {
        def database = new Database([new TextColumn("v.segm"), new TextColumn("j.segm"), new SequenceColumn("cdr3")])

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"))

        assert database.rows[0]["cdr3"].value == "CASSLTSGSPYNEQF"
        assert database.rows[0]["v.segm"].value == "TRBV27"
        assert database.rows[0]["j.segm"].value == "TRBJ2-1"
    }

    @Test
    public void createTest3() {
        def database = new Database([new TextColumn("v.segm"), new TextColumn("j.segm"), new SequenceColumn("cdr3")])

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"))

        def clone = Database.create(database.search([], []))
        
        assert clone.rows.size() == database.rows.size()
    }

    @Test
    public void fillTest1() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"))

        assert database.rows.size() == resourceAsStream("vdjdb_legacy.txt").readLines().size() - 1
    }

    @Test
    public void fillTest2() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(
                [["CASSLTSGSPYNEQF", "TRBV27", "TRBJ2-1", "TRB",
                  "human", "HLA-A*02", "CMV", "NLVPMVATV", ".", ".", "19017975"]])

        assert database.rows.size() == 1
        assert database.rows[0]["cdr3"].value == "CASSLTSGSPYNEQF"
    }

    @Test
    public void fillTest3() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"), [new ExactTextFilter("origin", "CMV", true)])

        assert ((TextColumn) database["origin"]).values.contains("EBV")
        assert !((TextColumn) database["origin"]).values.contains("CMV")
    }
}
