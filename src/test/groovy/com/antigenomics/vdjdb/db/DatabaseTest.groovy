package com.antigenomics.vdjdb.db

import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.TextColumn
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.*
import static com.antigenomics.vdjdb.Util.resourceAsStream
import static com.antigenomics.vdjdb.impl.ClonotypeDatabase.*

class DatabaseTest {
    @Test
    public void createTest1() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        assert database.columns.size() == resourceAsStream("vdjdb_legacy.meta").readLines().size() - 1
    }

    @Test
    public void createTest2() {
        def database = new Database([new TextColumn(ID_COL),
                                     new TextColumn(V_COL_DEFAULT), new TextColumn(J_COL_DEFAULT),
                                     new SequenceColumn(CDR3_COL_DEFAULT)])

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"))

        def row = database.rows.find { it[ID_COL].value == "VDJDB000.1" }

        assert row != null
        assert row[CDR3_COL_DEFAULT].value == "CASSLGQAYEQYF"
        assert row[V_COL_DEFAULT].value == "TRBV7-8"
        assert row[J_COL_DEFAULT].value == "TRBJ2-7"
    }

    @Test
    public void createTest3() {
        def database = new Database([new TextColumn(V_COL_DEFAULT), new TextColumn(J_COL_DEFAULT),
                                     new SequenceColumn(CDR3_COL_DEFAULT)])

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
                [["VDJDB000.3",
                  "CASSLAPGATNEKLFF",
                  "TRBV1",
                  "TRBJ1-4",
                  "TRB",
                  "human",
                  "HLA-A*0201",
                  "non-self",
                  "infection",
                  "viral",
                  "CMV",
                  "NLVPMVATV",
                  "CMVpp65 495-503",
                  "restimulation",
                  "public",
                  "known",
                  "21555537"]])

        assert database.rows.size() == 1
        assert database.rows[0][CDR3_COL_DEFAULT].value == "CASSLAPGATNEKLFF"
    }

    @Test
    public void fillTest3() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"), [new ExactTextFilter(SOURCE_COL, "CMV", true)])

        assert ((TextColumn) database[SOURCE_COL]).values.contains("EBV")
        assert !((TextColumn) database[SOURCE_COL]).values.contains("CMV")
    }

    @Test
    public void fillTest4() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"), "__${SOURCE_COL}__=~/(EBV|influenza)/")

        assert ((TextColumn) database[SOURCE_COL]).values.contains("EBV")
        assert !((TextColumn) database[SOURCE_COL]).values.contains("CMV")
        assert ((TextColumn) database[SOURCE_COL]).values.contains("influenza")
    }

    @Test
    public void fillTest5() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"),
                "__${SOURCE_COL}__==\"EBV\" && __${PEPTIDE_COL}__==\"NLVPMVATV\"")

        assert database.rows.empty

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"),
                "__${SOURCE_COL}__==\"CMV\" && __${PEPTIDE_COL}__!=\"NLVPMVATV\"")

        assert !database.rows.empty
    }
}
