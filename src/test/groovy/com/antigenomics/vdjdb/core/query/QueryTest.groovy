package com.antigenomics.vdjdb.core.query

import com.antigenomics.vdjdb.core.db.CdrDatabase
import org.junit.Test


class QueryTest {
    @Test
    public void infamousLapgatTest() {
        println "Testing database search"

        def seq = "CASSLAPGATNEKLFF"

        def db = new CdrDatabase()

        def searcher = new CdrDatabaseSearcher(db)

        searcher.search(seq).each {
            println it
        }

        println "Testing exact match"
        assert searcher.lucky(seq).cdrEntrySet == db[seq]
        println "fine"
    }
}
