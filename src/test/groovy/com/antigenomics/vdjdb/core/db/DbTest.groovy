package com.antigenomics.vdjdb.core.db

import com.antigenomics.vdjdb.util.Util
import org.junit.Test


class DbTest {
    @Test
    public void loadTest() {
        println "Testing default database load procedure"

        def cdrDatabase = new CdrDatabase()

        println "MD5=" + cdrDatabase.md5sum
        println "#ENTRIES=" + cdrDatabase.entryCount
        println "#CDRS=" + cdrDatabase.cdrCount

        def storedMd5 = Util.resourceStreamReader("${CdrDatabase.DEFAULT_DB_NAME}.md5").readLines()[0]
        assert cdrDatabase.md5sum == storedMd5
    }
}
