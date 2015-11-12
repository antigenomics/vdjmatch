package com.antigenomics.vdjdb.core.db

import com.antigenomics.vdjdb.core.Util
import org.junit.Test

class DbTest {
    @Test
    public void loadTest() {
        println "Testing default database load procedure"

        def cdrDatabase = new CdrDatabase()

        //println "MD5=" + cdrDatabase.md5sum
        println "#ENTRIES=" + cdrDatabase.entryCount
        println "#CDRS=" + cdrDatabase.cdrCount

        //def storedMd5 = Util.resourceStreamReader("${CdrDatabase.DEFAULT_DB_NAME}.md5").readLines()[0]
        //assert cdrDatabase.md5sum == storedMd5

        assert cdrDatabase.entryCount == (Util.resourceStreamReader("${CdrDatabase.DEFAULT_DB_NAME}.txt").readLines().size() - 1)
        
        new File("src/main/resources/db/vdjdb.md5").write(cdrDatabase.md5sum)
    }
}
