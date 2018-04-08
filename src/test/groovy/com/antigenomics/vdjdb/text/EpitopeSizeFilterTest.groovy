package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import org.junit.Test

class EpitopeSizeFilterTest {
    @Test
    void testCount() {
        final VdjdbInstance vdjdbInst = new VdjdbInstance(Util.resourceAsStream("vdjdb_17.meta.txt"),
                Util.resourceAsStream("vdjdb_17.txt"))

        def counts = EpitopeSizeFilterUtil.generateEpitopeCounts(vdjdbInst)

        println counts.sort { -it.value }

        assert counts["GLCTLVAML"] == 804

        def counts2 = EpitopeSizeFilterUtil.generateEpitopeCounts(vdjdbInst, "HomoSapiens", "TRB")

        assert counts2["GLCTLVAML"] == 696
    }

    @Test
    void testFilter() {
        final VdjdbInstance vdjdbInst = new VdjdbInstance(Util.resourceAsStream("vdjdb_17.meta.txt"),
                Util.resourceAsStream("vdjdb_17.txt"))

        def filter = EpitopeSizeFilterUtil.createEpitopeSizeFilter(vdjdbInst, null, null, 956)

        def gil = "GILGFVFTL"

        assert filter.values == [gil] as Set<String>

        vdjdbInst.filter([filter]).dbInstance.rows.each { it[ClonotypeDatabase.EPITOPE_COL_DEFAULT].value == gil }
    }
}
