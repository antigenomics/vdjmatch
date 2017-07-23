package com.antigenomics.vdjdb.web

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.VdjdbInstance
import org.junit.Test

/**
 * Created by mikesh on 7/23/17.
 */
class EpitopeSuggestionTest {
    @Test
    void test() {
        final VdjdbInstance vdjdbInst = new VdjdbInstance(Util.resourceAsStream("vdjdb_17.meta.txt"),
                Util.resourceAsStream("vdjdb_17.txt"))

        def z = EpitopeSuggestionGenerator.generateSuggestions(vdjdbInst)

        assert EpitopeSuggestionGenerator.generateSuggestions(vdjdbInst)["EAAGIGILTV"].size() > 0
    }
}
