package com.antigenomics.vdjdb.sequence

import org.junit.Test

class SearchScopeTest {
    @Test
    void test1() {
        new SearchScope(3,2,3).treeSearchParameters
    }

    @Test
    void test2() {
        new SearchScope(3,1,1,3).treeSearchParameters
    }
}
