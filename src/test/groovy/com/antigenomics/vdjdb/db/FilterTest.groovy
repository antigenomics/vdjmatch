/*
 * Copyright 2013-{year} Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjdb.db

import com.antigenomics.vdjdb.TestUtil
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.PatternTextFilter
import com.antigenomics.vdjdb.text.SubstringTextFilter
import com.milaboratory.core.tree.TreeSearchParameters
import org.junit.Test

class FilterTest {
    @Test
    public void filterTest1() {
        def database = TestUtil.loadLegacyDb()

        assert database.search([], []).size() == database.rows.size()

        assert database.search([new ExactTextFilter("antigen.seq", "NLVPMVATV", false)], []).every {
            it.row["origin"].value == "CMV"
        }

        assert database.search([new SubstringTextFilter("antigen.seq", "LVPMVATV", false)], []).every {
            it.row["origin"].value == "CMV"
        }

        assert database.search([new PatternTextFilter("antigen.seq", "NL.PMV.TV", false)], []).every {
            it.row["origin"].value == "CMV"
        }
    }

    @Test
    public void filterTest2() {
        def database = TestUtil.loadLegacyDb()

        assert database.search([new ExactTextFilter("antigen.seq", "NLVPMVATV", false),
                                new ExactTextFilter("origin", "CMV", false)], []).size() > 0

        assert database.search([new ExactTextFilter("antigen.seq", "NLVPMVATV", false),
                                new ExactTextFilter("origin", "EBV", false)], []).size() == 0
    }

    @Test
    public void filterTest3() {
        def database = TestUtil.loadLegacyDb()

        assert new HashSet<>(database.search([new ExactTextFilter("cdr3", "CASSLAPGATNEKLFF", false)], [])*.row)
                .containsAll(database.search([], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF", new TreeSearchParameters(0, 0, 0))])*.row)

        assert database.search([new ExactTextFilter("antigen.seq", "NLVPMVATV", false),
                                new ExactTextFilter("origin", "CMV", false)], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF")]).size() > 0

        assert database.search([], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF")]).size() >
                database.search([], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF", new TreeSearchParameters(0, 0, 0))]).size()
    }

    @Test
    public void filterTest4() {
        def database = TestUtil.loadLegacyDb()

        assert new HashSet<>(database.search([new ExactTextFilter("cdr3", "CASSLAPGATNEKLFF", false)], [])*.row)
                .containsAll(database.search([], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF", new TreeSearchParameters(0, 0, 0))])*.row)

        assert database.search([], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF"),
                                    new SequenceFilter("antigen.seq", "NLVPMVATV")]).size() > 0
        
        assert database.search([], [new SequenceFilter("cdr3", "CASSLAPGATNEKLFF"),
                                    new SequenceFilter("antigen.seq", "TPRVTGGGAM")]).size() == 0
    }
}
