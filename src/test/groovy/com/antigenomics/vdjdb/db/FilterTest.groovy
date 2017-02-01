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

import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.PatternTextFilter
import com.antigenomics.vdjdb.text.SubstringTextFilter
import com.milaboratory.core.tree.TreeSearchParameters
import org.junit.Test

import static com.antigenomics.vdjdb.TestUtil.*
import static com.antigenomics.vdjdb.impl.ClonotypeDatabase.CDR3_COL_DEFAULT

class FilterTest {
    @Test
    public void filterTest1() {
        def database = loadLegacyDb()

        assert database.search([], []).size() == database.rows.size()

        assert !database.search([new ExactTextFilter(PEPTIDE_COL, "NLVPMVATV", false)], []).empty
        assert !database.search([new SubstringTextFilter(PEPTIDE_COL, "LVPMVATV", false)], []).empty
        assert !database.search([new PatternTextFilter(PEPTIDE_COL, "NL.PMV.TV", false)], []).empty

        assert database.search([new ExactTextFilter(PEPTIDE_COL, "NLVPMVATV", false)], []).every {
            it.row[SOURCE_COL].value == "CMV"
        }

        assert database.search([new SubstringTextFilter(PEPTIDE_COL, "LVPMVATV", false)], []).every {
            it.row[SOURCE_COL].value == "CMV"
        }

        assert database.search([new PatternTextFilter(PEPTIDE_COL, "NL.PMV.TV", false)], []).every {
            it.row[SOURCE_COL].value == "CMV"
        }
    }

    @Test
    public void filterTest2() {
        def database = loadLegacyDb()

        assert database.search([new ExactTextFilter(PEPTIDE_COL, "NLVPMVATV", false),
                                new ExactTextFilter(SOURCE_COL, "CMV", false)], []).size() > 0

        assert database.search([new ExactTextFilter(PEPTIDE_COL, "NLVPMVATV", false),
                                new ExactTextFilter(SOURCE_COL, "EBV", false)], []).size() == 0
    }

    @Test
    public void filterTest3() {
        def database = loadLegacyDb()

        assert new HashSet<>(database.search([new ExactTextFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF", false)], [])*.row)
                .containsAll(database.search([], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF")])*.row)

        assert database.search([new ExactTextFilter(PEPTIDE_COL, "NLVPMVATV", false),
                                new ExactTextFilter(SOURCE_COL, "CMV", false)], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF")]).size() > 0

        assert database.search([], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF", new TreeSearchParameters(2, 1, 1))]).size() >
                database.search([], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF")]).size()
    }

    @Test
    public void filterTest4() {
        def database = loadLegacyDb()

        assert new HashSet<>(database.search([new ExactTextFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF", false)], [])*.row)
                .containsAll(database.search([], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF", new TreeSearchParameters(0, 0, 0))])*.row)

        assert database.search([], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF"),
                                    new SequenceFilter(PEPTIDE_COL, "NLVPMVATV")]).size() > 0

        assert database.search([], [new SequenceFilter(CDR3_COL_DEFAULT, "CASSLAPGATNEKLFF"),
                                    new SequenceFilter(PEPTIDE_COL, "TPRVTGGGAM")]).size() == 0
    }
}
