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

import com.antigenomics.vdjdb.text.ExactTextFilter
import com.antigenomics.vdjdb.text.PatternTextFilter
import com.antigenomics.vdjdb.text.SubstringTextFilter
import org.junit.Test

import static com.antigenomics.vdjdb.Util.resourceAsStream

class FilterTest {
    @Test
    public void filterTest1() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"))

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

        assert database.search([new ExactTextFilter("antigen.seq", "NLVPMVATV", false)], []).size() ==
                database.search([new ExactTextFilter("origin", "CMV", false)], []).size()
    }
}
