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

package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjtools.io.InputStreamFactory
import com.antigenomics.vdjtools.io.SampleStreamConnection
import com.antigenomics.vdjtools.sample.Sample

import java.util.zip.GZIPInputStream

import static com.antigenomics.vdjdb.Util.resourceAsStream

class TestUtil {
    static final String ID_COL = "id", SOURCE_COL = "source", PEPTIDE_COL = "antigen.seq"

    static Database loadLegacyDb() {
        def database = new Database(resourceAsStream("vdjdb_legacy.meta.txt"))

        database.addEntries(resourceAsStream("vdjdb_legacy.txt"))

        database
    }

    static final Sample TEST_SAMPLE = SampleStreamConnection.load([
            create: {
                new GZIPInputStream(resourceAsStream("sergey_anatolyevich.gz"))
            },
            getId : { "sergey_anatolyevich.gz" }
    ] as InputStreamFactory)
}
