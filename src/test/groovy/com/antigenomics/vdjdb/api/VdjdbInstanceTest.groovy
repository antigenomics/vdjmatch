/*
 * Copyright 2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjdb.api

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import org.junit.Test


class VdjdbInstanceTest {
    final VdjdbInstance vdjdbInst = new VdjdbInstance(Util.resourceAsStream("vdjdb_new.meta.txt"),
            Util.resourceAsStream("vdjdb_new.txt"))

    @Test
    void headerTest() {
        println vdjdbInst.header.collect { it.name }
        assert ['cdr3', 'antigen'].every { name -> vdjdbInst.header.any { Column col -> col.name == name } }
    }

    @Test
    void filterTest() {
        assert !vdjdbInst.filter().columns.empty
        assert !vdjdbInst.filter().rows.empty
    }

    @Test
    void clonotypeDbTest() {
        def row = vdjdbInst.dbInstance.rows[0]
        assert !vdjdbInst.asClonotypeDatabase().search(
                row[ClonotypeDatabase.V_COL_DEFAULT].value,
                row[ClonotypeDatabase.J_COL_DEFAULT].value,
                row[ClonotypeDatabase.CDR3_COL_DEFAULT].value
        ).empty
    }
}