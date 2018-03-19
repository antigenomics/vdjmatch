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

package com.antigenomics.vdjdb.db

import groovy.transform.CompileStatic
import com.fasterxml.jackson.annotation.*

/**
 * A database row 
 */
@CompileStatic
class Row {
    /**
     * Parent database 
     */
    @JsonIgnore
    final Database parent
    /**
     * Row index 
     */
    final int index
    /**
     * Entries in the row 
     */
    final Entry[] entries

    /**
     * Creates a new row 
     * @param parent parent database
     */
    Row(Database parent) {
        this.parent = parent
        this.index = parent.rows.size()
        this.entries = new Entry[parent.columns.size()]
    }

    /**
     * Gets an entry by index 
     * @param index entry index
     * @return entry
     */
    Entry getAt(int index) {
        entries[index]
    }

    /**
     * Gets an entry by column indentifier
     * @param name column identifier
     * @return entry
     */
    Entry getAt(String name) {
        entries[parent.getColumnIndex(name)]
    }

    @Override
    boolean equals(o) {
        index == ((Row) o).index
    }

    @Override
    int hashCode() {
        index
    }

    String toTabDelimitedString() {
        entries.collect { Entry it -> it.value }.join("\t")
    }

    @Override
    String toString() {
        entries.collect { Entry it -> it.value }.toString()
    }
}