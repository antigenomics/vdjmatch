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

/**
 * A database entry 
 */
class Entry {
    /**
     * Parent column 
     */
    final Column column
    /**
     * Parent row 
     */
    final Row row
    /**
     * Entry value 
     */
    final String value

    /**
     * Creates a new entry
     * @param column parent column
     * @param row parent row
     * @param value entry value
     */
    Entry(Column column, Row row, String value) {
        this.column = column
        this.row = row
        this.value = value
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Entry entry = (Entry) o

        if (row != entry.row) return false
        if (value != entry.value) return false

        return true
    }

    @Override
    int hashCode() {
        int result
        result = row.hashCode()
        result = 31 * result + value.hashCode()
        return result
    }

    @Override
    String toString() {
        "$row.index:$column.name:$value"
    }
}