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
 * A database column 
 */
abstract class Column {
    /**
     * Column identifier 
     */
    final String name
    /**
     * Column metadata
     */
    final Map<String, String> metadata

    /**
     * Creates an empty database column
     * @param name column identifier
     * @param metadata column metadata map
     */
    Column(String name, Map<String, String> metadata) {
        this.name = name
        this.metadata = metadata
    }

    /**
     * Adds an entry to this column. Depending on implementation indexing can be performed on this step. 
     * @param entry database entry
     */
    abstract void append(Entry entry)

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Column column = (Column) o

        if (name != column.name) return false

        return true
    }

    @Override
    int hashCode() {
        return name.hashCode()
    }
}