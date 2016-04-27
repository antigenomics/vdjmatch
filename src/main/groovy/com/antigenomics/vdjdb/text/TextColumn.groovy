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

package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.Entry
import com.fasterxml.jackson.annotation.JsonIgnore

/**
 * A column holding plain-text values, this column is hashed, filters are applied on-the-fly 
 */
class TextColumn extends Column {
    private final Map<String, List<Entry>> map = new HashMap<>()

    /**
     * Creates an empty plain-text column 
     * @param name column name
     * @param metadata column metadata
     */
    TextColumn(String name, Map<String, String> metadata = [:]) {
        super(name, metadata)
    }

    /**
     * Gets the set of all possible values in the column 
     * @return a set of unique values in the column
     */
    @JsonIgnore
    Set<String> getValues() {
        map.keySet()
    }

    /**
     * Adds an entry to the column. Entry value will be hashed
     * @param entry entry to add
     */
    @Override
    void append(Entry entry) {
        if (entry.value.length() > 0) {
            def entries = map[entry.value]
            if (entries == null) {
                map.put(entry.value, entries = new ArrayList<Entry>())
            }
            entries.add(entry)
        }
    }
}
