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

import com.antigenomics.vdjdb.db.Entry


/**
 * Internal
 */
class FieldSubstringFilter extends TextFilter {
    private String value
    private String field

    FieldSubstringFilter(String columnId, String field, String value) {
        super(columnId, value, false)
        this.value = value.toLowerCase()
        this.field = field
    }

    @Override
    protected boolean passInner(Entry entry) {
        def entryFields = getEntryFields(entry.value)

        for (String entryField : entryFields) {
            if (entryField.contains(field) && entryField.contains(value)) {
                return true
            }
        }
        return false
    }

    private static Collection<String> getEntryFields(String value) {
        value.toLowerCase().split(",").collect { it.toLowerCase() }
    }
}
