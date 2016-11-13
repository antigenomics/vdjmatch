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
class JSONSubstringFilter extends TextFilter {
    private Collection<Collection<String>> values

    JSONSubstringFilter(String columnId, String value) {
        super(columnId, value, false)
        this.values = parse_value(value)
    }

    @Override
    protected boolean passInner(Entry entry) {
        def jsonFields = parse_json(entry.value)

        for (Collection<String> c: values) {
            def nameOnly = false
            if (c.size() < 2 || c.getAt(1).isEmpty()) nameOnly = true

            def found = false

            for (Collection<String> field: jsonFields) {
                def fieldName = field.getAt(0)
                def fieldValue = field.getAt(1)

                if (nameOnly) {
                    if (fieldName.contains(c.getAt(0)) && fieldValue != null && !fieldValue.isEmpty()) {
                        found = true
                    }
                } else {
                    if (fieldName.contains(c.getAt(0)) && fieldValue != null && !fieldValue.isEmpty() && fieldValue.contains(c.getAt(1))) {
                        found = true
                    }
                }
            }

            if (!found) {
                return false
            }

        }
        return true
    }

    private static Collection<Collection<String>> parse_json(String json) {
        json.split("\",\"").collect { it.split("\":\"").collect { it.trim().replaceAll(/^[{,},\"]*$/, "") } }
    }

    private static Collection<Collection<String>> parse_value(String value) {
        value.split(",").collect { it.split(":").collect { it.trim() } }
    }
}
