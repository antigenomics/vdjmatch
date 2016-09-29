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
class FieldLevelFilter extends TextFilter {
    private double value
    private String field

    FieldLevelFilter(String columnId, String field, String value) {
        super(columnId, value, false)
        if (!value.isDouble()) {
            throw new RuntimeException("Numeric value should be provided to level filter ($value)")
        }
        this.value = value.toDouble()
        this.field = field
    }

    @Override
    protected boolean passInner(Entry entry) {
        def entryFields = getEntryFields(entry.value)

        for (String entryField : entryFields) {
            if (entryField.contains(field)) {
                def split = entryField.split(":")
                def number = split[1].replace("\"", "")
                if (!number.isEmpty()) {
                    def entryValue = parseNumber(number)
                    return entryValue >= value
                } else {
                    return false
                }
            }
        }
        return false
    }

    private static Collection<String> getEntryFields(String value) {
        value.toLowerCase().split(",").collect { it.toLowerCase() }
    }

    private static double parseNumber(String number) {
        if (number.contains("/")) {
            def rat = number.split("/");
            if (rat[0].isEmpty() || rat[1].isEmpty()) {
                return 0.0
            }
            return rat[0].toDouble() / rat[1].toDouble();
        } else {
            return number.toDouble()
        }
    }
}
