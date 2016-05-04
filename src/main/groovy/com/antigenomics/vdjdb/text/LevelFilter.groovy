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
 * A filter for plain-text columns containing numeric values. The filter passes if the value is greater than the
 * threshold. Non-numeric values fail filter by default.
 */
class LevelFilter extends TextFilter {
    private double value

    /**
     * Creates a new numeric value filter. Non-numeric values pass by default.
     * @param columnId column identifier
     * @param value a lower threshold for a value to pass the filter
     * @param negative inverse the filter
     */
    LevelFilter(String columnId, String value, boolean negative) {
        super(columnId, value, negative)
        if (!value.isDouble()) {
            throw new RuntimeException("Numeric value should be provided to level filter ($value)")
        }
        this.value = value.toDouble()
    }

    @Override
    protected boolean passInner(Entry entry) {
        entry.value.isDouble() ? value <= entry.value.toDouble() : false
    }
}
