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
import com.antigenomics.vdjdb.db.Filter


/**
 * A base class for filtering based on plain-text entries. Only rows that have entry that match the filter are retained.
 */
abstract class TextFilter implements Filter {
    /**
     * Identifier of the column this filter will be applied to
     */
    final String columnId
    /**
     * Value to be matched
     */
    final String value
    /**
     * True if the filter is negated, false otherwise 
     */
    final boolean negative

    /**
     * Creates a new entry filtering rule
     * @param columnId column identifier
     * @param value value to be matched
     * @param negative invert filter
     */
    TextFilter(String columnId, String value, boolean negative) {
        this.columnId = columnId
        this.value = value
        this.negative = negative
    }

    protected abstract boolean passInner(Entry entry)

    /**
     * Checks if an entry is passing the filter, accounting for {@link #negative}
     * @param entry entry to check
     * @return true if entry passes the filter, false otherwise
     */
    boolean pass(Entry entry) {
        negative ^ passInner(entry)
    }

    @Override
    boolean isSequenceFilter() {
        false
    }
}
