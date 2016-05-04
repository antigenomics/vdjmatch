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
 * An entry filtering rule based on exact text matching. Works on values containing comma-separated items.
 * Only rows that have at least one item that match an item from the set supplied to filter using "value"
 * parameter are retained. If either the set searched for or the set supplied to filter is empty (or contain an
 * item with a value of ".") the filter passes automatically.
 */
class ExactSetTextFilter extends TextFilter {
    private final Set<String> splitValue

    /**
     * Creates a new entry filtering rule
     * @param columnId column identifier
     * @param value value to be matched
     * @param negative invert filter
     */
    ExactSetTextFilter(String columnId, String value, boolean negative) {
        super(columnId, value, negative)
        splitValue = value.toLowerCase().split(",") as Set<String>
    }

    @Override
    protected boolean passInner(Entry entry) {
        if (checkAutoPass(splitValue))
            return true
        def otherSplitValue = entry.value.split(",") as Collection<String>
        if (checkAutoPass(otherSplitValue))
            return true
        otherSplitValue.any { splitValue.contains(it.toLowerCase()) }
    }

    private static boolean checkAutoPass(Collection<String> set) {
        set.empty || set.contains(".")
    }
}
