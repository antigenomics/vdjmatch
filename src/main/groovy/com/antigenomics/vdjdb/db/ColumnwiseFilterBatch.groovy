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

import com.antigenomics.vdjdb.text.TextFilter

class ColumnwiseFilterBatch implements FilterBatch {
    final int[] filterColIds
    final List<TextFilter> filters

    ColumnwiseFilterBatch(final Database database, final List<TextFilter> filters) {
        this.filterColIds = database.getFilterColIds(filters)
        this.filters = filters
    }

    @Override
    boolean pass(Row row) {
        for (int i = 0; i < filters.size(); i++) {
            if (!filters[i].pass((Entry) row[filterColIds[i]])) {
                return false
            }
        }
        true
    }
}
