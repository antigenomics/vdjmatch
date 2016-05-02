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
 * *INTERNAL* A filter that evaluates an expression for each row
 */
class ExpressionFilterBatch implements FilterBatch {
    final static String FILTER_MARK = "__"

    final String filter

    ExpressionFilterBatch(Database database, String filter) {
        filter.split(FILTER_MARK).each { token ->
            if (database.hasColumn(token)) {
                def columnIndex = database.getColumnIndex(token)

                filter = filter.replaceAll("$FILTER_MARK$token$FILTER_MARK", "x[$columnIndex].value")
            }
        }

        if (filter.contains(FILTER_MARK)) {
            throw new RuntimeException("Failed to parse filter, " +
                    "perhaps some of the columns do not exist in the database")
        }

        this.filter = filter
    }

    @Override
    boolean pass(Row row) {
        (boolean) Eval.x(row, filter)
    }
}
