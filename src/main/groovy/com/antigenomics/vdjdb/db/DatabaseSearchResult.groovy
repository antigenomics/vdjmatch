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

import com.antigenomics.vdjdb.sequence.Hit
import groovy.transform.CompileStatic

/**
 * A database search result holding database row and
 * corresponding alignments/scores if search was performed using 'sequence' columns
 */
@CompileStatic
class DatabaseSearchResult implements SearchResult {
    final Row row
    final Hit[] hits

    DatabaseSearchResult(Row row, Hit[] hits) {
        this.row = row
        this.hits = hits
    }
}
