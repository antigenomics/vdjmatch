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

package com.antigenomics.vdjdb.stat

import com.antigenomics.vdjdb.AtomicDouble
import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.db.Row
import com.antigenomics.vdjdb.sequence.SequenceColumn
import com.antigenomics.vdjdb.text.TextColumn

import java.util.concurrent.atomic.AtomicInteger

/**
 * A base class for accumulation, reporting and accessing database search statistics 
 */
class SearchSummary {
    /**
     * Database that was queries 
     */
    final Database database
    protected final List<String> columnNames
    protected final Map<String, SummaryStatisticsCounter> labeledCounters = new HashMap<>()
    private final SummaryStatisticsCounter total = new SummaryStatisticsCounter(),
                                           foundOnce = new SummaryStatisticsCounter()

    /**
     * Creates an empty database search summary. Column name order only matters in case {@link #listCombinations}
     * is planned to be used (e.g. for constructing antigen phylogeny sunburst chart or tree). In case a user 
     * would like to simply infer some top categories with {@link #listTopCombinations}, only corresponding 
     * column names should be specified, in any order
     * @param database database that will be queried
     * @param columnNames database column names to be summarized, order specifies column hierarchy
     */
    SearchSummary(Database database, List<String> columnNames) {
        columnNames.each {
            def column = database[it]
            if (column == null)
                throw new RuntimeException("Error creating search summary. Column $it does not exist")
            if (column instanceof SequenceColumn)
                throw new RuntimeException("Error creating search summary. Sequence-valued columns are not allowed ($it)")
        }

        this.database = database
        this.columnNames = columnNames

        database.rows.each { row ->
            def key = columnNames.collect { row[it].value }.join("\t")
            labeledCounters.put(key, new SummaryStatisticsCounter())
        }
    }

    /**
     * Append search results
     * @param rows found rows
     * @param weight query weight
     */
    void append(Collection<Row> rows, double weight) {
        total.update(weight)
        def combinations = rows.collect { row ->
            if (row.parent != database) {
                throw new RuntimeException("Rows should come from the same database that was used during initialization")
            }
            columnNames.collect { row[it].value }.join("\t")
        }.unique()

        if (combinations.size() == 1)
            foundOnce.update(weight)

        combinations.each {
            labeledCounters[it].update(weight)
        }
    }

    /**
     * List all possible combinations of entry values that correspond to column identifier specified 
     * during class creation. Order and hierarchy is preserved. Unique combinations of all lengths are returned,
     * e.g. database [[A, B, C], [1, 2, 3], [1, 2, 4], [2, 3, 4]]
     * will result in [[1],[2],[1,2],[2,3],[1,2,3],[2,3,4],[1,2,4]] .
     * @return combinations of entry values for specified hierarchy
     */
    List<List<String>> listCombinations() {
        def keyVariants = new HashSet<>(
                labeledCounters.keySet().collect {
                    def splitLine = it.split("\t")
                    (0..<splitLine.length).collect {
                        splitLine[0..it].join("\t")
                    }
                }.flatten()
        )

        keyVariants.collect { String it ->
            it.split("\t") as List<String>
        }
    }

    /**
     * List all possible top-level combinations of entry values that correspond to column identifier specified 
     * during class creation. Order and hierarchy is preserved. Unique combinations are returned,
     * e.g. database [[A, B, C], [1, 2, 3], [1, 2, 4], [2, 3, 4]]
     * will result in [[1,2,3],[2,3,4],[1,2,4]] 
     * @return top-level combinations of entry values for specified hierarchy
     */
    List<List<String>> listTopCombinations() {
        listCombinations().findAll { it.size() == columnNames.size() }
    }

    private Counter searchInner(String combination) {
        def counter = new SummaryStatisticsCounter()
        labeledCounters.entrySet().findAll { it.key.startsWith(combination) }.each { counter.update(it.value) }
        counter
    }

    /**
     * Gets a counter for a given entry value combination. Use {@link #listCombinations} or 
     * {@link #listTopCombinations} to get the list of allowed combinations
     * @param values combination of entry values
     * @return search summary counter
     */
    Counter getCombinationCounter(List<String> values) {
        if (columnNames.size() < values.size()) {
            throw new RuntimeException("More values provided than the number of columns")
        }

        values.eachWithIndex { it, ind ->
            def column = database[columnNames[ind]]
            if (!(column as TextColumn).values.contains(it)) {
                throw new RuntimeException("Bad value list, '$it' not found among selected columns. " +
                        "Note that list ordering matters here.")
            }
        }

        def key = values.join("\t")

        if (columnNames.size() == values.size()) {
            return labeledCounters[key]
        } else {
            return searchInner(key)
        }
    }

    /**
     * Get the counter for queries found in the database 
     * @return counter for queries found in the database
     */
    Counter getFoundCounter() {
        total
    }

    /**
     * Get the counter for queries found once in the database 
     * @return counter for queries found once in the database
     */
    Counter getFoundOnce() {
        foundOnce
    }

    /**
     * Get the counter for queries found two or more times in the database
     * @return counter for queries found two or more times in the database
     */
    Counter getFoundTwiceAndMore() {
        new CounterImpl(total.uniqueCount - foundOnce.uniqueCount,
                total.weightedCount - foundOnce.weightedCount)
    }

    private static class SummaryStatisticsCounter implements Counter {
        private final AtomicInteger uniqueCounter = new AtomicInteger()
        private final AtomicDouble weightedCounter = new AtomicDouble()

        void update(double weight) {
            uniqueCounter.incrementAndGet()
            weightedCounter.addAndGet(weight)
        }

        void update(SummaryStatisticsCounter statisticsCounter) {
            uniqueCounter.addAndGet(statisticsCounter.uniqueCount)
            weightedCounter.addAndGet(statisticsCounter.weightedCount)
        }

        @Override
        int getUniqueCount() {
            uniqueCounter.get()
        }

        @Override
        double getWeightedCount() {
            weightedCounter.get()
        }
    }
}
