/*
 * Copyright 2013-{year} Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
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

class SummaryStatistics {
    final Database database
    final List<String> columnNames
    protected final Map<String, SummaryStatisticsCounter> labeledCounters = new HashMap<>()
    private final SummaryStatisticsCounter total = new SummaryStatisticsCounter(),
                                           foundOnce = new SummaryStatisticsCounter()

    SummaryStatistics(Database database, List<String> columnNames) {
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

    List<List<String>> listCombinations() {
        def keyVariants = new HashSet<>(labeledCounters.keySet().collect {
            def splitLine = it.split("\t")
            (0..<splitLine.length).collect {
                splitLine[0..it].join("\t")
            }
        }.flatten())

        keyVariants.collect { it.split("\t") }
    }

    private Counter searchInner(String combination) {
        def counter = new SummaryStatisticsCounter()
        labeledCounters.entrySet().findAll { it.key.startsWith(combination) }.each { counter.update(it.value) }
        counter
    }

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

    Counter getFoundCounter() {
        total
    }

    Counter getFoundOnce() {
        foundOnce
    }

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
