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
import com.antigenomics.vdjtools.sample.Clonotype

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong

class ClonotypeCounter {
    private final AtomicInteger uniqueCounter
    private final AtomicLong readCounter
    private final AtomicDouble frequencyCounter
    private final AtomicDouble weightCounter
    private final Set<Clonotype> clonotypes = ConcurrentHashMap.newKeySet()
    private final int databaseUnique

    ClonotypeCounter(int databaseUnique = -1) {
        this(0, 0, 0, databaseUnique, 0)
    }

    ClonotypeCounter(int unique, long reads, double frequency, double databaseUnique, double weight) {
        this.uniqueCounter = new AtomicInteger(unique)
        this.readCounter = new AtomicLong(reads)
        this.frequencyCounter = new AtomicDouble(frequency)
        this.weightCounter = new AtomicDouble(weight)
        this.databaseUnique = databaseUnique
    }

    /**
     * Updates the counter with abundance of current clonotype
     * @param clonotype clonotype to append
     */
    void update(Clonotype clonotype, double weight = 1.0) {
        if (!clonotypes.contains(clonotype)) {
            clonotypes.add(clonotype)

            uniqueCounter.incrementAndGet()
            readCounter.addAndGet(clonotype.count)
            frequencyCounter.addAndGet(clonotype.freq)
            weightCounter.addAndGet(weight) // todo: perhaps we need to compute max weight..
        }
    }

    /**
     * Number of unique clonotypes in a given category
     * @return number of clonotypes
     */
    int getUnique() {
        uniqueCounter.get()
    }

    /**
     * Total frequency of clonotypes in a given category
     * @return frequency (from 0 to 1)
     */
    double getFrequency() {
        frequencyCounter.get()
    }

    /**
     * Total read count of clonotypes in a given category
     * @return number of reads
     */
    long getReads() {
        readCounter.get()
    }

    /**
     * Number of unique clonotypes in a given category _in the original database_
     * @return number of clonotypes
     */
    long getDatabaseUnique() {
        databaseUnique
    }

    /**
     * Total weight of DB hits
     * @return weight
     */
    double getWeight() {
        weightCounter.get()
    }

    Set<Clonotype> getClonotypes() {
        Collections.unmodifiableSet(clonotypes)
    }

    static final String HEADER = "unique\tfrequency\treads\tdb.unique\tweight"

    @Override
    String toString() {
        unique + "\t" + frequency + "\t" + reads + "\t" + databaseUnique + "\t" + weight
    }
}
