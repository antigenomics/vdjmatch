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

import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong

class ClonotypeCounter {
    private final AtomicInteger uniqueCounter
    private final AtomicLong readCounter
    private final AtomicDouble frequencyCounter

    ClonotypeCounter() {
        this(0, 0, 0)
    }

    ClonotypeCounter(int unique, long reads, double frequency) {
        this.uniqueCounter = new AtomicInteger(unique)
        this.readCounter = new AtomicLong(reads)
        this.frequencyCounter = new AtomicDouble(frequency)
    }

    /**
     * Updates the counter with abundance of current clonotype
     * @param clonotype clonotype to append
     */
    void update(Clonotype clonotype) {
        uniqueCounter.incrementAndGet()
        readCounter.addAndGet(clonotype.count)
        frequencyCounter.addAndGet(clonotype.freq)
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

    static final String HEADER = "unique\tfrequency\treads"

    @Override
    String toString() {
        unique + "\t" + frequency + "\t" + reads
    }
}
