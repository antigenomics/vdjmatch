/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 *
 * Last modified on 17.11.2014 by mikesh
 */


package com.antigenomics.vdjdb.core.db

import com.antigenomics.vdjdb.util.Util
import groovy.transform.PackageScope

import java.security.MessageDigest

class CdrDatabase implements Iterable<CdrEntrySet> {
    public static final String FILTER_MARK = "__",
                               DEFAULT_DB_NAME = "db/vdjdb"

    private final BigInteger md5sum
    private final String fileName

    private final Map<String, CdrEntrySet> entriesByCdr = new HashMap<>()
    private final Map<String, Integer> field2Index = new HashMap<>()
    public final String[] annotationHeader
    public final String ANNOTATION_HEADER
    private int entryCount = 0

    public CdrDatabase() {
        this(null)
    }

    public CdrDatabase(String filter) {
        this(Util.resourceStreamReader("${DEFAULT_DB_NAME}.txt"), "${DEFAULT_DB_NAME}.txt", filter)
    }

    public CdrDatabase(String fileName, String filter) {
        this(new InputStreamReader(new FileInputStream(fileName)), fileName, filter)
    }

    public CdrDatabase(InputStreamReader dbReader, String fileName, String filter) {

        def headerLine = dbReader.readLine()
        if (!headerLine.startsWith("#"))
            throw new Exception("Header line SHOULD start with '#'")
        headerLine = headerLine[1..-1]
        def header = headerLine.split("\t")

        if (filter) {
            filter.split(FILTER_MARK).each { token ->
                def columnIndex = header.findIndexOf { it.toUpperCase() == token.toUpperCase() }
                if (columnIndex >= 0)
                    filter = filter.replaceAll("$FILTER_MARK$token$FILTER_MARK", "x[$columnIndex]")
            }

            if (filter.contains(FILTER_MARK))
                throw new Exception("Failed to parse filter, " +
                        "perhaps some of the columns do not exist in the database")
        }

        def cdr3aaInd = -1, vInd = -1, jInd = -1
        def annotationHeader = new ArrayList<String>()
        def annotationIndices = new ArrayList<Integer>()

        header.eachWithIndex { String it, int ind ->
            switch (it.toUpperCase()) {
                case "CDR3AA":
                    cdr3aaInd = ind
                    break
                case "V":
                    vInd = ind
                    break
                case "J":
                    jInd = ind
                    break
                default:
                    annotationHeader.add(it)
                    annotationIndices.add(ind)
                    break
            }
        }

        if (cdr3aaInd < 0 || vInd < 0 || jInd < 0)
            throw new Exception("The following columns are MANDATORY: cdr3aa, v and j columns")

        this.fileName = fileName
        this.annotationHeader = annotationHeader as String[]
        this.ANNOTATION_HEADER = annotationHeader.join("\t")
        annotationHeader.eachWithIndex { it, ind -> field2Index.put(it, ind) }

        MessageDigest digest = MessageDigest.getInstance("MD5")

        def line
        while ((line = dbReader.readLine()) != null) {
            digest.update(line.bytes);
            def splitLine = line.split("\t") as String[]

            if (filter && !Eval.x(splitLine, filter))
                continue

            String cdr3aa = splitLine[cdr3aaInd], v, j
            (v, j) = Util.extractVDJ(splitLine[[vInd, jInd]])

            def entrySet = entriesByCdr[cdr3aa]
            if (!entrySet)
                entriesByCdr.put(cdr3aa, entrySet = new CdrEntrySet(this, cdr3aa))
            entrySet.addEntry(v, j, splitLine[annotationIndices] as String[])

            entryCount++
        }

        this.md5sum = new BigInteger(1, digest.digest())
    }

    String getMd5sum() {
        md5sum.toString(16).padLeft(32, '0')
    }

    @PackageScope
    int getIndex(String field) {
        field2Index[field] ?: -1
    }

    public CdrEntrySet getAt(String cdr3aa) {
        entriesByCdr[cdr3aa]
    }

    @Override
    public Iterator iterator() {
        entriesByCdr.values().iterator()
    }

    public int getCdrCount() {
        entriesByCdr.size()
    }

    public int getEntryCount() {
        entryCount
    }

    String getFileName() {
        fileName
    }

    @Override
    public boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        CdrDatabase that = (CdrDatabase) o

        ANNOTATION_HEADER == that.ANNOTATION_HEADER && fileName == that.fileName
    }

    @Override
    public int hashCode() {
        31 * fileName.hashCode() +
                ANNOTATION_HEADER.hashCode()
    }
}
