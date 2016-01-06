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

package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.fixer.Cdr3Fixer
import com.antigenomics.vdjdb.fixer.FixerResult

def DEFAULT_SPECIES_COL = "species", DEFAULT_V_COL = "v.segm", DEFAULT_J_COL = "j.segm",
    DEFAULT_CDR3_COL = "cdr3", DEFAULT_MIN_MATCH = "2", DEFAULT_MAX_MISMATCH = "1"

def cli = new CliBuilder(usage: "FixCdr3 [options] " +
        "input_file output_file")
cli.s(longOpt: "species-col", argName: "name", args: 1,
        "Species column name. [default = $DEFAULT_SPECIES_COL]")
cli.v(longOpt: "v-col", argName: "name", args: 1,
        "Variable segment column name. [default = $DEFAULT_V_COL]")
cli.j(longOpt: "j-col", argName: "name", args: 1,
        "Joining segment column name. [default = $DEFAULT_J_COL]")
cli.c(longOpt: "cdr3-col", argName: "name", args: 1,
        "CDR3 column name. [default = $DEFAULT_CDR3_COL]")
cli._(longOpt: "min-match-size", argName: "int", args: 1,
        "Minimal match between CDR3 and V/J germline. [default = $DEFAULT_MIN_MATCH]")
cli._(longOpt: "max-mismatch-size", argName: "int", args: 1,
        "Maximal number of mismatches between CDR3 and V/J germline. [default = $DEFAULT_MAX_MISMATCH]")

def opt = cli.parse(args)


if (opt == null) {
    System.exit(1)
}

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(1)
}

def inputFile = opt.arguments()[0], outputFile = opt.arguments()[1]
def speciesCol = (opt.s ?: DEFAULT_SPECIES_COL),
    vCol = (opt.v ?: DEFAULT_V_COL), jCol = (opt.j ?: DEFAULT_J_COL),
    cdr3Col = (opt.c ?: DEFAULT_CDR3_COL),
    minMatch = (opt.'min-match-size' ?: DEFAULT_MIN_MATCH).toInteger(),
    maxMismatch = (opt.'max-mismatch-size' ?: DEFAULT_MAX_MISMATCH).toInteger()

boolean first = true
def vInd = -1, jInd = -1, speciesInd = -1, cdr3Ind = -1
def fixer = new Cdr3Fixer(maxMismatch, minMatch)

new File(outputFile).withPrintWriter { pw ->
    new File(inputFile).splitEachLine('\t') { splitLine ->
        if (first) {
            first = false
            pw.println(splitLine.join("\t") + "\t" + FixerResult.HEADER)
            vInd = splitLine.findIndexOf { it.toLowerCase() == vCol.toLowerCase() }
            jInd = splitLine.findIndexOf { it.toLowerCase() == jCol.toLowerCase() }
            speciesInd = splitLine.findIndexOf { it.toLowerCase() == speciesCol.toLowerCase() }
            cdr3Ind = splitLine.findIndexOf { it.toLowerCase() == cdr3Col.toLowerCase() }

            if ([vInd, jInd, speciesInd, cdr3Ind].any { it < 0 }) {
                println "[ERROR] One or more required columns were not found."
                System.exit(0)
            }

            return
        }

        def fixResult = fixer.fix(splitLine[cdr3Ind], splitLine[vInd], splitLine[jInd], splitLine[speciesInd])

        pw.println(splitLine.join("\t") + "\t" + fixResult.toString())
    }
}