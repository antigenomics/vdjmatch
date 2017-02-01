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

import com.antigenomics.vdjdb.scoring.AlignmentScoringProvider
import com.antigenomics.vdjdb.scoring.SequenceSearcherPreset
import com.antigenomics.vdjdb.stat.ClonotypeCounter
import com.antigenomics.vdjdb.stat.ClonotypeSearchSummary
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.misc.ExecUtil
import com.milaboratory.core.tree.TreeSearchParameters

if (args.length > 0 && args[0].toLowerCase() == "update") {
    Util.checkDatabase(true)
    System.exit(0)
}

def DEFAULT_SEARCH_SCOPE = "0,0,0,0",
    DEFAULT_SEARCH_SCORING = AlignmentScoringProvider.DUMMY_SCORING_ID,
    DEFAULT_CONFIDENCE_THRESHOLD = "0",
    ALLOWED_SPECIES_ALIAS = ["human" : "homosapiens", "mouse": "musmusculus",
                             "monkey": "macacamulatta"],
    ALLOWED_GENES = ["TRA", "TRB"]

def cli = new CliBuilder(usage: "vdjdb [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix\n" +
        "Input samples should be provided in VDJtools format if --software is not specified/supported.")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli._(longOpt: "search-scope", argName: "s,i,d,t", args: 1,
        "Sets CDR3 sequence search parameters: " +
                "allowed number of substitutions (s), insertions (i), deletions (d) and total number of mutations. " +
                "[default=$DEFAULT_SEARCH_SCOPE]")
cli._(longOpt: "search-scoring", argName: "name", args: 1,
        "Sets CDR3 alignment scoring preset. " +
                "Either '${AlignmentScoringProvider.DUMMY_SCORING_ID}' that counts number the number of " +
                "matched bases, or 'v1' that uses a pre-computed scoring matrix." +
                "[default='$DEFAULT_SEARCH_SCORING']")
cli._(longOpt: "database", argName: "string", args: 1, "Path and prefix of an external database. " +
        "The prefix should point to a '.txt' file (database itself) and '.meta.txt' (database column metadata).")
cli._(longOpt: "use-fat-db", "Use a more redundant database version, with extra fields (meta, method, etc). " +
        "Fat database can contain several records for a " +
        "TCR:pMHC pair corresponding to different replicates/tissue sources/targets.")
cli._(longOpt: "summary-columns", argName: "col1,col2,...", args: 1,
        "Table columns for summarizing, see DB metadata for column names. " +
                "[default=${ClonotypeSearchSummary.FIELDS_PLAIN_TEXT.join(",")}]. ")
cli._(longOpt: "filter", argName: "logical expression(__field__,...)", args: 1,
        "Logical filter evaluated for database columns. Supports Regex, .contains(), .startsWith(), etc.")
cli.S(longOpt: "species", argName: "name", args: 1, required: true,
        "Species of input sample(s), allowed values: ${ALLOWED_SPECIES_ALIAS.keySet()}.")
cli._(longOpt: "software", argName: "string", args: 1,
        "Input RepSeq data format. Currently supported: ${Software.values().join(", ")}. " +
                "[default = ${Software.VDJtools}]")
cli.R(longOpt: "gene", argName: "name", args: 1, required: true,
        "Receptor gene of input sample(s), allowed values: $ALLOWED_GENES.")
cli._(longOpt: "vdjdb-conf", argName: "[0,3]", args: 1,
        "VDJdb confidence level threshold, from lowest to highest. [default=$DEFAULT_CONFIDENCE_THRESHOLD]")
cli.v(longOpt: "v-match", "Require V segment matching.")
cli.j(longOpt: "j-match", "Require J segment matching.")
cli.c("Compressed output")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(1)
}

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(1)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "[ERROR] Only output prefix should be provided in case of -m"
    else
        println "[ERROR] At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(1)
}

// Basic arguments

def dbPrefix = (String) (opt.'database' ?: null),
    compress = (boolean) opt.c,
    vMatch = (boolean) opt."v-match", jMatch = (boolean) opt."j-match",
    species = (String) opt.S, gene = (String) opt.R,
    q = (opt.'vdjdb-conf' ?: DEFAULT_CONFIDENCE_THRESHOLD).toInteger(),
    filterStr = opt.'filter',
    useFatDb = (boolean) opt.'use-fat-db',
    outputPrefix = opt.arguments()[-1],
    software = opt.'software' ? Software.byName(opt.'software') : Software.VDJtools,
    summaryColumns = (opt.'summary-columns' ?: ClonotypeSearchSummary.FIELDS_PLAIN_TEXT.join(",")).split(",") as List<String>

def allowedSpecies = [ALLOWED_SPECIES_ALIAS.keySet(), ALLOWED_SPECIES_ALIAS.values()].flatten()
if (!allowedSpecies.any { species.equalsIgnoreCase(it) }) {
    println "Wrong species name, use one of ${allowedSpecies} (case-insensitive)"
    System.exit(1)
}

species = ALLOWED_SPECIES_ALIAS[species.toLowerCase()] ?: species

if (!ALLOWED_GENES.any { gene.equalsIgnoreCase(it) }) {
    println "Wrong gene name, use one of $ALLOWED_GENES (case-insensitive)"
    System.exit(1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Search parameters

def searchScope = opt.'search-scope' ?: DEFAULT_SEARCH_SCOPE,
    searchScoring = opt.'search-scoring' ?: DEFAULT_SEARCH_SCORING

def ss = searchScope.split(",").collect { it.toInteger() }

parameterPreset = new SequenceSearcherPreset(
        new TreeSearchParameters(ss[0], ss[1], ss[2], ss[3]),
        AlignmentScoringProvider.loadScoring(searchScoring))

println "[${new Date()} $scriptName] Loading database..."

// Either load db from specified path, or use built-in database

def database

if (dbPrefix) {
    def metaStream = new FileInputStream("${dbPrefix}.meta.txt"),
        dataStream = new FileInputStream("${dbPrefix}.txt")
    database = new VdjdbInstance(metaStream, dataStream)
} else {
    database = new VdjdbInstance(useFatDb)
}

def missingSummaryCols = summaryColumns.findAll { !database.header*.name.contains(it) }
if (!missingSummaryCols.empty) {
    println "Columns $missingSummaryCols specified for summary generation are missing in the database."
    System.exit(1)
}

println "[${new Date()} $scriptName] Loaded database. \n${database.dbInstance}"

// Expression filtering if specified
if (filterStr) {
    println "[${new Date()} $scriptName] Filtering using $filterStr."
    database = database.filter(filterStr)
    println "[${new Date()} $scriptName] Done. \n${database.dbInstance}"
}

println "[${new Date()} $scriptName] Preparing clonotype database for $species $gene."

database = database.asClonotypeDatabase(vMatch, jMatch, parameterPreset, species, gene, q)

println "[${new Date()} $scriptName] Done. \n$database"

if (database.rows.empty) {
    println "No records present in filtered database"
    System.exit(1)
}

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)..."

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) to process."

//
// Main loop
//

println "[${new Date()} $scriptName] Annotating sample(s) & writing results."

def sw = new SampleWriter(compress)

new File(ExecUtil.formOutputPath(outputPrefix, "annot", "summary")).withPrintWriter { pwSummary ->
    pwSummary.println([MetadataTable.SAMPLE_ID_COLUMN,
                       sampleCollection.metadataTable.columnHeader,
                       "db.column.name",
                       "db.column.value",
                       ClonotypeCounter.HEADER].flatten().join("\t"))

    sampleCollection.eachWithIndex { Sample sample, int ind ->
        println "[${new Date()} $scriptName] Annotating..."

        def sampleId = sample.sampleMetadata.sampleId

        def results = database.search(sample)

        def writer = sw.getWriter(ExecUtil.formOutputPath(outputPrefix, sampleId, "annot"))

        writer.println(sw.getFullHeader(sample) + "\tid.in.sample\tscore\t" + database.header)

        results.sort { -it.key.count }.each { result ->
            result.value.each {
                writer.println(sw.getFullClonotypeString(result.key) + "\t" +
                        it.id + "\t" +
                        it.result.score + "\t" + it.row.toString())
            }
        }

        writer.close()

        println "[${new Date()} $scriptName] Summarizing..."

        def summary = new ClonotypeSearchSummary(results, sample, summaryColumns)

        def summaryPrefix = sampleId + "\t" + sample.sampleMetadata.toString()

        summary.fieldCounters.each { kvp1 ->
            def columnName = kvp1.key

            kvp1.value.each { kvp2 ->
                pwSummary.println(summaryPrefix + "\t" + columnName + "\t" + kvp2.key + "\t" + kvp2.value)
            }
        }

        pwSummary.println(summaryPrefix + "\tsummary\tfound\t" + summary.totalCounter)
        pwSummary.println(summaryPrefix + "\tsummary\tnot.found\t" + summary.notFoundCounter)

        println "[${new Date()} $scriptName] ${ind + 1} sample(s) done."
    }
}

sampleCollection.metadataTable.storeWithOutput(outputPrefix, compress,
        "vdjdb:${filterStr ?: "all"}")

println "[${new Date()} $scriptName] Finished."