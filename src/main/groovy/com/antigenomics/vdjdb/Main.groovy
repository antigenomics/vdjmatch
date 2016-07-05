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

import com.antigenomics.vdjdb.stat.ClonotypeSearchSearchSummary
import com.antigenomics.vdjdb.stat.Counter
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.misc.ExecUtil

if (args.length > 0 && args[0].toLowerCase() == "update") {
    Util.checkDatabase(true)
    System.exit(0)
}

def DEFAULT_PARAMETERS = "2,1,1,2", DEFAULT_CONFIDENCE_THRESHOLD = "2",
    ALLOWED_SPECIES_ALIAS = ["human": "homosapiens", "mouse": "musmusculus",
                             "rat"  : "rattusnorvegicus", "monkey": "macacamulatta"],
    ALLOWED_GENES = ["TRA", "TRB"]

def cli = new CliBuilder(usage: "vdjdb [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix\n" +
        "Input samples should be provided in VDJtools format. See VDJtools/Convert utility.")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli._(longOpt: "search-params", argName: "s,i,d,t", args: 1,
        "CDR3 sequence search parameters: " +
                "allowed number of substitutions (s), insertions (i), deletions (d) and total number of mutations. " +
                "[default=$DEFAULT_PARAMETERS]")
cli._(longOpt: "database", argName: "string", args: 1, "Path and prefix of an external database. " +
        "The prefix should point to a '.txt' file (database itself) and '.meta.txt' (database column metadata).")
cli._(longOpt: "use-fat-db", "Use a more redundant database version, with extra fields (meta, method, etc). " +
        "Fat database can contain several records for a TCR:pMHC pair corresponding to different replicates/tissue sources/targets.")
//cli._(longOpt: "summary", argName: "col1,col2,...", args: 1,
//        "Table columns for summarizing, e.g. origin,disease.type,disease,source for default database.")
cli._(longOpt: "filter", argName: "logical expression(__field__,...)", args: 1,
        "Logical filter evaluated for database columns. Supports Regex, .contains(), .startsWith(), etc.")
cli.S(longOpt: "species", argName: "name", args: 1, required: true,
        "Species of input sample(s), allowed values: ${ALLOWED_SPECIES_ALIAS.keySet()}.")
cli._(longOpt: "software", argName: "string", required: true, args: 1,
        "Input RepSeq data format. Currently supported: ${Software.values().join(", ")}. Default is VDJtools format.")
cli.R(longOpt: "gene", argName: "name", args: 1, required: true,
        "Receptor gene of input sample(s), allowed values: $ALLOWED_GENES.")
cli._(longOpt: "vdjdb-conf-threshold", argName: "[0,7]", args: 1,
        "VDJdb confidence level threshold, [default=$DEFAULT_CONFIDENCE_THRESHOLD]")
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

// Remaining arguments

def dbPrefix = (String) (opt.'database' ?: null),
    p = (opt.'search-params' ?: DEFAULT_PARAMETERS).split(",").collect { it.toInteger() },
    summaryCols = ["antigen.species", "antigen.gene"],//opt.'summary' ? ((String) opt.'summary').split(",") as List<String> : [],
    compress = (boolean) opt.c,
    vMatch = (boolean) opt."v-match", jMatch = (boolean) opt."j-match",
    species = (String) opt.S, gene = (String) opt.R,
    q = (opt.'vdjdb-conf-threshold' ?: DEFAULT_CONFIDENCE_THRESHOLD).toInteger(),
    filterStr = opt.'filter',
    useFatDb = (boolean) opt.'use-fat-db',
    outputPrefix = opt.arguments()[-1],
    software = opt.'software' ? Software.byName(opt.'software') : Software.VDJtools

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

println "[${new Date()} $scriptName] Loaded database. \n${database.dbInstance}"

// Expression filtering if specified
if (filterStr) {
    println "[${new Date()} $scriptName] Filtering using $filterStr."
    database = database.filter(filterStr)
    println "[${new Date()} $scriptName] Done. \n${database.dbInstance}"
}

println "[${new Date()} $scriptName] Preparing clonotype database for $species $gene."

database = database.asClonotypeDatabase(vMatch, jMatch, p[0], p[1], p[2], p[3], species, gene, q)

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

new File(ExecUtil.formOutputPath(outputPrefix, "annot", "stats")).withPrintWriter { pwStats ->
    new File(ExecUtil.formOutputPath(outputPrefix, "annot", "summary")).withPrintWriter { pwSummary ->
        def headerPrefix = [MetadataTable.SAMPLE_ID_COLUMN,
                            sampleCollection.metadataTable.columnHeader]
        pwSummary.println([headerPrefix,
                           "counter.type", summaryCols, "summary.count"].
                flatten().join("\t"))
        pwStats.println([headerPrefix,
                         "database",
                         "species", "gene", "counter.type",
                         "not.found", "found.once", "found.twice.and.more"].
                flatten().join("\t"))

        sampleCollection.eachWithIndex { Sample sample, int ind ->
            def sampleId = sample.sampleMetadata.sampleId

            def results = database.search(sample)

            def writer = sw.getWriter(ExecUtil.formOutputPath(outputPrefix, sampleId, "annot"))

            writer.println(sw.header + "\tsample_id\tscore\t" + database.header)

            results.sort { -it.key.count }.each { result ->
                result.value.each {
                    writer.println(sw.getClonotypeString(result.key) + "\t" +
                            it.result.id + "\t" +
                            it.result.score + "\t" + it.row.toString())
                }
            }

            writer.close()


            def summary = new ClonotypeSearchSearchSummary(database, summaryCols as List<String>, sample)
            summary.append(results)

            def prefix = [sampleId, sample.sampleMetadata.toString()],
                prefix1 = [prefix, dbPrefix ?: "default", species, gene].flatten()

            ["unique", "weighted"].each { counterType ->
                summary.listTopCombinations().each { List<String> combination ->
                    pwSummary.println([
                            prefix, counterType, combination,
                            summary.getCombinationCounter(combination).collect { Counter it -> it."${counterType}Count" }
                    ].flatten().join("\t"))
                }

                pwStats.println([
                        prefix1, counterType,
                        [summary.notFound, summary.foundOnce,
                         summary.foundTwiceAndMore].collect { Counter it -> it."${counterType}Count" }
                ].flatten().join("\t"))
            }

            println "[${new Date()} $scriptName] ${ind + 1} sample(s) done."
        }
    }
}

println "[${new Date()} $scriptName] Finished."