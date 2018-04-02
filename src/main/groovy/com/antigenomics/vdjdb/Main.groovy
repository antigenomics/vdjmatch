/*
 * Copyright 2015-2017 Mikhail Shugay
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

import com.antigenomics.vdjdb.impl.ClonotypeSearchResult
import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.ScoringProvider
import com.antigenomics.vdjdb.impl.filter.DummyResultFilter
import com.antigenomics.vdjdb.impl.filter.MaxScoreResultFilter
import com.antigenomics.vdjdb.impl.filter.ScoreThresholdResultFilter
import com.antigenomics.vdjdb.impl.filter.TopNResultFilter
import com.antigenomics.vdjdb.impl.weights.DegreeWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjdb.stat.ClonotypeCounter
import com.antigenomics.vdjdb.stat.ClonotypeSearchSummary
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.misc.ExecUtil

def scriptName = getClass().canonicalName.split("\\.")[-1]

if (args.length > 0 && args[0].toLowerCase() == "update") {
    Util.checkDatabase(true)
    System.exit(0)
}

//
// Command line options
//

def cli = new CliBuilder(usage: "vdjdb [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix\n" +
        "Input samples should be provided in VDJtools format if --software is not " +
        "specified/supported.")

def DEFAULT_SEARCH_SCOPE = "0,0,0",
    DEFAULT_EXHAUSTIVE = "1",
    DEFAULT_SCORING_MODE = "1",
    DEFAULT_CONFIDENCE_THRESHOLD = "0",
    ALLOWED_SPECIES_ALIAS = ["human" : "homosapiens", "mouse": "musmusculus",
                             "monkey": "macacamulatta"],
    ALLOWED_GENES = ["TRA", "TRB"]

cli.h("Displays help message")


cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli._(longOpt: "software", argName: "string", args: 1,
        "Input RepSeq data format. Currently supported: ${Software.values().join(", ")}. " +
                "[default = ${Software.VDJtools}]")
cli.c("Compressed output")


cli.S(longOpt: "species", argName: "name", args: 1, required: true,
        "Species of input sample(s), allowed values: ${ALLOWED_SPECIES_ALIAS.keySet()}.")
cli.R(longOpt: "gene", argName: "name", args: 1, required: true,
        "Receptor gene of input sample(s), allowed values: $ALLOWED_GENES.")
cli._(longOpt: "vdjdb-conf", argName: "0..3", args: 1,
        "VDJdb confidence level threshold, from lowest to highest. [default=$DEFAULT_CONFIDENCE_THRESHOLD]")
cli._(longOpt: "filter", argName: "logical expression(__field__,...)", args: 1,
        "[advanced] Logical filter evaluated for database columns. " +
                "Supports Regex, .contains(), .startsWith(), etc.")


cli._(longOpt: "database", argName: "string", args: 1,
        "[advanced] Path and prefix of an external database. " +
                "The prefix should point to a '.txt' file (database itself) and " +
                "'.meta.txt' (database column metadata).")
cli._(longOpt: "use-fat-db",
        "[advanced] Use a more redundant database version, with extra fields (meta, method, etc). " +
                "Fat database can contain several records for a " +
                "TCR:pMHC pair corresponding to different replicates/tissue sources/targets.")


cli._(longOpt: "v-match", "Require exact (up to allele) V segment id matching.")
cli._(longOpt: "j-match", "Require exact (up to allele) J segment id matching.")
cli.O(longOpt: "search-scope", argName: "s,id,t or s,i,d,t", args: 1,
        "Sets CDR3 sequence matching parameters aka 'search scope': " +
                "allowed number of substitutions (s), insertions (i), deletions (d) / or indels (id) and " +
                "total number of mutations (t). [default=$DEFAULT_SEARCH_SCOPE]")
cli._(longOpt: "search-exhaustive", argName: "0..2", args: 1,
        "Perform exhaustive CDR3 alignment: 0 - no (fast), " +
                "1 - check and select best alignment for smallest edit distance, " +
                "2 - select best alignment across all edit distances within search scope (slow). " +
                "[default=$DEFAULT_EXHAUSTIVE]")


cli.A(longOpt: "scoring-vdjmatch",
        "Use VDJMATCH algorithm that computes full alignment score as a function of " +
                "CDR3 mutations (weighted with VDJAM scoring matrix) and pre-computed V/J segment " +
                "match scores. If not set, will just count the number of mismatches.")
cli._(longOpt: "scoring-mode", argName: "0..1", args: 1,
        "Either '0': scores mismatches only (faster) or '1': compute scoring for whole sequences (slower). " +
                "[default=$DEFAULT_SCORING_MODE]")


cli.T(longOpt: "hit-filter-score", argName: "threshold", args: 1,
        "Drops hits with a score less than the specified threshold.")
cli.X(longOpt: "hit-filter-max",
        "Only select hit with maximal score for a given query clonotype " +
                "(will consider all max score hits in case of ties).")
cli._(longOpt: "hit-filter-topn", argName: "n", args: 1,
        "Select best 'n' hits by score " +
                "(can randomly drop hits in case of ties).")


cli._(longOpt: "hit-weight-inf",
        "Weight database hits by their 'informativeness', i.e. the log probability of them " +
                "being matched by chance.")

cli._(longOpt: "summary-columns", argName: "col1,col2,...", args: 1,
        "Table columns for summarizing, see DB metadata for column names. " +
                "[default=${ClonotypeSearchSummary.FIELDS_PLAIN_TEXT.join(",")}]. ")

//
// Parse arguments & run checks
//

def opt = cli.parse(args)
if (opt == null) {
    System.exit(1)
}
if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(1)
}

/* software type, metadata, input files & output options */

def optMetadataFileName = opt.m,
    optSoftware = opt.'software' ? Software.byName((String) opt.'software') : Software.VDJtools,
    optCompress = (boolean) opt.c,
    outputPrefix = opt.arguments()[-1]
if (optMetadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (optMetadataFileName)
        println "[ERROR] Only output prefix should be provided in case of -m"
    else
        println "[ERROR] At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(1)
}

if (optMetadataFileName) {
    println "Will use metadata from " + optMetadataFileName
} else {
    println "Will use provided sample files: " + opt.arguments()[0..-2]
}

/* species, gene, pre-filtering */

def optSpecies = (String) opt.S, optGene = (String) opt.R,
    optVdjdbConf = (opt.'vdjdb-conf' ?: DEFAULT_CONFIDENCE_THRESHOLD).toInteger(),
    optFilterString = (String) (opt.'filter' ?: null)
def allowedSpecies = [ALLOWED_SPECIES_ALIAS.keySet(), ALLOWED_SPECIES_ALIAS.values()].flatten()
if (!allowedSpecies.any { optSpecies.equalsIgnoreCase((String) it) }) {
    println "Wrong species name, use one of ${allowedSpecies} (case-insensitive)"
    System.exit(1)
}
optSpecies = ALLOWED_SPECIES_ALIAS[optSpecies.toLowerCase()] ?: optSpecies
if (!ALLOWED_GENES.any { optGene.equalsIgnoreCase(it) }) {
    println "Wrong gene name, use one of $ALLOWED_GENES (case-insensitive)"
    System.exit(1)
}

/* advanced db setup */

def dbPrefix = (String) (opt.'database' ?: null),
    useFatDb = (boolean) opt.'use-fat-db'

/* initial search */

def optVMatch = (boolean) opt.'v-match',
    optJMatch = (boolean) opt.'j-match',
    optSearchScope = (opt.'search-scope' ?: DEFAULT_SEARCH_SCOPE).split(",").collect { it.toInteger() },
    optExhaustive = (opt.'search-exhaustive' ?: DEFAULT_EXHAUSTIVE).toInteger()

if (!opt.'scoring-vdjmatch') {
    optExhaustive = 0 // has no effect in case no scoring is used
}

def searchScope = optSearchScope.size() == 3 ?
        new SearchScope(optSearchScope[0], optSearchScope[1], optSearchScope[2],
                optExhaustive > 0, optExhaustive < 2)
        :
        new SearchScope(optSearchScope[0], optSearchScope[1], optSearchScope[2], optSearchScope[3],
                optExhaustive > 0, optExhaustive < 2)

//println

/* scoring */

def optVdjmatchScoring = (boolean) opt.'scoring-vdjmatch',
    optScoringMode = (opt.'scoring-mode' ?: DEFAULT_SCORING_MODE).toInteger()
def scoringBundle = optVdjmatchScoring ?
        ScoringProvider.loadScoringBundle(optSpecies, optGene,
                optScoringMode == 0) :
        ScoringBundle.DUMMY

/* filtering */

def optScoreThreshold = (opt.'hit-filter-score' ?: "-Infinity").toFloat()
def resultFilter
if (opt.'hit-filter-max') {
    resultFilter = new MaxScoreResultFilter(optScoreThreshold)
} else if (opt.'hit-filter-topn') {
    resultFilter = new TopNResultFilter(optScoreThreshold,
            (int) (opt.'hit-filter-topn').toInteger())
} else if (opt.'hit-filter-score') {
    resultFilter = new ScoreThresholdResultFilter(optScoreThreshold)
} else {
    resultFilter = DummyResultFilter.INSTANCE
}

/* weighting */

def optWeightByInfo = opt.'hit-weight-inf'
def weightFunctionFactory = optWeightByInfo ?
        DegreeWeightFunctionFactory.DEFAULT :
        DummyWeightFunctionFactory.INSTANCE

/* summary */

def summaryColumns = (opt.'summary-columns' ?: ClonotypeSearchSummary.FIELDS_PLAIN_TEXT.join(",")).split(",") as List<String>

//
// Database loading and filtering
//

println "[${new Date()} $scriptName] Loading database..."

def vdjdbInstance

if (dbPrefix) {
    /* load from specified path */
    def metaStream = new FileInputStream("${dbPrefix}.meta.txt"),
        dataStream = new FileInputStream("${dbPrefix}.txt")
    vdjdbInstance = new VdjdbInstance(metaStream, dataStream)
} else {
    /* load local */
    vdjdbInstance = new VdjdbInstance(useFatDb)
}

/* Re-check summary columns */

def missingSummaryCols = summaryColumns.findAll { !vdjdbInstance.header*.name.contains(it) }
if (!missingSummaryCols.empty) {
    println "[ERROR] Columns $missingSummaryCols specified for summary generation are missing in the database."
    System.exit(1)
}
println "[${new Date()} $scriptName] Loaded database. \n${vdjdbInstance.dbInstance}"

/* Expression filtering if specified */

if (optFilterString) {
    println "[${new Date()} $scriptName] Filtering using $optFilterString."
    vdjdbInstance = vdjdbInstance.filter(optFilterString)
    println "[${new Date()} $scriptName] Done. \n${vdjdbInstance.dbInstance}"
}

//
// Initialize clonotype database - filter by species and gene, specify search parameters
//

println "[${new Date()} $scriptName] Preparing clonotype database for $optSpecies $optGene."

def clonotypeDatabase = vdjdbInstance.asClonotypeDatabase(
        optSpecies, optGene,
        searchScope,
        scoringBundle,
        weightFunctionFactory,
        resultFilter,
        optVMatch, optJMatch,
        optVdjdbConf)

println "[${new Date()} $scriptName] Done. \n$clonotypeDatabase"

/* Re-check if we have any records */
if (clonotypeDatabase.rows.empty) {
    println "[ERROR] No records present in filtered database"
    System.exit(1)
}

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)..."

def sampleCollection = optMetadataFileName ?
        new SampleCollection((String) optMetadataFileName, optSoftware) :
        new SampleCollection(opt.arguments()[0..-2], optSoftware)

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) to process."

//
// Main loop - processing, summarizing, writing output
//

println "[${new Date()} $scriptName] Annotating sample(s) & writing results."

def sw = new SampleWriter(optCompress)

new File(ExecUtil.formOutputPath(outputPrefix, "annot", "summary")).withPrintWriter { pwSummary ->
    pwSummary.println([MetadataTable.SAMPLE_ID_COLUMN,
                       sampleCollection.metadataTable.columnHeader,
                       "db.column.name",
                       "db.column.value",
                       ClonotypeCounter.HEADER].flatten().join("\t"))

    sampleCollection.eachWithIndex { Sample sample, int ind ->
        println "[${new Date()} $scriptName] Annotating..."

        def sampleId = sample.sampleMetadata.sampleId

        def results = clonotypeDatabase.search(sample)

        def writer = sw.getWriter(ExecUtil.formOutputPath(outputPrefix, sampleId))

        writer.println(sw.getFullHeader(sample) +
                "\tid.in.sample\t" +
                "match.score\tmatch.weight\t" +
                clonotypeDatabase.header)

        results.sort { -it.key.count }.each { matchListEntry ->
            matchListEntry.value.each { ClonotypeSearchResult match ->
                writer.println(sw.getFullClonotypeString(matchListEntry.key) + "\t" +
                        match.id + "\t" +
                        match.score + "\t" +
                        match.weight + "\t" +
                        match.row.toTabDelimitedString())
            }
        }

        writer.close()

        println "[${new Date()} $scriptName] Summarizing..."

        def summary = new ClonotypeSearchSummary(results, sample, summaryColumns, clonotypeDatabase)

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

sampleCollection.metadataTable.storeWithOutput(outputPrefix, optCompress,
        "vdjdb:${optFilterString ?: "all"}")

println "[${new Date()} $scriptName] Finished."