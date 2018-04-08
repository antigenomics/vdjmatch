package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.cluster.ClonotypeDistanceCalculator
import com.antigenomics.vdjdb.cluster.ClonotypeEmbeddingCalculator
import com.antigenomics.vdjdb.cluster.ClonotypeGraph
import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.ScoringProvider
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection

def scriptName = getClass().canonicalName.split("\\.")[-1]

def cli = new CliBuilder(usage: "cluster [options] " +
        "sample output_prefix\n" +
        "Input samples should be provided in VDJtools format if --software is not " +
        "specified/supported.")

def DEFAULT_SEARCH_SCOPE = "0,0,0",
    DEFAULT_EXHAUSTIVE = "1",
    DEFAULT_SCORING_MODE = "1",
    ALLOWED_SPECIES_ALIAS = ["human" : "homosapiens", "mouse": "musmusculus",
                             "monkey": "macacamulatta"],
    ALLOWED_GENES = ["TRA", "TRB"]


cli._(longOpt: "software", argName: "string", args: 1,
        "Input RepSeq data format. Currently supported: ${Software.values().join(", ")}. " +
                "[default = ${Software.VDJtools}]")

cli.S(longOpt: "species", argName: "name", args: 1, required: true,
        "Species of input sample(s), allowed values: ${ALLOWED_SPECIES_ALIAS.keySet()}.")
cli.R(longOpt: "gene", argName: "name", args: 1, required: true,
        "Receptor gene of input sample(s), allowed values: $ALLOWED_GENES.")

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

cli._(longOpt: "hit-weight-inf",
        "Weight database hits by their 'informativeness', i.e. the log probability of them " +
                "being matched by chance.")

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

def optSpecies = (String) opt.S, optGene = (String) opt.R
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

/* initial search */

def optVMatch = (boolean) opt.'v-match', // todo
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

/* scoring */

def optVdjmatchScoring = (boolean) opt.'scoring-vdjmatch',
    optScoringMode = (opt.'scoring-mode' ?: DEFAULT_SCORING_MODE).toInteger()
def scoringBundle = optVdjmatchScoring ?
        ScoringProvider.loadScoringBundle(optSpecies, optGene,
                optScoringMode == 0) :
        ScoringBundle.DUMMY

def optWeightByInfo = (boolean) opt.'hit-weight-inf'

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)..."

def sampleCollection = optMetadataFileName ?
        new SampleCollection((String) optMetadataFileName, optSoftware) :
        new SampleCollection(opt.arguments()[0..-2], optSoftware)

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) to process."

//
// Run pairwise distance calc, graph construction & embedding for all samples
//

println "[${new Date()} $scriptName] Clustering sample(s) & writing results."

def sw = new SampleWriter(optCompress)

//new File(ExecUtil.formOutputPath(outputPrefix, "embed")).withPrintWriter { pwSummary ->
sampleCollection.eachWithIndex { Sample sample, int ind ->
    println "[${new Date()} $scriptName] Clustering..."

    def sampleId = sample.sampleMetadata.sampleId

    ///////// Clustering

    println "[${new Date()} $scriptName] Computing pairwise distances / clonotype graph"

    def distanceCalc = new ClonotypeDistanceCalculator(searchScope, scoringBundle)

    def distances = distanceCalc.computeDistances(sample)

    println "[${new Date()} $scriptName] Constructing clonotype graph"

    def graph = new ClonotypeGraph(sample, distances)

    println "[${new Date()} $scriptName] Embedding clonotypes"

    def embeddings = ClonotypeEmbeddingCalculator.isoMap(graph, 5, optWeightByInfo)

    ///////// Output

    def writer = sw.getWriter(ExecUtil.formOutputPath(outputPrefix, sampleId))

    writer.println(sw.getFullHeader(sample) +
            "\tid.in.sample\tcomponent\tx\ty")

    embeddings.each {
        writer.println([sw.getFullClonotypeString(it.clonotype), it.id, (int) it.coordinates[0],
                        it.coordinates[1], it.coordinates[2]].join("\t"))
    }

    writer.close()

    println "[${new Date()} $scriptName] ${ind + 1} sample(s) done."
}
//}

sampleCollection.metadataTable.storeWithOutput(outputPrefix, optCompress,
        "vdjdb:embedding")

println "[${new Date()} $scriptName] Finished."