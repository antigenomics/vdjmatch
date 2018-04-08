package com.antigenomics.vdjdb.cli

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjtools.misc.Software

abstract class CliBase<T extends OptBase> {
    static def DEFAULT_SEARCH_SCOPE = "0,0,0", // todo: different defaults?
               DEFAULT_EXHAUSTIVE = "1",
               DEFAULT_SCORING_MODE = "1",
               ALLOWED_SPECIES_ALIAS = ["human" : "homosapiens", "mouse": "musmusculus",
                                        "monkey": "macacamulatta"],
               ALLOWED_GENES = ["TRA", "TRB"]

    final String routineName
    final CliBuilder cli

    CliBase(String routineName) {
        this.routineName = routineName

        this.cli = new CliBuilder(usage: "$routineName [options] " +
                "[sample1 sample2 sample3 ... if -m is not specified] output_prefix\n" +
                "Input samples should be provided in VDJtools format if --software is not " +
                "specified/supported.")

        cli.h("Displays help message")

        /* software type, metadata, input files & output options */

        cli.m(longOpt: "metadata", argName: "filename", args: 1,
                "Metadata file. First and second columns should contain file name and sample id. " +
                        "Header is mandatory and will be used to assign column names for metadata.")
        cli._(longOpt: "software", argName: "string", args: 1,
                "Input RepSeq data format. Currently supported: ${Software.values().join(", ")}. " +
                        "[default = ${Software.VDJtools}]")
        cli.c("Compressed output")

        /* species, gene */

        cli.S(longOpt: "species", argName: "name", args: 1, required: true,
                "Species of input sample(s), allowed values: ${ALLOWED_SPECIES_ALIAS.keySet()}.")
        cli.R(longOpt: "gene", argName: "name", args: 1, required: true,
                "Receptor gene of input sample(s), allowed values: $ALLOWED_GENES.")

        /* initial search */

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

        /* scoring */

        cli.A(longOpt: "scoring-vdjmatch",
                "Use VDJMATCH algorithm that computes full alignment score as a function of " +
                        "CDR3 mutations (weighted with VDJAM scoring matrix) and pre-computed V/J segment " +
                        "match scores. If not set, will just count the number of mismatches.")
        cli._(longOpt: "scoring-mode", argName: "0..1", args: 1,
                "Either '0': scores mismatches only (faster) or '1': compute scoring for whole sequences (slower). " +
                        "[default=$DEFAULT_SCORING_MODE]")

        /* filtering */

        cli.T(longOpt: "hit-filter-score", argName: "threshold", args: 1,
                "Drops hits with a score less than the specified threshold.")
        cli.X(longOpt: "hit-filter-max",
                "Only select hit with maximal score for a given query clonotype " +
                        "(will consider all max score hits in case of ties).")
        cli._(longOpt: "hit-filter-topn", argName: "n", args: 1,
                "Select best 'n' hits by score " +
                        "(can randomly drop hits in case of ties).")

        /* weighting */

        cli._(longOpt: "hit-weight-inf",
                "Weight query hits by their 'informativeness', i.e. the log probability of them " +
                        "being matched by chance.")
    }

    abstract T parseArguments(String[] args)

    void progress(String message) {
        Util.sout(routineName, message)
    }
}
