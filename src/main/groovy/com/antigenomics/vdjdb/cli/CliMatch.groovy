package com.antigenomics.vdjdb.cli

import com.antigenomics.vdjdb.stat.ClonotypeSearchSummary

class CliMatch extends CliBase<OptMatch> {
    static def DEFAULT_CONFIDENCE_THRESHOLD = "0", DEFAULT_MIN_EPI_SIZE = "10"

    CliMatch() {
        super("match")

        /* pre-filtering */

        cli._(longOpt: "vdjdb-conf", argName: "0..3", args: 1,
                "VDJdb confidence level threshold, from lowest to highest. [default=$DEFAULT_CONFIDENCE_THRESHOLD]")
        cli._(longOpt: "min-epi-size", argName: "integer", args: 1,
                "Minimal number of unique CDR3 sequences per epitope in VDJdb, " +
                        "filters underrepresented epitopes. [default=$DEFAULT_MIN_EPI_SIZE]")
        cli._(longOpt: "filter", argName: "logical expression(__field__,...)", args: 1,
                "[advanced] Logical filter evaluated for database columns. " +
                        "Supports Regex, .contains(), .startsWith(), etc.")

        /* advanced db setup */

        cli._(longOpt: "database", argName: "string", args: 1,
                "[advanced] Path and prefix of an external database. " +
                        "The prefix should point to a '.txt' file (database itself) and " +
                        "'.meta.txt' (database column metadata).")
        cli._(longOpt: "use-fat-db",
                "[advanced] Use a more redundant database version, with extra fields (meta, method, etc). " +
                        "Fat database can contain several records for a " +
                        "TCR:pMHC pair corresponding to different replicates/tissue sources/targets.")

        /* summary */

        cli._(longOpt: "summary-columns", argName: "col1,col2,...", args: 1,
                "Table columns for summarizing, see DB metadata for column names. " +
                        "[default=${ClonotypeSearchSummary.FIELDS_PLAIN_TEXT.join(",")}]. ")
    }

    @Override
    OptMatch parseArguments(String[] args) {
        new OptMatch(this, args)
    }
}
