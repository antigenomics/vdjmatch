package com.antigenomics.vdjdb.cli

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.VdjdbInstance
import com.antigenomics.vdjdb.stat.ClonotypeSearchSummary

class OptMatch extends OptBase {
    final VdjdbInstance vdjdbInstance
    final int optVdjdbConf
    final String optFilterString
    final List<String> summaryColumns

    OptMatch(CliMatch cliBase, String[] args) {
        super(cliBase, args)

        // Use alternative database?

        def dbPrefix = (String) (opt.'database' ?: null),
            useFatDb = (boolean) opt.'use-fat-db'

        // Load database

        Util.info("Loading specificity database...")

        if (dbPrefix) {
            /* load from specified path */
            def metaStream = new FileInputStream("${dbPrefix}.meta.txt"),
                dataStream = new FileInputStream("${dbPrefix}.txt")
            vdjdbInstance = new VdjdbInstance(metaStream, dataStream)
        } else {
            /* load local */
            vdjdbInstance = new VdjdbInstance(useFatDb)
        }

        Util.info("Loaded database.\n${vdjdbInstance.dbInstance}")

        // Additional VDJdb record filtering options

        optVdjdbConf = (opt.'vdjdb-conf' ?: CliMatch.DEFAULT_CONFIDENCE_THRESHOLD).toInteger()
        optFilterString = (String) (opt.'filter' ?: null)

        // List of summary columns

        summaryColumns = (opt.'summary-columns' ?: ClonotypeSearchSummary.FIELDS_PLAIN_TEXT.join(","))
                .split(",") as List<String>

        /* Re-check summary columns */

        def missingSummaryCols = summaryColumns.findAll { !vdjdbInstance.header*.name.contains(it) }
        if (!missingSummaryCols.empty) {
            Util.error("Columns $missingSummaryCols specified for summary generation are missing in the database.")
        }
    }
}
