package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.cli.CliMatch
import com.antigenomics.vdjdb.impl.ClonotypeSearchResult
import com.antigenomics.vdjdb.stat.ClonotypeCounter
import com.antigenomics.vdjdb.stat.ClonotypeSearchSummary
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

//
// Parse command line arguments, prepare options
//

def opts = new CliMatch().parseArguments(args)

//
// Expression filtering of VDJdb records if specified
//

def vdjdbInstance = opts.vdjdbInstance
if (opts.optFilterString) {
    opts.cliBase.progress("Filtering using $opts.optFilterString")
    vdjdbInstance = vdjdbInstance.filter(opts.optFilterString)
    opts.cliBase.progress("Done.\n${vdjdbInstance.dbInstance}")
}

//
// Initialize clonotype database - filter by species and gene, specify search parameters
//

opts.cliBase.progress("Preparing clonotype database for $opts.optSpecies $opts.optGene.")
def clonotypeDatabase = vdjdbInstance.asClonotypeDatabase(
        opts.optSpecies, opts.optGene,
        opts.searchScope,
        opts.scoringBundle,
        opts.weightFunctionFactory,
        opts.resultFilter,
        opts.optVMatch, opts.optJMatch,
        opts.optVdjdbConf, opts.optMinEpiSize)
opts.cliBase.progress("Done.\n$clonotypeDatabase")

/* Re-check if we have any records */
if (clonotypeDatabase.rows.empty) {
    Util.error("No records present in filtered database.")
}

//
// Main loop - processing, summarizing, writing output
//

opts.cliBase.progress("Annotating sample(s) & writing results.")

def sw = new SampleWriter(opts.optCompress)

new File(ExecUtil.formOutputPath(opts.outputPrefix, "annot", "summary")).withPrintWriter { pwSummary ->
    pwSummary.println([MetadataTable.SAMPLE_ID_COLUMN,
                       opts.sampleCollection.metadataTable.columnHeader,
                       "db.column.name",
                       "db.column.value",
                       ClonotypeCounter.HEADER].flatten().join("\t"))

    opts.sampleCollection.eachWithIndex { Sample sample, int ind ->
        opts.cliBase.progress("Annotating sample #${ind + 1} of ${opts.sampleCollection.size()}")

        def sampleId = sample.sampleMetadata.sampleId

        // Annotate

        def results = clonotypeDatabase.search(sample)

        // Write results

        def writer = sw.getWriter(ExecUtil.formOutputPath(opts.outputPrefix, sampleId))

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

        // Summarise

        opts.cliBase.progress("Summarising..")

        def summary = new ClonotypeSearchSummary(results, sample, opts.summaryColumns, clonotypeDatabase)

        def summaryPrefix = sampleId + "\t" + sample.sampleMetadata.toString()

        summary.fieldCounters.each { kvp1 ->
            def columnName = kvp1.key

            kvp1.value.each { kvp2 ->
                pwSummary.println(summaryPrefix + "\t" + columnName + "\t" + kvp2.key + "\t" + kvp2.value)
            }
        }

        pwSummary.println(summaryPrefix + "\tsummary\tfound\t" + summary.totalCounter)
        pwSummary.println(summaryPrefix + "\tsummary\tnot.found\t" + summary.notFoundCounter)

        opts.cliBase.progress("Done.")
    }
}

opts.sampleCollection.metadataTable.storeWithOutput(opts.outputPrefix,
        opts.optCompress,
        "vdjdb:${opts.optFilterString ?: "all"}")

opts.cliBase.progress("Finished.")