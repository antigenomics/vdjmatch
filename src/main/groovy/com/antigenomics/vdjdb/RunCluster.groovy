package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.cli.CliCluster
import com.antigenomics.vdjdb.cluster.ClonotypeDistanceCalculator
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.sample.Sample

//
// Parse command line arguments, prepare options
//

def opts = new CliCluster().parseArguments(args)

//
// Initialize distance calculator
//

def metric = new ClonotypeDistanceCalculator(opts.searchScope,
        opts.scoringBundle, opts.weightFunctionFactory,
        opts.resultFilter, opts.optVMatch, opts.optJMatch)

//
// Compute pairwise distances - build graph - embed
//

def sw = new SampleWriter(opts.optCompress)

opts.sampleCollection.eachWithIndex { Sample sample, int ind ->
    opts.cliBase.progress("Clustering sample #${ind + 1} of ${opts.sampleCollection.size()}")

    def sampleId = sample.sampleMetadata.sampleId
    def writer = sw.getWriter(ExecUtil.formOutputPath(opts.outputPrefix, sampleId, "clusters"))

    writer.println([sw.getFullHeader(sample).split("\t").collect { "from." + it }.join("\t"),
                    sw.getFullHeader(sample).split("\t").collect { "to." + it }.join("\t"),
                    "from.id",
                    "to.id",
                    "match.score",
                    "match.weight",
                    "dissimilarity"].join("\t"))

    // Distances

    opts.cliBase.progress("Computing distances..")
    def distances = metric.computeDistances(sample)

    distances.each {
        writer.println([sw.getFullClonotypeString(it.query),
                        sw.getFullClonotypeString(it.target),
                        it.from,
                        it.to,
                        it.score,
                        it.weight,
                        it.dissimilarity].join("\t"))
    }

    writer.close()

    // Todo: embedding

    opts.cliBase.progress("Done.")
}

opts.cliBase.progress("Finished.")