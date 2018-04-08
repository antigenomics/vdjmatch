package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.cli.CliCluster
import com.antigenomics.vdjdb.cluster.ClonotypeDistanceCalculator
import com.antigenomics.vdjdb.cluster.ClonotypeEmbeddingCalculator
import com.antigenomics.vdjdb.cluster.ClonotypeGraph
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
    def writer = sw.getWriter(ExecUtil.formOutputPath(opts.outputPrefix, sampleId, "distances"))

    writer.println([sw.getFullHeader(sample).split("\t").collect { "from." + it }.join("\t"),
                    sw.getFullHeader(sample).split("\t").collect { "to." + it }.join("\t"),
                    "from.id",
                    "to.id",
                    "match.score",
                    "match.weight",
                    "dissimilarity"].join("\t"))

    // Distances

    opts.cliBase.progress("Computing and writing distances..")
    def distances = metric.computeDistances(sample)

    opts.cliBase.progress("Got " + distances.size() + " edges")

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

    if (opts.optIsoMds) {
        // Compute graph
        opts.cliBase.progress("Computing graph..")
        def graph = new ClonotypeGraph(sample, distances)
        opts.cliBase.progress("Got " + graph.connectedComponents.size() + " connected components")

        // Compute embeddings and write them
        opts.cliBase.progress("Computing and writing embeddings using ISOMAP..")
        def embeddings = ClonotypeEmbeddingCalculator.isoMap(graph, opts.optIsoMdsMinCompSz, opts.optIsoMdsD)

        writer = sw.getWriter(ExecUtil.formOutputPath(opts.outputPrefix, sampleId, "isomds"))
        writer.println([sw.getFullHeader(sample),
                        "component", "x", "y"].join("\t"))
        embeddings.each {
            writer.println([sw.getFullClonotypeString(it.clonotype),
                            (int) it.coordinates[0],
                            it.coordinates[1],
                            it.coordinates[2]].join("\t"))
        }
        writer.close()
    }

    opts.cliBase.progress("Done.")
}

opts.cliBase.progress("Finished.")