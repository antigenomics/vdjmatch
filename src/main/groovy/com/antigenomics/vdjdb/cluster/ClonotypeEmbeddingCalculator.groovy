package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjtools.misc.ExecUtil
import groovyx.gpars.GParsPool
import smile.graph.Graph
import smile.math.Math

import java.util.concurrent.ConcurrentLinkedQueue

class ClonotypeEmbeddingCalculator {

    private ClonotypeEmbeddingCalculator() {
    }

    /**
     *
     * @param clonotypeGraph
     * @param minComponentSize
     * @param infoWeight
     * @param d
     * @return
     */
    static Collection<ClonotypeEmbedding> isoMap(ClonotypeGraph clonotypeGraph,
                                           int minComponentSize = 5, int d = 2) {
        def clonotypeEmbeddings = new ConcurrentLinkedQueue<ClonotypeEmbedding>()

        // iterate over connected components
        GParsPool.withPool ExecUtil.THREADS, {
            clonotypeGraph.getConnectedComponents().eachParallel { component ->
                int n = component.index.length
                Graph graph = component.graph

                double[][] coords

                if (n >= minComponentSize) {
                    // a-la C-Isomap
                    /*if (infoWeight) {
                        for (Graph.Edge edge : graph.getEdges()) {
                            double d1 = graph.getDegree(edge.v1),
                                   d2 = graph.getDegree(edge.v2)
                            edge.weight *= (Math.log2(d1 + 1) + Math.log2(d2 + 1)) / 2d
                        }
                    }*/

                    // isoMap algorithm
                    coords = EmbeddingHelper.isoMap(graph, d)
                } else {
                    // return dummy coordinates for components that are too small
                    coords = new double[n][d]
                }

                for (int i = 0; i < n; i++) {
                    // id as in sample
                    int id = clonotypeGraph.index[i]

                    // add component index to coordinate vector
                    double[] vector = new double[d + 1]
                    vector[0] = (double) component.componentIndex
                    for (int j = 0; j < d; j++) {
                        vector[j + 1] = coords[i][j]
                    }

                    // create embedding
                    clonotypeEmbeddings.add(new ClonotypeEmbedding(vector,
                            clonotypeGraph.sample[id],
                            id
                    ))
                }
            }
        }

        clonotypeEmbeddings
    }
}