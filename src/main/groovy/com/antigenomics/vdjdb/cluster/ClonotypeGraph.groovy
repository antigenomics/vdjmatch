package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjtools.sample.Sample
import smile.graph.AdjacencyList

/**
 * Graph representation of sparse clonotype dissimilarity matrix
 */
class ClonotypeGraph {
    /**
     * Connected component index. Zero indicates full graph, its sub-graphs are indexed from 1 to N
     */
    final int componentIndex
    /**
     * Parent sample of clonotype graph
     */
    final Sample sample
    /**
     * Internal graph object
     */
    final AdjacencyList graph
    /**
     * Array of clonotype IDs. Stores 1-to-1 correspondence between vertex index (array element index) and
     * the ID of clonotype in the sample.
     */
    final int[] index

    /**
     * List of pairwise dissimilarities between clonotypes
     */
    final List<ClonotypeDistance> clonotypeDistanceList

    private List<ClonotypeGraph> components = null

    /**
     * INTERNAL - for subgraph initialization
     */
    private ClonotypeGraph(Sample sample, int componentIndex, int[] index, AdjacencyList graph) {
        this.componentIndex = componentIndex
        this.sample = sample
        this.index = index
        this.graph = graph
    }

    /**
     * Create a clonotype graph for a given sample and distances (dissimilarities) between clonotypes
     * @param sample parent sample
     * @param clonotypeDistanceList list of pairwise distances
     */
    ClonotypeGraph(Sample sample,
                   List<ClonotypeDistance> clonotypeDistanceList) {
        this.clonotypeDistanceList = clonotypeDistanceList
        this.componentIndex = 0
        this.sample = sample
        this.index = new int[sample.diversity]
        this.graph = new AdjacencyList(sample.diversity)

        for (int i = 0; i < sample.diversity; i++) {
            index[i] = i
        }

        clonotypeDistanceList.each {
            graph.setWeight(it.from, it.to, it.dissimilarity)
        }
    }

    /**
     * Get connected components of the graph (subgraphs as list)
     * @return
     */
    List<ClonotypeGraph> getConnectedComponents() {
        if (components == null) {
            if (componentIndex > 0) {
                throw new IllegalStateException("Cannot call connected components from subgraph")
            }

            components = new ArrayList<ClonotypeGraph>()
            int[][] cc = graph.bfs()

            for (int i = 0; i < cc.length; i++) {
                int[] copyIndex = cc[i]

                components.add(new ClonotypeGraph(sample,
                        i + 1, copyIndex, graph.subgraph(copyIndex)))
            }
        }

        components
    }
}
