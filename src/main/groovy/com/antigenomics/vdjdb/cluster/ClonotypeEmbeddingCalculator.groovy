package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.sample.Sample
import smile.graph.Graph
import smile.math.Math

class ClonotypeEmbeddingCalculator {
    final ClonotypeDistanceCalculator clonotypeDistanceCalculator

    ClonotypeEmbeddingCalculator(Sample sample,
                                 SearchScope searchScope = SearchScope.EXACT,
                                 ScoringBundle scoringBundle = ScoringBundle.DUMMY) {
        this.clonotypeDistanceCalculator = new ClonotypeDistanceCalculator(sample,
                searchScope, scoringBundle, DummyWeightFunctionFactory.INSTANCE)
    }

    List<ClonotypeEmbedding> isoMap(int minComponentSize = 5, boolean infoWeight = true, int d = 2) {
        def clonotypeEmbeddings = new ArrayList<ClonotypeEmbedding>()

        // compute full clonotype graph
        def clonotypeGraph = clonotypeDistanceCalculator.computeClonotypeGraph()

        // iterate over connected components
        clonotypeGraph.getConnectedComponents().each { component ->
            int n = component.index.length
            def graph = component.graph

            double[][] coords

            if (n >= minComponentSize) {
                // a-la C-Isomap
                if (infoWeight) {
                    for (Graph.Edge edge : graph.getEdges()) {
                        double d1 = graph.getDegree(edge.v1),
                               d2 = graph.getDegree(edge.v2);
                        edge.weight *= (Math.log2(d1 + 1) + Math.log2(d2 + 1)) / 2;
                    }
                }

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
                        clonotypeDistanceCalculator.sample[id],
                        id
                ))
            }
        }

        clonotypeEmbeddings
    }
}