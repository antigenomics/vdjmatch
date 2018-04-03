package com.antigenomics.vdjdb.cluster;

/*

Re-implemented to handle sparse distance matrix/graph with dissimilarities from smile/Isomap

Copyright (c) 2010 Haifeng Li

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

import smile.graph.Graph;
import smile.math.Math;
import smile.math.matrix.Matrix;
import smile.math.matrix.DenseMatrix;
import smile.math.matrix.EVD;

class EmbeddingHelper {
    public static double[][] isoMap(Graph graph, int d) {
        int n = graph.getNumVertices();

        double[][] D = graph.dijkstra();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                D[i][j] = -0.5 * D[i][j] * D[i][j];
                D[j][i] = D[i][j];
            }
        }

        double[] mean = Math.rowMeans(D);
        double mu = Math.mean(mean);

        DenseMatrix B = Matrix.zeros(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double b = D[i][j] - mean[i] - mean[j] + mu;
                B.set(i, j, b);
                B.set(j, i, b);
            }
        }

        B.setSymmetric(true);

        EVD eigen = B.eigen(d);

        DenseMatrix V = eigen.getEigenVectors();
        double[][] coordinates = new double[n][d];

        if (eigen.getEigenValues().length < d) {
            d = eigen.getEigenValues().length;
        }

        for (int j = 0; j < d; j++) {
            if (eigen.getEigenValues()[j] < 0) {
                throw new IllegalArgumentException(String.format("Some of the first %d eigenvalues are < 0.", d));
            }

            double scale = Math.sqrt(eigen.getEigenValues()[j]);
            for (int i = 0; i < n; i++) {
                coordinates[i][j] = V.get(i, j) * scale;
            }
        }

        return coordinates;
    }
}
