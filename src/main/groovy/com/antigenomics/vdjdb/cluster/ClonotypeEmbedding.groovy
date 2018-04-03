package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjtools.sample.Clonotype

class ClonotypeEmbedding {
    final double [] coordinates
    final Clonotype clonotype
    final int id

    ClonotypeEmbedding(double[] coordinates, Clonotype clonotype, int id) {
        this.coordinates = coordinates
        this.clonotype = clonotype
        this.id = id
    }
}
