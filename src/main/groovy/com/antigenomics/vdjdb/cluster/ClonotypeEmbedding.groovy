package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjtools.sample.Clonotype

/**
 * Clonotype real vector representation
 */
class ClonotypeEmbedding {
    /**
     * Clonotype vector
     */
    final double[] coordinates
    /**
     * Clonotype
     */
    final Clonotype clonotype
    /**
     * Clonotype id in parent sample
     */
    final int id

    /**
     * Constructor
     * @param coordinates real vector
     * @param clonotype clonotype
     * @param id id in parent sample
     */
    ClonotypeEmbedding(double[] coordinates, Clonotype clonotype, int id) {
        this.coordinates = coordinates
        this.clonotype = clonotype
        this.id = id
    }
}
