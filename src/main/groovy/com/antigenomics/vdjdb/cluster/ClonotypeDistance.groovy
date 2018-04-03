package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjtools.sample.Clonotype

/**
 * Distance (dissimilarity) between a pair of clonotypes
 */
class ClonotypeDistance {
    /**
     * Id of 'from' clonotype in parent sample
     */
    final int from
    /**
     * Id of 'to' clonotype in parent sample
     */
    final int to
    /**
     * 'from' clonotype
     */
    final Clonotype query
    /**
     * 'to' clonotype
     */
    final Clonotype target
    /**
     * TCR alignment score / probability of TCR specificity match
     */
    final float score
    /**
     * Weight aka informativeness of 'from' clonotype (typically undefined)
     */
    final float weight
    /**
     * TCR alignment score transformed into dissimilarity
     */
    final double dissimilarity

    /**
     * Constructor
     * @param from id of 'from' clonotype in parent sample
     * @param to id of 'to' clonotype in parent sample
     * @param query 'from' clonotype
     * @param target 'to' clonotype
     * @param score TCR alignment score / probability of TCR specificity match
     * @param weight informativeness of 'from' clonotype
     * @param dissimilarity TCR alignment score transformed into dissimilarity
     */
    ClonotypeDistance(int from, int to, Clonotype query, Clonotype target,
                      float score, float weight, double dissimilarity) {
        this.from = from
        this.to = to
        this.query = query
        this.target = target
        this.score = score
        this.weight = weight
        this.dissimilarity = dissimilarity
    }


    @Override
    String toString() {
        "(" + from + "," + to + "):" + dissimilarity
    }
}
