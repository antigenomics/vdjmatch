package com.antigenomics.vdjdb.cluster

import com.antigenomics.vdjdb.misc.DistanceEdge
import com.antigenomics.vdjtools.sample.Clonotype

class ClonotypeDistance implements DistanceEdge {
    final int idQuery, idTarget
    final Clonotype query, target
    final float score
    final float weight
    final float weightedScore

    ClonotypeDistance(int idQuery, int idTarget, Clonotype query, Clonotype target,
                      float score, float weight, float weightedScore) {
        this.idQuery = idQuery
        this.idTarget = idTarget
        this.query = query
        this.target = target
        this.score = score
        this.weight = weight
        this.weightedScore = weightedScore
    }


    @Override
    String toString() {
        "(" + idQuery + "," + idTarget + "):" + weightedScore
    }

    @Override
    int getFrom() {
        idQuery
    }

    @Override
    int getTo() {
        idTarget
    }

    @Override
    double getDissimilarity() {
        (1.0d - score) / (double) weight
    }
}
