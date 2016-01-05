package com.antigenomics.vdjdb.fixer

import com.antigenomics.vdjdb.Util

class Segment {
    final String seq

    Segment(String seq, int refPoint, boolean variable) {
        seq =  variable ? seq.substring(refPoint + 1) : seq.substring(0, refPoint)
        this.seq = Util.translateLinear(seq.substring(refPoint + 1))
    }
}
