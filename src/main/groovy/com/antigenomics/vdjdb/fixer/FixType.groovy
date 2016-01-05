package com.antigenomics.vdjdb.fixer

enum FixType {
    NoFixNeeded(false, true),
    FixAdd(true, true), FixTrim(true, true),
    FixReplace(true, true), FailedReplace(true, false),
    FailedBadSegment(false, false),
    FailedNoAlignment(true, false)

    final boolean fixAttempted, good

    FixType(boolean fixAttempted, boolean good) {
        this.fixAttempted = fixAttempted
        this.good = good
    }
}