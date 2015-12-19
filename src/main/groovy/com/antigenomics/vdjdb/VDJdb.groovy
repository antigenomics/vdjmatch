package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.db.Database
import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjdb.sequence.SequenceFilter
import com.antigenomics.vdjdb.text.TextFilter

class VDJdb {
    static final String DEFAULT_DB_RESOURCE_NAME = "vdjdb_legacy.txt",
                        DEFAULT_META_RESOURCE_NAME = "vdjdb_legacy.meta"

    private static final Database dbInstance = new Database(Util.resourceAsStream(DEFAULT_META_RESOURCE_NAME))

    static {
        dbInstance.addEntries(Util.resourceAsStream(DEFAULT_DB_RESOURCE_NAME))
    }

    static Database get(List<TextFilter> textFilters = [],
                        List<SequenceFilter> sequenceFilters = []) {
        Database.create(dbInstance.search(textFilters, sequenceFilters), dbInstance)
    }

    static ClonotypeDatabase asClonotypeDatabase(Database db, boolean matchV = false, boolean matchJ = false,
                                                 int maxMismatches = 2, int maxInsertions = 1,
                                                 int maxDeletions = 1, int maxMutations = 2, int depth = -1) {
        new ClonotypeDatabase(db.columns, matchV, matchJ,
                maxMismatches, maxInsertions,
                maxDeletions, maxMutations, depth)
    }
}