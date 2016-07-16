package com.antigenomics.vdjdb.scoring

import org.junit.Test

class MetadataTest {
    @Test
    void metadataLoadTest() {
        new ScoringMetadataTable()
    }

    @Test
    void metadataFetchTest() {
        [0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f].each {
            SequenceSearcherPreset.byPrecision(it)
            SequenceSearcherPreset.byRecall(it)
        }

        SequenceSearcherPreset.optimal
    }

    @Test
    void presetTest() {
        SequenceSearcherPreset.ALLOWED_PRESETS.each {
            SequenceSearcherPreset.byName(it)
        }
    }
}
