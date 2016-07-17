/*
 * Copyright 2016 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjdb.scoring

import com.milaboratory.core.tree.TreeSearchParameters

class SequenceSearcherPreset {
    final AlignmentScoring scoring
    final TreeSearchParameters parameters

    static final List<String> ALLOWED_PRESETS = ["hamming", "high-recall", "optimal", "high-precision"]

    static SequenceSearcherPreset byName(String name) {
        switch (name.toLowerCase()) {
            case ALLOWED_PRESETS[0]:
                return new SequenceSearcherPreset(DummyAlignmentScoring.INSTANCE,
                        new TreeSearchParameters(2, 1, 1, 2))
            case ALLOWED_PRESETS[1]:
                return byRecall(0.8f)
            case ALLOWED_PRESETS[2]:
                return getOptimal()
            case ALLOWED_PRESETS[3]:
                return byPrecision(0.8f)
            default:
                throw new RuntimeException("Unknown parameter preset '$name'")
        }
    }

    static SequenceSearcherPreset getOptimal() {
        def scoringMetadata = new ScoringMetadataTable().optimal

        System.err.println "[SearchPresetProvider] Requested scoring with highest F-score, " +
                "closest scoring schema selected is {$scoringMetadata}"

        byMetadata(scoringMetadata)
    }

    static SequenceSearcherPreset byPrecision(float precision) {
        def scoringMetadata = new ScoringMetadataTable().getByPrecision(precision)

        System.err.println "[SearchPresetProvider] Requested scoring with precision=$precision, " +
                "closest scoring schema selected is {$scoringMetadata}"

        byMetadata(scoringMetadata)
    }

    static SequenceSearcherPreset byRecall(float recall) {
        def scoringMetadata = new ScoringMetadataTable().getByRecall(recall)

        System.err.println "[SearchPresetProvider] Requested scoring with recall=$recall, " +
                "closest scoring schema selected is {$scoringMetadata}"

        byMetadata(scoringMetadata)
    }

    static SequenceSearcherPreset byMetadata(ScoringMetadata scoringMetadata) {
        new SequenceSearcherPreset(AlignmentScoringProvider.loadScoring(scoringMetadata.scoringId),
                new TreeSearchParameters(5, 2, 2, 7))
    }

    SequenceSearcherPreset(AlignmentScoring scoring, TreeSearchParameters parameters) {
        this.scoring = scoring
        this.parameters = parameters
    }

    SequenceSearcherPreset withScoringThreshold(float scoringThreshold) {
        new SequenceSearcherPreset(scoring.withScoreThreshold(scoringThreshold),
                parameters);
    }

    SequenceSearcherPreset withSearchParameters(int maxSubstitutions,
                                                int maxInsertions,
                                                int maxDeletions,
                                                int maxMismatches) {
        new SequenceSearcherPreset(scoring,
                new TreeSearchParameters(maxSubstitutions, maxDeletions, maxInsertions, maxMismatches));
    }
}
