/*
 * Copyright 2015 Mikhail Shugay
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

    static final List<String> ALLOWED_PRESETS = ["dummy", "high-recall", "balanced", "high-precision"]

    static SequenceSearcherPreset byName(String name) {
        switch (name.toLowerCase()) {
            case "dummy":
                return new SequenceSearcherPreset(DummyAlignmentScoring.INSTANCE,
                        new TreeSearchParameters(2, 1, 1, 2))
            case "high-recall":
            case "balanced":
            case "high-precision":
                return new SequenceSearcherPreset(AlignmentScoringProvider.loadScoring(name.toLowerCase()),
                        new TreeSearchParameters(5, 2, 2, 7))
            default:
                throw new RuntimeException("Unknown parameter preset '$name'")
        }
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
