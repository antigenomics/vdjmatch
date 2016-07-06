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

package com.antigenomics.vdjdb.scoring;

import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.Sequence;

public class DummyAlignmentScoring implements AlignmentScoring {
    public static DummyAlignmentScoring INSTANCE = new DummyAlignmentScoring();

    private DummyAlignmentScoring() {

    }

    @Override
    public float computeScore(Mutations mutations, float baseScore, int refLength) {
        return refLength - mutations.size();
    }

    @Override
    public float computeBaseScore(Sequence reference) {
        return 0f;
    }

    @Override
    public float getScoreThreshold() {
        return Float.MIN_VALUE;
    }

    @Override
    public AlignmentScoring withScoreThreshold(float scoreThreshold) {
        return this;
    }

}
