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

import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.Sequence;

import static com.milaboratory.core.mutations.Mutation.*;
import static com.milaboratory.core.mutations.Mutation.getPosition;

public class VdjdbAlignmentScoring implements AlignmentScoring {
    private final float[][] substitutionMatrix;
    private final float[] gapPenalties;

    public VdjdbAlignmentScoring(float[][] substitutionMatrix,
                                 float[] gapPenalties) {
        this.substitutionMatrix = substitutionMatrix;
        this.gapPenalties = gapPenalties;
    }

    @Override
    public float computeBaseScore(Sequence reference) {
        float score = 0;
        for (int i = 1; i < reference.size() - 1; ++i) { // exclude C and F/W
            byte aa = reference.codeAt(i);
            score += substitutionMatrix[aa][aa];
        }
        return score;
    }

    @Override
    public float computeScore(Alignment alignment) {
        Sequence reference = alignment.getSequence1();
        return computeScore(alignment.getAbsoluteMutations(), computeBaseScore(reference), reference.size());
    }

    @Override
    public float computeScore(Mutations mutations, float baseScore, int refLength) {
        float score = baseScore;

        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i);
            int pos = getPosition(mutation);
            if (pos > 0 & pos < refLength -1) { // exclude 1st and last AAs of CDR3, C and F/W
                double deltaScore = 0;          // This helps with TRAJ having W/F at the end

                if (isInsertion(mutation)) {
                    byte to = getTo(mutation);
                    deltaScore += gapPenalties[getTo(mutation)];
                    deltaScore -= substitutionMatrix[to][to];
                } else {
                    byte from = getFrom(mutation);
                    deltaScore += isDeletion(mutation) ? gapPenalties[from] :
                            substitutionMatrix[from][getTo(mutation)];
                    deltaScore -= substitutionMatrix[from][from];
                }

                score += deltaScore;
            }
        }

        return score;
    }

    @Override
    public float computePValue(float score) {
        return 1f;
    }
}
