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
    private final LinearGapAlignmentScoring scoring;
    private final float[] positionWeights;
    private final int posCenterBin;
    private float scoreThreshold;

    public VdjdbAlignmentScoring(LinearGapAlignmentScoring scoring,
                                 float[] positionWeights,
                                 float scoreThreshold) {
        this.scoring = scoring;
        this.positionWeights = positionWeights;
        this.posCenterBin = positionWeights.length / 2;
        this.scoreThreshold = scoreThreshold;
    }

    private int getBin(int centeredPos) {
        int k = centeredPos + posCenterBin;
        return k < 0 ? 0 : (k < positionWeights.length ? k : (positionWeights.length - 1));
    }

    private float getPositionWeight(int pos, int cdr3Length) {
        int center = cdr3Length / 2;
        if (cdr3Length % 2 == 0) {
            return 0.5f * (positionWeights[getBin(pos - center)] + positionWeights[getBin(pos - center + 1)]);
        } else {
            return positionWeights[getBin(pos - center)];
        }
    }

    @Override
    public float computeBaseScore(Sequence reference) {
        float score = 0;
        for (int i = 0; i < reference.size(); ++i) {
            byte aa = reference.codeAt(i);
            score += scoring.getScore(aa, aa) * getPositionWeight(i, reference.size());
        }
        return score;
    }

    public float computeScore(Alignment alignment) {
        Sequence reference = alignment.getSequence1();
        return computeScore(alignment.getAbsoluteMutations(), computeBaseScore(reference), reference.size());
    }

    @Override
    public float computeScore(Mutations mutations, float baseScore, int refLength) {
        float score = baseScore;

        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i);

            double deltaScore = 0;
            if (isInsertion(mutation)) {
                deltaScore += scoring.getGapPenalty();
            } else {
                byte from = getFrom(mutation);
                deltaScore += isDeletion(mutation) ? scoring.getGapPenalty() :
                        (scoring.getScore(from, getTo(mutation)));
                deltaScore -= scoring.getScore(from, from);
            }

            score += deltaScore * getPositionWeight(getPosition(mutation), refLength);
        }

        return score;
    }

    public LinearGapAlignmentScoring getScoring() {
        return scoring;
    }

    public float[] getPositionWeights() {
        return positionWeights;
    }

    public int getPosCenterBin() {
        return posCenterBin;
    }

    @Override
    public float getScoreThreshold() {
        return scoreThreshold;
    }

    @Override
    public AlignmentScoring withScoreThreshold(float scoreThreshold) {
        return new VdjdbAlignmentScoring(scoring,
                positionWeights, scoreThreshold);
    }
}
