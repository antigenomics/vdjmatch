package com.antigenomics.vdjdb.sequence;

import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.Sequence;

import static com.milaboratory.core.mutations.Mutation.*;
import static com.milaboratory.core.mutations.Mutation.getPosition;

public class AlignmentScoring {
    private final LinearGapAlignmentScoring scoring;
    private final double[] positionWeights;
    private final int posCenterBin;
    private double scoreThreshold;

    public AlignmentScoring(LinearGapAlignmentScoring scoring,
                            double[] positionWeights,
                            double scoreThreshold) {
        this.scoring = scoring;
        this.positionWeights = positionWeights;
        this.posCenterBin = positionWeights.length / 2;
        this.scoreThreshold = scoreThreshold;
    }

    private int getBin(int centeredPos) {
        int k = centeredPos + posCenterBin;
        return k < 0 ? 0 : (k < positionWeights.length ? k : (positionWeights.length - 1));
    }

    private double getPositionWeight(int pos, int cdr3Length) {
        int center = cdr3Length / 2;
        if (cdr3Length % 2 == 0) {
            return 0.5 * positionWeights[getBin(pos - center)] + 0.5 * positionWeights[getBin(pos - center + 1)];
        } else {
            return positionWeights[getBin(pos - center)];
        }
    }

    public double computeScore(Alignment alignment) {
        Sequence reference = alignment.getSequence1();
        Mutations mutations = alignment.getAbsoluteMutations();
        double score = 0;

        for (int i = 0; i < reference.size(); i++) {
            byte aa = reference.codeAt(i);
            score += scoring.getScore(aa, aa) * getPositionWeight(i, reference.size());
        }

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

            score += deltaScore * getPositionWeight(getPosition(mutation), reference.size());
        }

        return score;
    }

    public LinearGapAlignmentScoring getScoring() {
        return scoring;
    }

    public double[] getPositionWeights() {
        return positionWeights;
    }

    public int getPosCenterBin() {
        return posCenterBin;
    }

    public double getScoreThreshold() {
        return scoreThreshold;
    }

    public void setScoreThreshold(double scoreThreshold) {
        this.scoreThreshold = scoreThreshold;
    }
}
