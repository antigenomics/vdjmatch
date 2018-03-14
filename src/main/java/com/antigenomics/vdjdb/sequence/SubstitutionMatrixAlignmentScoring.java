package com.antigenomics.vdjdb.sequence;

import com.milaboratory.core.alignment.BLASTMatrix;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutation;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class SubstitutionMatrixAlignmentScoring implements AlignmentScoring {
    private static final int N = AminoAcidSequence.ALPHABET.size();
    private final float[][] substitutionPenalties = new float[N][N];
    private final float[] gapPenalties = new float[N];
    private final float gapFactor;

    public static final SubstitutionMatrixAlignmentScoring DEFAULT_BLOSUM62;

    static {
        float[][] sm = new float[N][N];
        LinearGapAlignmentScoring<AminoAcidSequence> s = LinearGapAlignmentScoring.getAminoAcidBLASTScoring(BLASTMatrix.BLOSUM62);
        int maxScore = 0;
        for (byte i = 0; i < N; i++) {
            for (byte j = 0; j < N; j++) {
                int score = s.getScore(i, j);
                sm[i][j] = score;
                maxScore = Math.max(maxScore, Math.abs(score));
            }
        }
        DEFAULT_BLOSUM62 = new SubstitutionMatrixAlignmentScoring(sm, (float) (maxScore + 1));
    }

    public SubstitutionMatrixAlignmentScoring(float[][] substitutionMatrix, float gapFactor) {
        for (int i = 0; i < N; i++) {
            gapPenalties[i] = substitutionMatrix[i][i];
            for (int j = 0; j < N; j++) {
                substitutionPenalties[i][j] = substitutionMatrix[i][j] -
                        Math.max(substitutionMatrix[i][i], substitutionMatrix[j][j]);
            }
        }
        this.gapFactor = gapFactor;
    }

    @Override
    public float computeScore(AminoAcidSequence query, Mutations<AminoAcidSequence> mutations) {
        float score = 0;
        int indels = 0;

        for (int i = 0; i < mutations.size(); i++) {
            int code = mutations.getMutation(i);

            if (Mutation.isSubstitution(code)) {
                score += substitutionPenalties[Mutation.getFrom(code)][Mutation.getTo(code)];
            } else {
                indels++;
                if (Mutation.isDeletion(code)) {
                    score -= gapPenalties[Mutation.getFrom(code)];
                } else {
                    score -= gapPenalties[Mutation.getTo(code)];
                }
            }
        }

        return score - indels * gapFactor;
    }
}
