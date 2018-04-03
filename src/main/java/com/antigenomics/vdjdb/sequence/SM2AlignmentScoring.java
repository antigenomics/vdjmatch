package com.antigenomics.vdjdb.sequence;

import com.milaboratory.core.alignment.BLASTMatrix;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutation;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class SM2AlignmentScoring implements AlignmentScoring {
    private static final int N = AminoAcidSequence.ALPHABET.size();
    private final float[][] substitutionMatrix;
    private final float gapFactor;

    public static final SM2AlignmentScoring DEFAULT_BLOSUM62;

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
        DEFAULT_BLOSUM62 = new SM2AlignmentScoring(sm, -(float) (maxScore + 1));
    }

    public SM2AlignmentScoring(float[][] substitutionMatrix, float gapFactor) {
        this.substitutionMatrix = substitutionMatrix;
        this.gapFactor = gapFactor;
        assert gapFactor <= 0;
    }

    @Override
    public float computeScore(AminoAcidSequence query, Mutations<AminoAcidSequence> mutations) {
        float queryScore = 0;
        for (int i = 0; i < query.size(); i++) {
            byte base = query.codeAt(i);
            queryScore += substitutionMatrix[base][base];
        }
        float targetScore = queryScore, score = queryScore;
        int indels = 0;

        for (int i = 0; i < mutations.size(); i++) {
            int code = mutations.getMutation(i);

            if (Mutation.isSubstitution(code)) {
                byte from = Mutation.getFrom(code),
                        to = Mutation.getTo(code);
                float matchScore = substitutionMatrix[from][from];
                // replace match by mismatch
                score -= matchScore;
                score += substitutionMatrix[from][to];
                // target has different base:
                // remove score for query base add score for this base
                targetScore -= matchScore;
                targetScore += substitutionMatrix[to][to];
            } else {
                indels++;
                if (Mutation.isDeletion(code)) {
                    // Query has additional base not found in target
                    byte from = Mutation.getFrom(code);
                    float matchScore = substitutionMatrix[from][from];
                    // remove base match score from alignment
                    score -= matchScore;
                    // also remove base score from target
                    targetScore -= matchScore;
                } else {
                    // Query has missing base found in target
                    byte to = Mutation.getTo(code);
                    float matchScore = substitutionMatrix[to][to];
                    // do not remove base match score from alignment

                    // add score to target
                    targetScore += matchScore;
                }
            }
        }

        //System.out.println(score);
        //System.out.println(queryScore);
        //System.out.println(targetScore);
        //System.out.println(indels);

        return score - Math.max(queryScore, targetScore) + gapFactor * indels;
    }

    @Override
    public ScoringType getScoringType() {
        return ScoringType.Probabilistic;
    }

    public float[][] getSubstitutionMatrix() {
        return substitutionMatrix;
    }

    public float getGapFactor() {
        return gapFactor;
    }
}
