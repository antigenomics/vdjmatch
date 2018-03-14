package com.antigenomics.vdjdb.sequence;

import com.milaboratory.core.Range;
import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.alignment.AlignmentHelper;
import com.milaboratory.core.alignment.AlignmentUtils;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.mutations.MutationsUtil;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class Hit {
    private final AminoAcidSequence query;
    private final int matchSequenceLength;
    private final Mutations<AminoAcidSequence> mutations;
    private final float alignmentScore;

    public Hit(AminoAcidSequence query,
               Mutations<AminoAcidSequence> mutations,
               int matchSequenceLength,
               float alignmentScore) {
        this.query = query;
        this.matchSequenceLength = matchSequenceLength;
        this.mutations = mutations;
        this.alignmentScore = alignmentScore;
    }

    public Alignment<AminoAcidSequence> computeAlignment() {
        return new Alignment<>(query, mutations,
                new Range(0, query.size()),
                new Range(0, matchSequenceLength),
                alignmentScore);
    }

    public AminoAcidSequence getQuery() {
        return query;
    }

    public int getMatchSequenceLength() {
        return matchSequenceLength;
    }

    public float getAlignmentScore() {
        return alignmentScore;
    }

    public Mutations<AminoAcidSequence> getMutations() {
        return mutations;
    }

    @Override
    public String toString() {
        return "<hit>\n" + AlignmentUtils.toStringSimple(query, mutations) +
                "s=" + alignmentScore + "\n" +
                "</hit>";
    }
}
