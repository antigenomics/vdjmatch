package com.antigenomics.vdjdb.impl.segment;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class PrecomputedSegmentScoring implements SegmentScoring {
    private final Map<String, Map<String, float[]>> pairedVScores = new HashMap<>();
    private final Map<String, float[]> vMarginalAvgScores = new HashMap<>();
    private final float[] vAvgScore = new float[3];

    private final Map<String, Map<String, Float>> pairedJScores = new HashMap<>();
    private final Map<String, Float> jMarginalAvgScores = new HashMap<>();
    private final float jAvgScore;

    public PrecomputedSegmentScoring(Map<String, Map<String, float[]>> vScores,
                                     Map<String, Map<String, Float>> jScores) {
        // Re-format segment naming
        for (String v1 : vScores.keySet()) {
            Map<String, float[]> inner = new HashMap<>();
            for (Map.Entry<String, float[]> v2Entry : vScores.get(v1).entrySet()) {
                inner.put(reformatSegmentId(v2Entry.getKey()), Arrays.copyOf(v2Entry.getValue(), 3));
            }
            pairedVScores.put(reformatSegmentId(v1), inner);
        }
        for (String j1 : jScores.keySet()) {
            Map<String, Float> inner = new HashMap<>();
            for (Map.Entry<String, Float> j2Entry : jScores.get(j1).entrySet()) {
                inner.put(reformatSegmentId(j2Entry.getKey()), j2Entry.getValue());
            }
            pairedJScores.put(reformatSegmentId(j1), inner);
        }

        // V score summaries
        int n = 0;
        for (String v : pairedVScores.keySet()) {
            // marginal average, in case one v is unknown
            float[] marginalScore = new float[3];
            int k = 0;

            for (float[] scores : pairedVScores.get(v).values()) {
                for (int i = 0; i < 3; i++) {
                    marginalScore[i] += scores[i];
                }
                k++;
            }

            vMarginalAvgScores.put(v, marginalScore);

            // global sum
            for (int i = 0; i < 3; i++) {
                vAvgScore[i] += marginalScore[i];
                marginalScore[i] = marginalScore[i] / k;
            }
            n += k;
        }
        for (int i = 0; i < 3; i++) {
            // overall average, in case both v are unknown
            vAvgScore[i] = vAvgScore[i] / n;
        }

        // J score summaries
        n = 0;
        float jAvgScoreTmp = 0;
        for (String j : pairedJScores.keySet()) {
            // marginal average, in case one j is unknown
            float marginalScore = 0;
            int k = 0;

            for (float score : pairedJScores.get(j).values()) {
                marginalScore += score;
                k++;
            }

            jMarginalAvgScores.put(j, marginalScore / k);

            // overall average, in case both j are unknown
            jAvgScoreTmp += marginalScore;
            n += k;
        }
        jAvgScore = jAvgScoreTmp / n;
    }

    private static String reformatSegmentId(String id) {
        return id.split("[*,]")[0];
    }

    @Override
    public SegmentScores computeScores(String v1, String v2,
                                       String j1, String j2) {
        // Re-format segment names
        v1 = reformatSegmentId(v1);
        v2 = reformatSegmentId(v2);
        j1 = reformatSegmentId(j1);
        j2 = reformatSegmentId(j2);

        // Get V scores
        float[] vResult;

        // get pair scores
        Map<String, float[]> vPairScores = pairedVScores.get(v1);
        if (vPairScores == null) {
            // not found v1, try v2
            vResult = vMarginalAvgScores.get(v2);
            if (vResult == null) {
                // no marginal - get average
                vResult = Arrays.copyOf(vAvgScore, 3);
            } else {
                // return marginal
                vResult = Arrays.copyOf(vResult, 3);
            }
        } else {
            // try get pair score
            vResult = vPairScores.get(v2);
            if (vResult == null) {
                // not found - there is still marginal for v1
                vResult = Arrays.copyOf(vMarginalAvgScores.get(v1), 3);
            } else {
                // return pair score
                vResult = Arrays.copyOf(vResult, 3);
            }
        }

        // Get J score
        Float jResult;
        Map<String, Float> jPairScores = pairedJScores.get(j1);
        if (jPairScores == null) {
            // not found j1, try j2
            jResult = jMarginalAvgScores.get(j2);
            if (jResult == null) {
                // no marginal - get average
                jResult = jAvgScore;
            }
        } else {
            // try get pair score
            jResult = jPairScores.get(j2);
            if (jResult == null) {
                // not found - there is still marginal for v1
                jResult = jMarginalAvgScores.get(j1);
            }
        }

        return new SegmentScores(vResult[0], vResult[1], vResult[2],
                jResult);
    }
}
