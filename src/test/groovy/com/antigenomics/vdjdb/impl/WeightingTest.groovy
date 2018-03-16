package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.impl.weights.DegreeWeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class WeightingTest {
    @Test
    void degreeWeightFunctionTest() {
        def dwff = new DegreeWeightFunctionFactory("ag")

        def createInfo = { String cdr3, String ag ->
            new DegreeWeightFunctionFactory.Cdr3Info(new AminoAcidSequence(cdr3), ag)
        }

        def cdr3InfoSet = new HashSet<DegreeWeightFunctionFactory.Cdr3Info>(
                [
                        createInfo("AAAAAAA", "1"),
                        createInfo("AAAAAWA", "1"),
                        createInfo("AAWAAAA", "2"),  // 1
                        createInfo("AAAWAAA", "2"),  // 2
                        createInfo("AAWWAAA", "2"),  // 3
                        createInfo("AAWWWAA", "2"),
                        createInfo("AAWWAA", "2"),   // 4
                        createInfo("AAAAW", "2")
                ]
        )

        def weighting = dwff.create(cdr3InfoSet, new SearchScope(2, 1, 3))

        assert (float) (Math.exp(-weighting.computeWeight("", "", "AAAAAAA")) - 1) == 4f
    }
}
