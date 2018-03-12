package com.antigenomics.vdjdb.impl.weights;

import com.antigenomics.vdjdb.db.Column;
import com.antigenomics.vdjdb.db.Row;
import com.antigenomics.vdjdb.impl.ClonotypeDatabase;
import com.antigenomics.vdjdb.sequence.SearchScope;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.SequenceTreeMap;

import java.util.*;

public class DegreeWeightFunctionFactory implements WeightFunctionFactory {
    // todo: to glossary / pre-set in ClonotypeDatabase
    public static final DegreeWeightFunctionFactory DEFAULT = new DegreeWeightFunctionFactory("antigen.epitope");

    private final String groupingColumnName;

    public DegreeWeightFunctionFactory(String groupingColumnName) {
        this.groupingColumnName = groupingColumnName;
    }

    public String getGroupingColumnName() {
        return groupingColumnName;
    }

    @Override
    public WeightFunction create(ClonotypeDatabase clonotypeDatabase) {
        int groupColumnIndex = -1, i = 0;
        for (Column column : clonotypeDatabase.getColumns()) {
            i++;
            if (column.getName().equals(groupingColumnName)) {
                groupColumnIndex = i;
                break;
            }
        }

        if (groupColumnIndex == -1) {
            throw new IllegalArgumentException("Clonotype database doesn't have specified grouping column '" +
                    groupingColumnName + "'");
        }

        // create sequence tree map holding epitope info & aux. set of cdrs
        // todo: here we don't handle/overwrite cross-reactive clonotypes
        SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);
        Set<Cdr3Info> cdr3InfoSet = new HashSet<>();
        for (Row row : clonotypeDatabase.getRows()) {
            String cdr3aa = row.getAt(clonotypeDatabase.getCdr3ColIdx()).getValue(),
                    group = row.getAt(groupColumnIndex).getValue();

            Cdr3Info cdr3Info = new Cdr3Info(new AminoAcidSequence(cdr3aa), group);
            stm.put(cdr3Info.cdr3, cdr3Info);
            cdr3InfoSet.add(cdr3Info);
        }

        // perform search for every cdr3 with a given scope, compute scores
        Map<String, Float> weights = new HashMap<>();
        cdr3InfoSet
                .parallelStream()
                .map(x -> new AbstractMap.SimpleEntry<>(x, computeMatches(clonotypeDatabase.getSearchScope(),
                        stm, x)))
                .forEach(x -> weights.put(x.getKey().cdr3.toString(), x.getValue()));


        return new DegreeWeightFunction(weights);
    }

    private static float computeMatches(SearchScope searchScope,
                                        SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm,
                                        Cdr3Info query) {
        return 1.0f;
    }

    private class Cdr3Info {
        final AminoAcidSequence cdr3;
        final String group;

        Cdr3Info(AminoAcidSequence cdr3, String group) {
            this.cdr3 = cdr3;
            this.group = group;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Cdr3Info cdr3Info = (Cdr3Info) o;

            if (!cdr3.equals(cdr3Info.cdr3)) return false;
            return group.equals(cdr3Info.group);
        }

        @Override
        public int hashCode() {
            int result = cdr3.hashCode();
            result = 31 * result + group.hashCode();
            return result;
        }
    }
}
