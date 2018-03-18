package com.antigenomics.vdjdb.impl.weights;

import com.antigenomics.vdjdb.db.Column;
import com.antigenomics.vdjdb.db.Row;
import com.antigenomics.vdjdb.impl.ClonotypeDatabase;
import com.antigenomics.vdjdb.sequence.SearchScope;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;

import java.util.*;

public class DegreeWeightFunctionFactory implements WeightFunctionFactory {
    // todo: to glossary / pre-set in ClonotypeDatabase
    public static final DegreeWeightFunctionFactory DEFAULT = new DegreeWeightFunctionFactory(ClonotypeDatabase.getEPITOPE_COL_DEFAULT());

    private final String groupingColumnName;

    public DegreeWeightFunctionFactory(String groupingColumnName) {
        this.groupingColumnName = groupingColumnName;
    }

    public String getGroupingColumnName() {
        return groupingColumnName;
    }

    public WeightFunction create(final Set<Cdr3Info> cdr3InfoSet, final SearchScope searchScope) {
        SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);
        for (Cdr3Info cdr3Info : cdr3InfoSet) {
            stm.put(cdr3Info.cdr3, cdr3Info);
        }

        // perform search for every cdr3 with a given scope, compute scores
        Map<String, Float> weights = new HashMap<>();
        cdr3InfoSet
                .parallelStream()
                .map(x -> new AbstractMap.SimpleEntry<>(x, computeMatches(searchScope,
                        stm, x)))
                .forEach(x -> weights.put(x.getKey().cdr3.toString(), x.getValue()));


        return new DegreeWeightFunction(weights);
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
        Set<Cdr3Info> cdr3InfoSet = new HashSet<>();
        for (Row row : clonotypeDatabase.getRows()) {
            String cdr3aa = row.getAt(clonotypeDatabase.getCdr3ColIdx()).getValue(),
                    group = row.getAt(groupColumnIndex).getValue();

            cdr3InfoSet.add(new Cdr3Info(new AminoAcidSequence(cdr3aa), group));
        }


        return create(cdr3InfoSet, clonotypeDatabase.getSearchScope());
    }

    private static float computeMatches(SearchScope searchScope,
                                        SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm,
                                        Cdr3Info query) {
        Set<AminoAcidSequence> matches = new HashSet<>();

        NeighborhoodIterator<AminoAcidSequence, Cdr3Info> ni = stm.getNeighborhoodIterator(query.cdr3,
                searchScope.getTreeSearchParameters());

        Cdr3Info target;
        while ((target = ni.next()) != null) {
            if (!Objects.equals(target.group, query.group) && // different groups
                    Math.abs(target.cdr3.size() - query.cdr3.size()) <= searchScope.getMaxIndels()) { // indel sum threshold
                matches.add(target.cdr3);
            }
        }

        return -(float) Math.log1p(matches.size());
    }

    public static class Cdr3Info {
        private final AminoAcidSequence cdr3;
        private final String group;

        public Cdr3Info(AminoAcidSequence cdr3, String group) {
            this.cdr3 = cdr3;
            this.group = group;
        }

        public AminoAcidSequence getCdr3() {
            return cdr3;
        }

        public String getGroup() {
            return group;
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
