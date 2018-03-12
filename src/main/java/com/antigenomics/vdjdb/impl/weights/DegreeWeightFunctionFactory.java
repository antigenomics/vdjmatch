package com.antigenomics.vdjdb.impl.weights;

import com.antigenomics.vdjdb.db.Column;
import com.antigenomics.vdjdb.db.Entry;
import com.antigenomics.vdjdb.db.Row;
import com.antigenomics.vdjdb.impl.ClonotypeDatabase;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.SequenceTreeMap;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

        // create sequence tree map holding epitope info
        SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);
        Set<AminoAcidSequence> cdr3aaSet = new HashSet<>();
        for (Row row : clonotypeDatabase.getRows()) {
            // todo: here we don't handle/overwrite cross-reactive clonotypes
            String cdr3aa = row.getAt(clonotypeDatabase.getCdr3ColIdx()).getValue(),
                    group = row.getAt(groupColumnIndex).getValue();

            AminoAcidSequence cdr3aaBin = new AminoAcidSequence(cdr3aa);
            cdr3aaSet.add(cdr3aaBin);
            stm.put(cdr3aaBin, new Cdr3Info(cdr3aaBin, group));
        }

        // perform search for every cdr3 with a given scope
        cdr3aaSet.parallelStream().

        // accumulate hit count,

        return new DegreeWeightFunction(new HashMap<>());
    }

    private class Cdr3Info {
        final AminoAcidSequence cdr3;
        final String group;

        public Cdr3Info(AminoAcidSequence cdr3, String group) {
            this.cdr3 = cdr3;
            this.group = group;
        }
    }
}
