package com.antigenomics.vdjdb.impl

import com.antigenomics.vdjdb.impl.model.AggregateScoring
import com.antigenomics.vdjdb.impl.model.DummyAggregateScoring
import com.antigenomics.vdjdb.impl.segment.DummySegmentScoring
import com.antigenomics.vdjdb.impl.segment.SegmentScoring
import com.antigenomics.vdjdb.sequence.AlignmentScoring
import com.antigenomics.vdjdb.sequence.DummyAlignmentScoring

class ScoringBundle {
    final AggregateScoring aggregateScoring
    final SegmentScoring segmentScoring
    final AlignmentScoring alignmentScoring

    ScoringBundle(AlignmentScoring alignmentScoring = DummyAlignmentScoring.INSTANCE,
                  SegmentScoring segmentScoring = DummySegmentScoring.INSTANCE,
                  AggregateScoring aggregateScoring = DummyAggregateScoring.INSTANCE) {
        this.aggregateScoring = aggregateScoring
        this.segmentScoring = segmentScoring
        this.alignmentScoring = alignmentScoring
    }
}
