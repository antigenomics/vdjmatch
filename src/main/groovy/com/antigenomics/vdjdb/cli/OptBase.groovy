package com.antigenomics.vdjdb.cli

import com.antigenomics.vdjdb.Util
import com.antigenomics.vdjdb.impl.ScoringBundle
import com.antigenomics.vdjdb.impl.ScoringProvider
import com.antigenomics.vdjdb.impl.filter.DummyResultFilter
import com.antigenomics.vdjdb.impl.filter.MaxScoreResultFilter
import com.antigenomics.vdjdb.impl.filter.ResultFilter
import com.antigenomics.vdjdb.impl.filter.ScoreThresholdResultFilter
import com.antigenomics.vdjdb.impl.filter.TopNResultFilter
import com.antigenomics.vdjdb.impl.weights.DegreeWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.DummyWeightFunctionFactory
import com.antigenomics.vdjdb.impl.weights.WeightFunctionFactory
import com.antigenomics.vdjdb.sequence.SearchScope
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.SampleCollection

class OptBase {
    final CliBase cliBase
    final OptionAccessor opt

    /* software type, metadata, input files & output options */
    final SampleCollection sampleCollection
    final boolean optCompress
    final String outputPrefix

    /* species, gene, pre-filtering */
    final String optSpecies, optGene

    /* initial search */
    final boolean optVMatch, optJMatch
    final SearchScope searchScope

    /* scoring */
    final ScoringBundle scoringBundle

    /* filtering */
    final ResultFilter resultFilter

    /* Weighting */
    final WeightFunctionFactory weightFunctionFactory

    OptBase(CliBase cliBase, String[] args) {
        this.cliBase = cliBase

        def cli = cliBase.cli

        Util.info("Executing {" + cliBase.routineName + "} with arguments {" + args.join(" ") + "}")
        this.opt = cli.parse(args)
        if (opt == null) {
            System.exit(1)
        }
        if (opt.h || opt.arguments().size() == 0) {
            cli.usage()
            System.exit(1)
        }

        /* software type, metadata, input files & output options */

        def optMetadataFileName = opt.m
        def optSoftware = opt.'software' ? Software.byName((String) opt.'software') : Software.VDJtools
        this.optCompress = (boolean) opt.c
        this.outputPrefix = opt.arguments()[-1]
        if (optMetadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
            if (optMetadataFileName) {
                Util.error("Only output prefix should be provided in case of -m")
            } else {
                Util.error("At least 1 sample files should be provided if not using -m")
            }
            cli.usage()
            System.exit(1)
        }

        // No actual loading here - lazy load
        if (optMetadataFileName) {
            this.sampleCollection = new SampleCollection((String) optMetadataFileName, optSoftware)
            Util.info(" Using metadata from: " + optMetadataFileName +
                    ", ${sampleCollection.size()} sample(s) to process.")
        } else {
            this.sampleCollection = new SampleCollection(opt.arguments()[0..-2], optSoftware)
            Util.info("Using sample files: " + opt.arguments()[0..-2] +
                    ", ${sampleCollection.size()} sample(s) to process.")
        }

        /* species, gene, pre-filtering */

        def optSpecies = (String) opt.S

        def allowedSpecies = [CliBase.ALLOWED_SPECIES_ALIAS.keySet(), CliBase.ALLOWED_SPECIES_ALIAS.values()].flatten()
        if (!allowedSpecies.any { optSpecies.equalsIgnoreCase((String) it) }) {
            Util.error("Wrong species name, use one of ${allowedSpecies} (case-insensitive)")
        }
        this.optSpecies = CliBase.ALLOWED_SPECIES_ALIAS[optSpecies.toLowerCase()] ?: optSpecies

        this.optGene = (String) opt.R
        if (!CliBase.ALLOWED_GENES.any { optGene.equalsIgnoreCase(it) }) {
            Util.error("Wrong gene name, use one of $CliBase.ALLOWED_GENES (case-insensitive)")
        }

        /* initial search */

        this.optVMatch = (boolean) opt.'v-match'
        this.optJMatch = (boolean) opt.'j-match'

        def optSearchScope = (opt.'search-scope' ?: CliBase.DEFAULT_SEARCH_SCOPE).split(",").collect { it.toInteger() },
            optExhaustive = (opt.'search-exhaustive' ?: CliBase.DEFAULT_EXHAUSTIVE).toInteger()

        if (!opt.'scoring-vdjmatch') {
            optExhaustive = 0 // has no effect in case no scoring is used
        }

        this.searchScope = optSearchScope.size() == 3 ?
                new SearchScope(optSearchScope[0], optSearchScope[1], optSearchScope[2],
                        optExhaustive > 0, optExhaustive < 2)
                :
                new SearchScope(optSearchScope[0], optSearchScope[1], optSearchScope[2], optSearchScope[3],
                        optExhaustive > 0, optExhaustive < 2)

        /* scoring */

        def optVdjmatchScoring = (boolean) opt.'scoring-vdjmatch',
            optScoringMode = (opt.'scoring-mode' ?: CliBase.DEFAULT_SCORING_MODE).toInteger()
        this.scoringBundle = optVdjmatchScoring ?
                ScoringProvider.loadScoringBundle(
                        this.optSpecies, this.optGene,
                        optScoringMode == 0) :
                ScoringBundle.DUMMY

        /* filtering */

        def optScoreThreshold = (opt.'hit-filter-score' ?: "-Infinity").toFloat()
        if (opt.'hit-filter-max') {
            this.resultFilter = new MaxScoreResultFilter(optScoreThreshold)
        } else if (opt.'hit-filter-topn') {
            this.resultFilter = new TopNResultFilter(optScoreThreshold,
                    (int) (opt.'hit-filter-topn').toInteger())
        } else if (opt.'hit-filter-score') {
            this.resultFilter = new ScoreThresholdResultFilter(optScoreThreshold)
        } else {
            this.resultFilter = DummyResultFilter.INSTANCE
        }

        /* weighting */

        def optWeightByInfo = opt.'hit-weight-inf'
        this.weightFunctionFactory = optWeightByInfo ?
                DegreeWeightFunctionFactory.DEFAULT :
                DummyWeightFunctionFactory.INSTANCE
    }
}
