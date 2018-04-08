package com.antigenomics.vdjdb.cli

class OptCluster extends OptBase {
    final boolean optIsoMds
    final int optIsoMdsD, optIsoMdsMinCompSz

    OptCluster(CliCluster cliBase, String[] args) {
        super(cliBase, args)

        this.optIsoMds = (boolean) opt.'isomds'
        this.optIsoMdsD = (opt.'isomds-d' ?: CliCluster.DEFAULT_ISOMDS_D).toInteger()
        this.optIsoMdsMinCompSz = (opt.'isomds-mincompsz' ?: CliCluster.DEFAULT_ISOMDS_MINCOMPSZ).toInteger()
    }
}
