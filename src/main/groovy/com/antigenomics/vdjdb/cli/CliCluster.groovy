package com.antigenomics.vdjdb.cli

class CliCluster extends CliBase<OptCluster> {
    CliCluster() {
        super("cluster")
    }

    @Override
    OptCluster parseArguments(String[] args) {
        new OptCluster(this, args)
    }
}
