package com.antigenomics.vdjdb.cli

class CliCluster extends CliBase<OptCluster> {
    static def DEFAULT_ISOMDS_MINCOMPSZ = "5", DEFAULT_ISOMDS_D = "2"

    CliCluster() {
        super("cluster")

        cli._(longOpt: "isomds",
                "Run isoMDS algorithm to embed clonotypes")
        cli._(longOpt: "isomds-mincompsz", argName: "integer", args: 1,
                "Minimal connected component size for isoMDS algorithm. [default = $DEFAULT_ISOMDS_MINCOMPSZ]")
        cli._(longOpt: "isomds-d", argName: "integer", args: 1,
                "Number of dimensions for isoMDS algorithm. [default = $DEFAULT_ISOMDS_D]")
    }

    @Override
    OptCluster parseArguments(String[] args) {
        new OptCluster(this, args)
    }
}
