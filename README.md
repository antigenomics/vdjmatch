## VDJdb: a software for functional annotation of T-cell repertoires

In a nutshell, VDJdb combines an API and simple CLI software implementation for browsing a curated databased of annotated V(D)J junctions. The software accepts output of commonly used immune repertoire sequencing processing tools using [VDJtools](http://vdjtools-doc.readthedocs.org/en/latest/index.html) engine. Therefore all input data should be converted to VDJtools [format](http://vdjtools-doc.readthedocs.org/en/latest/input.html#vdjtools-format). Note that the software also accepts [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) symantics introduced by VDJtools.

VDJdb is distributed as a JAR executable (see [releases section]()) and can be incorporated into projects using [Maven](https://maven.apache.org/). The software is cross-platform and requires [JRE 1.8]() to run.

### Running VDJdb

Standalone VDJdb annotation utility can be executed by running

```
java -jar -Xmx4G vdjdb.jar vdjdb [options] [sample1 sample2 sample3 ... if -m is not specified] output_prefix
```

First part of the command runs the JAR file and sets the memory limit to 4GB. The second part includes options, input samples and prefix of output files.
The list of the options is the following:

| Shorthand | Long name       | Required | Argument example                               |  Description                                                |
|-----------|-----------------|----------|------------------------------------------------|-------------------------------------------------------------|
| -h        |                 |          |                                                |  Display help message                                       |
| -m        | --metadata      |          | metadata.txt                                   |  A [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) file, holding paths to samples and user-provided information.
| -S        | --species       | Yes      | human,mouse,...                                |  Name of the species. All samples should belong to the same species, only one species is allowed.
| -R        | --chain         | Yes      | TRA,TRB,...                                    |  Name of the receptor chain. All samples should contain to the same receptor chain, only one chain is allowed.
| -v        | --v-match       | Yes      |                                                |  Require Variable segment matching when searching the database
| -j        | --j-match       | Yes      |                                                |  Require Joining segment matching when searching the database
|           | --summary       |          | origin,disease.type,disease,source             |  A comma-separated list of database column names on which summary statistics will be computed
|           | --search-params |          | 2,1,1,2                                        |  CDR3 sequence search parameters: allowed number of substitutions (s), insertions (i), deletions (d) and total number of mutations.
|           | --database      |          | vdjdb_v10                                      |  Path and prefix of an external database.
| -c        | --compress      |          |                                                |  Compress sample-level summary output with GZIP.

The following output files will be generated:

-
-
-

### Some API examples

Coming soon... also check [VDJdb-server]()


