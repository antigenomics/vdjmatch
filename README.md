[![Build Status](https://travis-ci.org/antigenomics/vdjdb-standalone.svg?branch=master)](https://travis-ci.org/antigenomics/vdjdb-standalone)

## VDJdb-standalone: a software for functional annotation of T-cell repertoires

In a nutshell, VDJdb-standalone is an API package with a command-line interface implementation allowing to query T-cell receptor (TCR) repertoires against a database of TCRs with known antigen specificity. The database itself is hosted in a separate [repository](https://github.com/antigenomics/vdjdb-db) together with dedicated database build and curation utilities. VDJdb-standalone will automatically download and use the latest version of the database, it is also possible to use custom user-provided databases.

The software accepts TCR clonotype table(s) as an input and relies on [VDJtools](http://vdjtools-doc.readthedocs.org/en/latest/index.html) framework to parse the output of commonly used immune repertoire sequencing processing tools. Therefore all input data should be converted to VDJtools [format](http://vdjtools-doc.readthedocs.org/en/latest/input.html#vdjtools-format). Note that the software also accepts [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) symantics introduced by VDJtools to facilitate running for multi-sample datasets.

### Installation

VDJdb is distributed as a JAR executable (see [releases section](https://github.com/antigenomics/vdjdb/releases)), the software is cross-platform and requires [JRE 1.8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) to run. Running executable jar is quite straightforward, just use ``java -jar vdjdb-standalone.jar`` command as described below.

To compile VDJdb-standalone from source: 
* Install VDJtools of appropriate version (see ``build.gradle`` in repository root folder) using Maven (``mvn clean install``)
* Build VDJdb-standalone with Gradle by navigating to the root folder and running ``gradle clean build``.

### Running VDJdb

Standalone VDJdb annotation utility can be executed by running the following command:

```
java -jar -Xmx4G vdjdb.jar vdjdb [options] [sample1 sample2 sample3 ... if -m is not specified] output_prefix
```

First part of the command runs the JAR file and sets the memory limit to 4GB (should be increased in case JVM drops with heap size exception). The second part includes options, input samples and prefix of output files.

The list of accepted options is the following:

| Shorthand | Long  name          | Required | Argument example |  Description |
|-----------|---------------------|----------|------------------|--------------|
| ``-h``    |                     |          |                  |  Display help message |
| ``-m``    | ``--metadata``      |          | metadata.txt     |  A [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) file, holding paths to samples and user-provided information. |
| ``-S``    | ``--species``       | Yes      | human,mouse,...  |  Species name. All samples should belong to the same species, only one species is allowed. |
| ``-R``    | ``--gene``          | Yes      | TRA,TRB,...      |  Name of the receptor gene. All samples should contain to the same receptor gene, only one gene is allowed. |
| ``-v``    | ``--v-match``       | Yes      |                  |  Require Variable segment matching when searching the database. |
| ``-j``    | ``--j-match``       | Yes      |                  |  Require Joining segment matching when searching the database. |
|           | ``--search-params`` |          | 2,1,1,2          |  CDR3 sequence search parameters in s/i/d/m format: allowed number of substitutions (s), insertions (i), deletions (d) and the total number of mismatches (m). |
|           | ``--database``      |          | /path/to/my_db   |  Path and prefix of an external database. Should point to files with '.txt', and '.meta.txt' suffices (the database itself and database metadata).|
|           | ``--use-fat-db``    |          |                  |  In case running with a built-in database, will use full database version instead of slim one. Full database contains extended info on method used to identify a given specific TCR and sample source, but has a higher degree of redundancy (several identical TCR:pMHC pairs from different publications, etc) that can complicate post-analysis |
|           |  ``--filter``       |          | ``__antigen.species__ =~ "EBV" | Logical filter expresstion that will be evaluated for database columns. Supports Java/Groovy syntax, Regex, .contains(), .startsWith(), etc. Parst of the expression marked with double underscore (``__``) will be subsituted with corresponding values from database rows. Those parts should be named exactly as columns in the database |
|           | ``--vdjdb-confidence-threshold`` | | 2 | VDJdb confidence level threshold, should lie in ``[0,7]`` (default is 2). See [database readme](https://github.com/antigenomics/vdjdb-db) for details on VDJdb confidence scores |
| ``-c``    | ``--compress``      |          |                  |  Compress sample-level summary output with GZIP.                                                                                                                                         |

The following output files will be generated:

* ``annot.stats.txt`` annotation statistics containing the number (``counter.type:unique``) and frequency (``counter.type:weighted``) of clonotypes that were not found, found once, and found 2+ times in the database.
* ``annot.summary.txt`` annotation summary containing the number (rows with ``unique`` record) and frequency (rows with ``weighted`` record) of clonotypes that match certain antigen annotation. Note that each clonotype:antigen pair is counted once, even if multiple matches are present.
* ``$sample_id.annot.txt`` annotation for each of the clonotypes found in database. This is an all-to-all merge table. Clonotype information from the sample (count, frequency, cdr3 sequence, v/d/j segments and v/d/j markup) is preserved. For the information on database columns that are appended see database schema in [VDJdb-db repository](https://github.com/antigenomics/vdjdb-db) readme. The ``score`` and ``index`` columns represent CDR3 alignment score and clonotype index respectively. The latter is useful as a clonotype can match several CDR3 records of the database.

### See also

A web-based application VDJdb-server that runs on top of [VDJdb-server](https://github.com/antigenomics/vdjdb-server) can be used to browse the database and annotated samples.