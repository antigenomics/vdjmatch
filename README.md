[![Build Status](https://travis-ci.org/antigenomics/vdjdb-standalone.svg?branch=master)](https://travis-ci.org/antigenomics/vdjdb-standalone)

## VDJdb-standalone: a software for functional annotation of T-cell repertoires

In a nutshell, VDJdb-standalone is an API package with a command-line interface implementation allowing to query T-cell receptor (TCR) repertoires against a database of TCRs with known antigen specificity. The database itself is hosted in a separate [repository](https://github.com/antigenomics/vdjdb-db) together with dedicated database build and curation utilities. VDJdb-standalone will automatically download and use the latest version of the database, it is also possible to use custom user-provided databases.

The software accepts TCR clonotype table(s) as an input and relies on [VDJtools](http://vdjtools-doc.readthedocs.org/en/latest/index.html) framework to parse the output of commonly used immune repertoire sequencing processing tools. See [format](http://vdjtools-doc.readthedocs.org/en/latest/input.html) section of VDJtools docs for the list of supported formats. Note that the software also accepts [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) semantics introduced by VDJtools to facilitate running annotation for multi-sample datasets.

### Installation

VDJdb is distributed as a JAR executable (see [releases section](https://github.com/antigenomics/vdjdb/releases)), the software is cross-platform and requires [JRE 1.8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) to run. Running executable jar is quite straightforward, just use ``java -jar vdjdb.jar`` command as described below.

To compile VDJdb-standalone from source: 
* Install VDJtools of appropriate version (see ``build.gradle`` in repository root folder) using Maven (``mvn clean install`` from VDJtools repository root folder).
* Build VDJdb-standalone with Gradle by running ``gradle clean build`` from the repository root folder.

### Running VDJdb

Standalone VDJdb annotation utility can be executed by running the following command:

```
java -jar -Xmx4G vdjdb.jar [options] [sample1 sample2 sample3 ... if -m is not specified] output_prefix
```

First part of the command runs the JAR file and sets the memory limit to 4GB (should be increased in case JVM drops with heap size exception). The second part includes options, input samples and prefix of output files.

The list of accepted options is the following:

| Shorthand | Long  name                 | Required | Argument example                 |  Description |
|-----------|----------------------------|----------|----------------------------------|--------------|
| ``-h``    |                            |          |                                  |  Display help message |
| ``-m``    | ``--metadata``             |          | ``/path/to/metadata.txt``        |  A [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) file, holding paths to samples and user-provided information. |
|           | ``--software``             |          | ``MITCR``,``MIGEC``,etc          |  Input RepSeq data format, see [formats supported for conversion](http://vdjtools-doc.readthedocs.io/en/latest/input.html#formats-supported-for-conversion). By default expects input in [VDJtools format](http://vdjtools-doc.readthedocs.io/en/latest/input.html#vdjtools-format).
| ``-S``    | ``--species``              | Yes      | ``human``,``mouse``,etc          |  Species name. All samples should belong to the same species, only one species is allowed. |
| ``-R``    | ``--gene``                 | Yes      | ``TRA``,``TRB``,etc              |  Name of the receptor gene. All samples should contain to the same receptor gene, only one gene is allowed. |
| ``-v``    | ``--v-match``              | Yes      |                                  |  Require Variable segment matching when searching the database. |
| ``-j``    | ``--j-match``              | Yes      |                                  |  Require Joining segment matching when searching the database. |
|           | ``--search-preset``        |          | ``balanced``,``high-recall``,etc |  Sets parameters for CDR3 match search and scoring according to specified preset. Default is ``balanced``. See **CDR3 matching** section below. |
|           | ``--search-scope``         |          | ``5,2,2,7``                      |  Initial parameters used to select CDR3 sequences from the database in s/i/d/m format: allowed number of substitutions (s), insertions (i), deletions (d) and the total number of mismatches (m). Depends on preset. |
|           | ``--search-threshold``     |          | ``-2000``,``2.5e3``,etc          |  Overrides CDR3 alignment score threshold. Score is computed according to scoring scheme (pre-optimized substitution matrix and gap penalty). Not applicable for 'dummy' preset. |
|           | ``--database``             |          | ``/path/to/my_db``               |  Path and prefix of an external database. Should point to files with '.txt', and '.meta.txt' suffices (the database itself and database metadata).|
|           | ``--use-fat-db``           |          |                                  |  In case running with a built-in database, will use full database version instead of slim one. Full database contains extended info on method used to identify a given specific TCR and sample source, but has a higher degree of redundancy (several identical TCR:pMHC pairs from different publications, etc) that can complicate post-analysis |
|           | ``--filter``               |          | ``__antigen.species__ =~ "EBV"`` |  Logical filter expresstion that will be evaluated for database columns. Supports Java/Groovy syntax, Regex, .contains(), .startsWith(), etc. Parst of the expression marked with double underscore (``__``) will be subsituted with corresponding values from database rows. Those parts should be named exactly as columns in the database |
|           | ``--vdjdb-conf-threshold`` |          | ``1``                            |  VDJdb confidence level threshold, from ``0`` (lowest) to ``3`` (highest), default is ``1``. See [database readme](https://github.com/antigenomics/vdjdb-db) for details on VDJdb confidence scoring procedure |
|           | ``summary-columns``        |          | ``antigen.species,antigen.gene`` |  Table columns for which a summary output is provided for each sample, see [VDJdb specification](https://github.com/antigenomics/vdjdb-db#database-specification) and database metadata file for more information on available columns. Default is ``mhc.class,mhc.a,mhc.b,antigen.species,antigen.gene``|
| ``-c``    | ``--compress``             |          |                                  |  Compress sample-level summary output with GZIP. |

The following output files will be generated:

* ``annot.summary.txt`` annotation summary containing the number of unique clonotypes (``unique``), their cumulative share of reads (``frequency``) and total read count (``reads``). Sample metadata will be appended to this table if provided via the ``-m`` option. Each row corresponds to a combination of database field values from the columns specified by the ``--summary-columns`` option. If a single clonotype is matched to several VDJdb records, its reads count and frequency and will be appended to all of them and the ``unique`` counter for each of the records will be incremented by ``1``. Each of database records is tagged as ``entry`` in ``counter.type`` column of summary table, statistics (total number of clonotypes, read share and count) of annotated and unannotated clonotypes is stored in rows tagged as ``found`` and ``not.found`` respectively.
* ``$sample_id.annot.txt`` annotation for each of the clonotypes found in database. This is an all-to-all merge table between the sample and database that includes all matches. Clonotype information from the sample (count, frequency, cdr3 sequence, v/d/j segments and v/d/j markup) is preserved. As a clonotype can be represented by multiple rows in the output, ``id.in.sample`` column can be used to make the correspondence between annotation record and 0-based index of clonotype in the original sample. For the information on database columns that are appended see database schema in [VDJdb-db repository](https://github.com/antigenomics/vdjdb-db) readme. The ``score`` column contains CDR3 alignment score that is computed as described below (not to be confused with [VDJdb record confidence score](https://github.com/antigenomics/vdjdb-db#vdjdb-scoring).

### CDR3 matching in VDJdb

Scanning the database with CDR3 sequences from the sample is performed in two steps:

* A fast and exhaustive search using a suffix tree that is controlled by number of amino acid substitutions and indels as set by ``--search-scope`` parameter.
* Each of resulting alignments is scored using a pre-trained substitution matrix, gap penalty and positional weights as specified by ``--search-preset`` parameter. Only records passing a certain threshold (set by ``--search-threshold`` or, by default, taken from the parameter preset) are reported in the final output.

VDJdb therefore uses optimized alignment scoring to rank and filter CDR3 matches [TBA](https://en.wikipedia.org/wiki/To_be_announced). The following presets are available:

Name               | Search scope | Description
-------------------|--------------|-------------------------------------------------------------------------------------------------
``dummy``          | ``2,1,1,2``  | No scoring threshold is applied, CDR3 length minus number of mismatches is reported as ``score``
``high-recall``    | ``5,2,2,7``  | Recall near ``~100%`` on training dataset, precision is ``~5%``
``balanced``       | ``5,2,2,7``  | Precision near ``~80%`` on training dataset, recall is ``~50%``
``high-precision`` | ``5,2,2,7``  | Precision near ``~100%`` on training dataset, recall is ``~5%``

Note that ``--search-scope`` argument can be used to speed up VDJdb-standalone and lower its memory requirements by decreasing the search space. This can severely affect the recall though.
The ``--search-threshold`` is encoded into the preset and selected during the optimization procedure. Changing this threshold can provide a mean to fine-tune VDJdb-standalone when a control set is available or manual inspection of annotation results reveal certain inconsistencies. However, a better choice would be to re-run CDR3 scoring optimization with a set of user-provided TCRs with known specificity.
Values for ``--search-scope`` and ``--search-threshold`` are initially determined from the ``--search-preset`` and can be overriden by using corresponding options.

### See also

A web-based application GUI for VDJdb-standalone called [VDJdb-server](https://github.com/antigenomics/vdjdb-server) can be used to browse the database and annotate the samples.