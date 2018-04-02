[![Build Status](https://travis-ci.org/antigenomics/vdjmatch.svg?branch=master)](https://travis-ci.org/antigenomics/vdjmatch)

## VDJmatch: a software for database-guided prediction of T-cell receptor antigen specificity

VDJmatch is a command-line tool designed for matching T-cell receptor (TCR) repertoires against a database of TCR sequences with known antigen specificity. VDJmatch implements an API for interacting and querying the [VDJdb database](https://github.com/antigenomics/vdjdb-db), and serves as a backend for [VDJdb web browser](vdjdb.cdr3.net). VDJmatch will automatically download and use the latest version of the VDJdb database, however, it is also possible to use a custom database provided by user if it matches [VDJdb format specification](https://github.com/antigenomics/vdjdb-db#database-specification).

VDJmatch accepts TCR clonotype table(s) as an input and relies on [VDJtools](http://vdjtools-doc.readthedocs.org/en/latest/index.html) framework to parse the output of commonly used immune repertoire sequencing (RepSeq) processing tools. See [format](http://vdjtools-doc.readthedocs.org/en/latest/input.html) section of VDJtools docs for the list of supported formats. Note that VDJmatch can be used with [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) semantics introduced by VDJtools to facilitate running annotation for multi-sample datasets.

### Installing and running

VDJdb is distributed as an executable JAR that can be downloaded from the [releases section](https://github.com/antigenomics/vdjdb/releases), the software is cross-platform and requires [Java v1.8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) or higher to run.

To run the executable JAR use the ``java -jar path/to/vdjmatch-version.jar [options]`` command as described below. Running without any ``[options]`` or with ``-h`` option will display the help message.

The latest version of VDJdb will be downloaded the first time you run VDJmatch. Note that in order to update to the most recent version next time, you will need to run ``java -jar path/to/vdjmatch-version.jar Update`` command.

### VDJmatch command line options

The following syntax should be used to run VDJmatch for RepSeq sample(s)

```
java -Xmx4G -jar path/to/vdjmatch-version.jar \
      [options] [sample1 sample2 sample3 ... if -m is not specified] output_prefix
```

First part of the command runs the JAR file and sets the memory limit to 4GB (should be increased in case JVM drops with heap size exception) and points to VDJmatch executable JAR (``version`` should be replaced with the software version). The second part includes options, input samples and the prefix of output files.

#### General

| Option name                 | Argument example                |  Description                                              |
|-----------------------------|---------------------------------|-----------------------------------------------------------|
| ``‑h``                      |                                 |  Display help message                                     |
| ``‑m``, ``‑‑metadata``      | ``/path/to/metadata.txt``       |  A [metadata](http://vdjtools-doc.readthedocs.org/en/latest/input.html#metadata) file, holding paths to samples and user-provided information.            |
| ``‑‑software``              | ``MITCR``,``MIGEC``,etc         |  Input RepSeq data format, see [formats supported for conversion](http://vdjtools-doc.readthedocs.io/en/latest/input.html#formats-supported-for-conversion). By default expects input in [VDJtools format](http://vdjtools-doc.readthedocs.io/en/latest/input.html#vdjtools-format).                        |
| ``‑c``, ``‑‑compress``      |                                 |  Compress sample-level summary output with GZIP.          |

If ``-m`` option is specified the list of sample file names should be omitted and the list of options should be followed by ``output_prefix``.

#### Database pre-filtering

| Option name                 | Argument example                 |  Description |
|-----------------------------|----------------------------------|--------------|
| ``‑S``, ``‑‑species``       | ``human``,``mouse``,etc          |  **(Required)** Species name. All samples should belong to the same species, only one species is allowed. |
| ``‑R``, ``‑‑gene``          | ``TRA``,``TRB``,etc              |  **(Required)** Name of the receptor gene. All samples should contain to the same receptor gene, only one gene is allowed. |
|  ``‑‑filter``               | ``__antigen.species__=~"EBV"`` |  **(Advanced)** Logical filter expresstion that will be evaluated for database columns. |
|  ``‑‑vdjdb‑conf``           | ``1``                            |  VDJdb confidence level threshold, from ``0`` (lowest) to ``3`` (highest), default is ``0``.  |

The ``--filter`` option supports Java/Groovy syntax, Regex, ``.contains()``, ``.startsWith()``, etc. Parts of the expression marked with double underscore (``__``, e.g. ``__antigen.epitope__``) will be substituted with corresponding values from database rows. Those parts should be named exactly as columns in the database, see [VDJdb specification](https://github.com/antigenomics/vdjdb-db#database-specification) for the list of column names.

VDJdb confidence level used by ``--vdjdb-conf`` is assigned based on the details of TCR specificity assay for each VDJdb record, see [VDJdb confidence scoring](https://github.com/antigenomics/vdjdb-db#vdjdb-scoring) for details on this procedure.

#### Using external database (advanced)

| Option name                 | Argument example                 |  Description |
|-----------------------------|----------------------------------|--------------|
|  ``‑‑database``             | ``/path/to/my_db``               |  Path and prefix of an external database. Should point to files with '.txt', and '.meta.txt' suffices (the database itself and database metadata).|
|  ``‑‑use‑fat‑db``           |                                  |  In case running with a built-in database, will use **full database** version instead of slim one. |

**Full database** contains extended info on method used to identify a given specific TCR and sample source, but has a higher degree of redundancy (several identical TCR:pMHC pairs from different publications, etc) that can complicate post-analysis

#### Search parameters

| Option name                 | Argument example                 |  Description |
|-----------------------------|----------------------------------|--------------|
|  ``‑‑v‑match``              |                                  |  Require exact Variable segment ID match (ignoring alleles) when searching the database. |
|  ``‑‑j‑match``              |                                  |  Require exact Joining segment ID match (ignoring alleles) when searching the database. |
|  ``‑O``, ``‑‑search‑scope`` |  ``2,1,2``, ``3,0,0,3``, ...     |  Sets CDR3 sequence search parameters aka *search scope*: allowed number of substitutions (``s``), insertions (``i``), deletions (``d``) / or indels (``id``) and total number of mutations (``t``). Default is ``0,0,0`` |
|  ``‑‑search‑exhaustive``    | ``0``, ``1`` or ``2``            | Perform exhaustive CDR3 alignment: ``0`` - no (fast), ``1`` - check and select best alignment for smallest edit distance, ``2`` - select best alignment across all edit distances within search scope (slow). Default is ``1``.

Search scope should be specified in either ``s,i,d,t`` or ``s,id,t`` form. While the second form is symmetric and counts the sum of insertions and deletions (indels), the first form is not symmetric - insertions and deletions are counted with respect to the query TCR sequence (i.e. clonotype records from input samples). Total number of mutations ``t`` specifies the edit distance threshold.

> Note that VDJmatch running time can greatly increase for large (wider than ``4,2,4``) search scopes.

With a ``--search-exhaustive 2`` the algorithm will compute an exact global alignment for CDR3 sequences, which is quite slow, for a small/moderate search scope (2 or less indels) ``--search-exhaustive 1`` is effectively the same as ``--search-exhaustive 2``. Exhaustive search will choose the best alignment based on the VDJAM scoring (see below), this option has no effect if full VDJMATCH scoring is not used.

#### Scoring parameters

| Option name                 | Argument example                 |  Description |
|-----------------------------|----------------------------------|--------------|
| ``‑A``, ``‑‑scoring‑vdjmatch``     |                                  |  Use full VDJMATCH algorithm that computes full alignment score as a function of CDR3 mutations (weighted with VDJAM scoring matrix) and pre-computed V/J segment match scores.  |
|  ``‑‑scoring‑mode``         |  ``0`` or ``1``                  |  Either ``0``: scores mismatches only (faster) or ``1``: compute scoring for whole sequences (slower). Default is ``1``. |

If ``--scoring-vdjmatch`` is not set, will just count the number of mismatches and ignore V/J segment alignment.

CDR3 alignment score is computed as:

* ``--scoring-mode 0`` $(CDR3_1, CDR3_2) = \\sum_s [M(s_1,s_2) - max(M(s_1,s_1), M(s_2,s_2))] - \\sum_{g} M(g_1,g_1)$ where $s: 1 \\rightarrow 2$ stands for substitution, $g: 1 \\rightarrow '-'$ stands for gap and $M(1,2)$ is the VDJAM matrix
* ``--scoring-mode 1`` $S(CDR3_1, CDR3_2) = aln(CDR3_1, CDR3_2) - max(aln(CDR3_1, CDR3_1), aln(CDR3_2, CDR3_2))$ where $aln(CDR3_1, CDR3_2)$ is the global alignment score without gap penalty between sequences $CDR3_1$ and $CDR3_2$ using VDJAM matrix

Full score / probability of matching the same antigen is computed using a Generalized Linear Model with cloglog link as $P \\sim  S(V_1, V_2) + S(CDR1_1, CDR1_2) + S(CDR2_1, CDR2_2) + S(J_1, J_2) + S(CDR3_1, CDR3_2) - G(CDR3_1, CDR3_2)$, i.e. sum of scores for germline regions, CDR3 region and a penalty for number of gaps in CDR3 alignment $G(CDR3_1, CDR3_2)$.

#### Hit filtering and weighting

| Option name                        | Argument example                 |  Description |
|------------------------------------|----------------------------------|--------------|
| ``‑T``, ``‑‑hit‑filter‑score``     |                                  |  Drops hits with a score less than the specified threshold.  |
| ``‑X``, ``‑‑hit‑filter‑max``       |                                  |  Only select hit with maximal score for a given query clonotype (will consider all max score hits in case of ties).  |
|    ``‑‑hit‑filter‑topn``           |       ``3``                      |  Select best ``n`` hits by score (can randomly drop hits in case of ties).  |
|    ``‑‑hit‑weight‑inf``            |                                  |  Weight database hits by their 'informativeness', i.e. the log probability of them being matched by chance.  |

Note that score threshold is applied to unweighted (see below) full scores, thus has little sense to use in case ``--scoring-vdjmatch`` is not set. For VDJmatch scoring the range of scores is ``[0, 1]`` and the recommended value for the threshold lies in the range of ``0.1-0.5``. Filtered records will be removed from output annotation files and will not affect the resulting summary statistics.

The idea behind ``--hit-weight-info`` is to give less score to more redundant 'public' clonotypes that are likely to be found in many donors simply by chance. The weight is computed as $I(dbCDR3, epitope) = -log_10 P(match dbCDR3 within search scope | database / epitope)$, i.e. based on the probability of random matching within a given search scope according to VDJdb. Note that matches between CDR3 sequences specific to the same epitope are not counted, a pseudocount of 1 is added to prevent undefined result.

#### Database search summary

| Option name                 | Argument example                 |  Description |
|-----------------------------|----------------------------------|--------------|
|  ``‑‑summary‑columns``      | ``antigen.species,antigen.gene`` |  Table columns for which a summary output is provided for each sample, see [VDJdb specification](https://github.com/antigenomics/vdjdb-db#database-specification) and database metadata file for more information on available columns. |

Default summary columns are ``mhc.class,antigen.species,antigen.gene,antigen.epitope``, see [VDJdb specification](https://github.com/antigenomics/vdjdb-db#database-specification) for the full list of column names.

### VDJmatch output

The following output files will be generated:

1. ``annot.summary.txt`` annotation summary containing the number of unique clonotypes (``unique``), their cumulative share of reads (``frequency``) and total read count (``reads``).
    * Sample metadata will be appended to this table if provided via the ``-m`` option.
    * Each row corresponds to a combination of database field values from the columns specified by the ``--summary-columns`` option (e.g. epitope and parent species, ``antigen.epitope + antigen.species``). If a single clonotype is matched to several VDJdb records, its reads count and frequency and will be appended to all of them and the ``unique`` counter for each of the records will be incremented by ``1``.
    * The weight/informativeness sum of database hits for each row is stored in the ``weight`` column and can be used to scale the results, together with the ``db.unique`` column, storing the total number of unique database TCR entries for a given combination of summary columns.
    * Each of database records is tagged as ``entry`` in ``counter.type`` column of summary table, statistics (total number of clonotypes, read share and count) of annotated and unannotated clonotypes is stored in rows tagged as ``found`` and ``not.found`` respectively.
2. ``$sample_id.annot.txt`` annotation for each of the clonotypes found in database, a separate file is generated for each input sample.
    * This is an all-to-all merge table between the sample and database that includes all matches.
    * Clonotype information from the sample (count, frequency, cdr3 sequence, v/d/j segments and v/d/j markup) is preserved.
    * As a clonotype can be represented by multiple rows in the output (i.e. match to several records in the database), ``id.in.sample`` column can be used to make the correspondence between annotation record and 0-based index of clonotype in the original sample. For the information on database columns that are appended see database schema in [VDJdb-db repository](https://github.com/antigenomics/vdjdb-db) readme.
    * The ``score`` column contains CDR3 alignment score that is computed as described **Scoring parameters** section (not to be confused with [VDJdb record confidence score](https://github.com/antigenomics/vdjdb-db#vdjdb-scoring).
    * The ``weight`` column contains the weight (or informativeness) of corresponding database records, see **Hit filtering and weighting** section for details.

### Some useful notes / tricks

When running with VDJtools output, all annotations generated by VDJtools e.g. NDN size, clonotype incidence for pooled samples, frequency vector in a joint sample, ... will be preserved and VDJdb-standalone annotation columns will be added after them. Vice-versa, VDJdb-annotated samples can be used in VDJtools analysis when keeping in mind the following **1)** no conversion to VDJtools format is needed, **2)** as a single clonotype can be reported several times many descriptive statistics are not applicable. An example usage of both VDJdb-standalone and VDJtools:

```bash
$VDJTOOLS PoolSample -m metadata.txt .
$VDJDB -S human -R TRB --filter=``__antigen.species__ =~ "CMV"`` pool.aa.table.txt .
$VDJTOOLS ApplySampleAsFilter -m metadata.txt pool.aa.table.annot.txt filtered/
```

The ``filtered/`` folder will now contain all samples (and corresponding ``metadata.txt``) with CMV-specific clonotypes that can be used for further VDJtools analysis.

A web-based GUI for querying VDJdb and annotating RepSeq samples can be both accessed at [vdjdb.cdr3.net](vdjdb.cdr3.net) (public server) and installed as a local server (see [VDJdb-web](https://github.com/antigenomics/vdjdb-web) repository for details). Local installation can be configured to remove file size/search parameter limits that are enforced in public server.

### Compiling from source (advanced)

Pre-requisites: Git, Maven and Gradle.

First install the dependencies:

* Check the ``build.gradle`` file in this repository for the version of *milib* library, clone the correct version from [this repository](https://github.com/milaboratory/milib) and install with ``mvn clean install``
* Do the same for *vdjtools* library, which can be obtained from [this repository](https://github.com/mikessh/vdjtools). **Important** the *milib* version in ``pom.xml`` in *vdjtools* folder should match the one from ``build.gradle`` of *vdjmatch*.

After installing both *milib* and *vdjtools* to your local maven repository using ``mvn clean install``, navigate to the *vdjmatch* folder and build it using ``gradle clean build``.
