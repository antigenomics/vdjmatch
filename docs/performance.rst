Laptop runtime, RAM and thread budget
=====================================

The worked examples use four native threads as a practical explicit setting.
The runtime default remains one: it is conservative for restricted CPU
allocations and tiny inputs. Choose a budget no larger than the process CPU
allocation. Samples run sequentially; native search batches the entire query
population, without an outer pool. Numerical helper pools are limited to one
thread before imports in the notebooks. No runtime/scorer default changed here.

Measured commodity laptop
-------------------------

The following whole-CLI observations were recorded on 14 October 2026 on an
Apple M3 Max laptop, 16 physical/logical cores and 48 GiB RAM (51.54 decimal GB).
All 16 cores were available to the host; the table varies the requested native
worker budget rather than asserting OS CPU affinity. The interpreter was
Python 3.13.7, with vdjmatch 0.4.0, seqtree 1.1.1, vdjtools 5.1.0 and Polars 1.41.2.
The pinned reference is ``2026-06-11-ZENODO/vdjdb_full.txt.gz`` from the benchmark
mirror; its SHA256 is
``2d585775222f7e16c722671aed3a8b59a37473c36e671a9f13d9ed12771f6f80``.

.. list-table:: Measured whole-invocation cost; decimal GB, one fresh invocation per row
   :header-rows: 1
   :widths: 32 19 8 10 11 10 10

   * - Workload
     - Observations
     - Native workers
     - Time (s)
     - Peak tree RSS (GB)
     - Mean occupied cores
     - Peak live threads
   * - Assayed TRB, NLV/ELA/LLL/LLW
     - 2,670 query rows
     - 1
     - 13.71
     - **0.841**
     - 0.99
     - 6
   * - Same complete TRB workload
     - 2,670 query rows
     - 4
     - **4.69**
     - 0.846
     - 3.01
     - 7
   * - Baseline repertoire, LLW target
     - 10,000 selected clonotypes
     - 4
     - 21.76
     - 3.530
     - 3.49
     - 9
   * - Three complete precursor sets
     - 22,164 usable unique junctions
     - 4
     - 70.81
     - 1.049
     - 3.39
     - 10

Bold identifies the better measured value only within the two comparable TRB
rows. These are measured examples, not fitted extrapolations or benchmark medians.
The precursor input contains 34,162 source observations and 22,171 distinct
raw junctions across three peptides; seven fail the model's junction check.
The vaccination example selects the 10,000 most abundant V/junction/J clonotypes
from 578,812 eligible unique clonotypes in one baseline table, representing
124,390 / 833,386 eligible read counts (14.93%). This abundance-enriched subset
does not establish cost or specificity prevalence for a complete cohort.

The one/four-thread TRB evidence tables are exactly equal, including source
identities, pMHC, all six component scores and unsupported observations.
The measured wall-time ratio is 2.92 (13.71 / 4.69 seconds). A second clean-kernel
verification gave 13.60 and 4.79 seconds, with equal evidence again and peak RSS
0.937 and 0.948 GB. Reference loading, startup and table export contribute to
elapsed time, so four workers do not imply fourfold speedup.

Are extra live threads oversubscription?
----------------------------------------

CPU seconds divided by elapsed seconds give occupied-core equivalents: 1.0
means an average 100% of one core. The four-worker examples averaged 3.01–3.49
occupied cores. Peak 50-ms sampled core equivalents were 4.19 for TRB, 4.46 for
baseline annotation and 4.47 for the precursor run. Short samples include timing
and CPU-accounting granularity; these measurements cannot establish a strict
instantaneous four-core ceiling. They show useful multi-core work without a
sustained whole-run excess over the requested native worker budget.

Live-thread counts also include main, runtime-helper and idle native/table-pool
threads. Seven to ten live threads therefore do not alone imply seven to ten
simultaneous CPU-bound workers. Record both CPU activity and live threads when
checking a new installation; do not infer either from a version string.

Measurement scope and limitations
---------------------------------

``examples/profile_cli.py`` samples the CLI and all live descendants every
50 ms, summing concurrent RSS and CPU time and recording elapsed wall time.
Wall time includes CLI startup, loading, calculation and export. Notebook
preparation, retained plotting arrays, the notebook kernel and figure rendering
are outside the CLI tree: total interactive-session RAM will be higher. Shared
pages may be counted in more than one process, and peaks between samples may
be missed. CPU is an approximate sampled process-tree total. No per-query
workers, computed caches or reused predictions are involved.

These observations establish the cost of the listed workloads on this host.
They do not establish the complete manuscript rebuild gate or guarantee the
same cost on an x86 laptop, a four-core host or a smaller-memory machine.
Re-measure on the intended hardware with the same raw inputs and CLI recipe.
Every biological notebook writes a fresh resource TSV beside its final tables;
complete local run logs and source-identifier mappings remain private.

* :download:`Exact aggregate resource observations <../examples/resources.tsv>`
* :download:`Process-tree measurement helper <../examples/profile_cli.py>`
* :doc:`Worked notebook gallery <notebooks>`
