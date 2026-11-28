Annotate single and paired receptors
====================================

Run one command and read the vX ``score``. Larger values indicate stronger
reference evidence. Use the matching calibration described in
:doc:`explanation` to decide whether to call a match.

The annotation workflow uses the curated manuscript reference population by default.
Pass ``--reference-population all-valid`` for the preserved all-valid baseline.
Exact evidence remains included in ordinary annotation.

Single-chain input
------------------

A tab-separated file contains one receptor per row. ``junction_aa`` includes
the conserved terminal residues; ``v_call`` and ``locus`` identify its V gene
and chain. These are the supplied synthetic example rows:

.. literalinclude:: ../examples/inputs/cohort/sample_a.tsv
   :language: text

Run
~~~

.. code-block:: sh

   vdjmatch annotate examples/inputs/cohort/sample_a.tsv \
     --vdjdb examples/inputs/reference.tsv --locus TRB \
     --threads 4 -o single

Read the output
~~~~~~~~~~~~~~~

``single.selection.tsv`` contains one row per input receptor. Selected columns
from the verified example are shown below; the full table preserves input IDs.

.. list-table:: Single-chain result
   :header-rows: 1

   * - sequence_id
     - selected_epitope
     - score
     - selection_status
   * - q1
     - DEMO_A
     - 0.01875
     - ranked_evidence
   * - q2
     - DEMO_A
     - 0.01875
     - ranked_evidence

Both receptors rank the synthetic label DEMO_A first. These scores rank
reference evidence. Add a matched null bank to apply a call threshold.

Paired input
------------

Put observed alpha and beta chains on the same row, with a stable ``clone_id``.
Do not manufacture pairs from independent single-chain receptors.

.. literalinclude:: ../examples/inputs/paired.tsv
   :language: text

Run
~~~

.. code-block:: sh

   vdjmatch annotate examples/inputs/paired.tsv \
     --vdjdb examples/inputs/reference.tsv --locus paired \
     --threads 4 -o paired

Read the output
~~~~~~~~~~~~~~~

``paired.selection.tsv`` reports evidence against genuine paired reference keys.

.. list-table:: Paired result
   :header-rows: 1

   * - clone_id
     - selected_epitope
     - score
     - selection_status
   * - pair1
     - DEMO_A
     - −0.45
     - ranked_evidence
   * - pair2
     - unavailable
     - unavailable
     - no_support

The complete pair ranks DEMO_A first. The second row lacks a beta chain and
keeps its identity and missing-chain disposition. Its alpha evidence is in
``paired.TRA.selection.tsv``. The marginal tables include unlinked reference
chains and report ranked evidence independently of paired calls.

Parameters that change the question
-----------------------------------

.. list-table:: Main options
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``--vdjdb reference.tsv``
     - Fix the reference used for all comparisons. Omit to fetch the current release.
   * - ``--locus TRA|TRB|paired``
     - Choose the receptor population; automatic routing is the default.
   * - ``--epitope LLWNGPMAV``
     - Restrict candidate evidence to the requested peptide.
   * - ``--mhc-a 'HLA-A*02' --mhc-match compatible``
     - Restrict to compatible HLA-A*02 annotations at the supplied resolution.
   * - ``--exclude-exact``
     - Remove exact query evidence for a benchmark. Ordinary annotation keeps it.
   * - ``--details``
     - Use the kernel workflow to export all candidate evidence as well as the selected label.
   * - ``--threads 4``
     - Use a budget of four CPU workers for the same score calculation.

Install and reproduce everything above
--------------------------------------

.. code-block:: sh

   python -m venv .venv
   . .venv/bin/activate
   git clone --depth 1 https://github.com/antigenomics/vdjmatch.git
   cd vdjmatch
   python -m pip install .
   sh examples/annotate.sh annotation-demo

The script runs the supplied single-chain cohort and paired inputs from scratch.
An invalid V gene in the second cohort file produces ``no_support`` rather than
losing the row. For real data, start with these same commands and replace the
input and reference paths. See :doc:`explanation` for significance and
:doc:`algorithms` for the scoring steps.

Testing a repertoire for enrichment
-----------------------------------

For the compiled bulk reference-panel calculation use ``repertoire-exact``.
Supply a sample manifest, VDJdb reference, distinct junction/V/J probes and
model-generated V/J counts. Optional ``joint_pgen`` probe values preserve a
specific probability model asset. Choose ``--sampler r --pool-mhc`` only when
reproducing the original single-epitope union; ordinary target labels retain
full pMHC. Declare correction groups with ``--fdr-group`` and inspect both
source-row and read-abundance endpoints. The complete command is in :doc:`cli`.

The flexible mode below uses marginal Pgen and distinct-receptor counting,
so it is a separate comparison rather than the same calculation.


Use ``vdjmatch repertoire-enrichment samples.tsv --null reference`` for a fixed
observed repertoire tested against Pgen-matched mock epitope corpora. The sample
manifest contains ``sample_id`` and ``path``; pass the VDJdb reference and raw
control population explicitly. ``--mock-policy pgen`` samples reference-size
controls, ``pgen-vj`` also preserves V/J usage, and ``pgen-different-v`` changes V
and junction identities. Set the policy, matching radius, replicate count and
significance family before examining outcome scores. See :doc:`cli` for a complete
command and paired-input requirements. A significant sample is enriched relative
to that null; its contributing receptors are not thereby experimentally validated.

A separate filtering-first comparison is available with ``repertoire-tcrnet``.
It runs the existing TCRnet neighbourhood test on samples and reference targets,
then flexible VDJmatch lookup. Use its raw-denominator counts, longitudinal
changes and linked assay evidence; it does not emit repertoire P-values. Keep
its filter FDR and matching scope explicit and separate from reference-null
enrichment settings.
