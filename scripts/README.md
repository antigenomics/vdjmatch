# match_background.py

Estimate epitope-specific TCR precursor frequencies by matching epitope-associated clonotypes against a large **random background repertoire** generated with OLGA (via `mir.basic.pgen.OlgaModel`).

The script:

- Reads a clonotype table (VDJdb-like or user-provided).
- Groups clonotypes by **(species, chain, epitope)**.
- Generates `N` random CDR3 amino-acid sequences (and optionally V/J genes) with OLGA.
- Counts how many background sequences match each group with **exactly 1 amino-acid substitution**.
- Optionally requires V and/or J gene agreement for a match.

## Input format

By default the script expects a TSV with these columns:

- `junction_aa` — CDR3 amino-acid sequence (uppercase letters)
- `v_call` — V gene
- `j_call` — J gene
- `antigen.epitope` — epitope label for grouping
- `species` — species label
- `gene` — chain label (e.g. `TRB`, `TRA`)

You can override column names via CLI flags (see the flags table below).

## Usage examples

### 1) Minimal run (no V/J constraints)

```bash
./match_background.py   -n 1000000   --seed 42 --batch-size 10000   --gene TRB   --group-col "antigen.epitope"   input_vdjdb.tsv
```

### 2) Require V-gene agreement (`--match-v`)

```bash
./match_background.py   -n 100000000   --seed 42   --threads 8   --batch-size 20000   --match-v   --gene TRB   input_vdjdb.tsv
```

### 3) Require both V and J (`--match-v --match-j`) and fail fast if OLGA lacks V/J

```bash
./match_background.py   -n 200000000   --seed 42   --threads 12   --batch-size 20000   --match-v   --match-j   --strict-vj   --gene TRB   input_vdjdb.tsv
```

### 4) Custom column names

```bash
./match_background.py   -n 50000000   --junction-col "CDR3"   --v-col "V"   --j-col "J"   --group-col "epitope"   --species-col "species"   --chain-col "chain"   input.tsv
```

---

## Full CLI flags

Notes:
- `input_tsv` is a required positional argument.

| Flag             | Type | Default                  | Description                                                                                                   |
|------------------|---:|--------------------------|---------------------------------------------------------------------------------------------------------------|
| `input_tsv`      | positional | —                        | Input clonotype table (TSV).                                                                                  |
| `--out`          | str | `background_matches.tsv` | Output TSV path.                                                                                              |
| `--gene`         | str | `TRB`                    | Chain passed to `OlgaModel(chain=...)` (e.g. TRB/TRA).                                                        |
| `-n`             | int | (required)               | Number of background sequences to generate and search.                                                        |
| `--seed`         | int | `42`                     | Random seed used in the OLGA producer process.                                                                |
| `--batch-size`   | int | `10000`                  | Number of background sequences per search batch per worker (controls how often the progress bar updates).     |
| `--junction-col` | str | `junction_aa`            | Column name for CDR3 amino-acid sequence in input.                                                            |
| `--v-col`        | str | `v_call`                 | Column name for V gene call in input.                                                                         |
| `--j-col`        | str | `j_call`                 | Column name for J gene call in input.                                                                         |
| `--group-col`    | str | `antigen.epitope`        | Column name used as the epitope/group label.                                                                  |
| `--species-col`  | str | `species`                | Column name for species label in input.                                                                       |
| `--chain-col`    | str | `gene`                   | Column name for chain label in input (e.g. TRB/TRA).                                                          |
| `--match-v`      | flag | off                      | Require background V gene to match when counting a match.                                                     |
| `--match-j`      | flag | off                      | Require background J gene to match when counting a match.                                                     |
| `--max-sub`      | int | `1`                      | Maximum substitutions allowed. Must be `1` for this script.                                                   |
| `--max-ins`      | int | `0`                      | Maximum insertions allowed. Must be `0` for this script.                                                      |
| `--max-del`      | int | `0`                      | Maximum deletions allowed. Must be `0` for this script.                                                       |
| `--max-edits`    | int | `0`                      | Maximum total number of edit operations allowed between a background CDR3 sequence and a database clonotype.  |

