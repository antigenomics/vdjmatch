#!/usr/bin/env python3
"""NLV/HLA-A*02:01 TCR-pMHC structures: TCRen CDR3beta-peptide contact maps + PyMOL renders.

Structural-basis material for Appendix A5 (Deliverable 2). Three native NLV (NLVPMVATV) complexes
(3GSN, 5D2L, 5D2N) recognise the same hydrophobic peptide ridge through CDR3beta loops of distinct
length and contact pattern. For each we:

  1. fetch the PDB from RCSB (cached under tmp/nlv_struct/),
  2. run ``tcren contacts --interface tcr_peptide`` (TCRen package, antigenomics) to get the
     CDR3beta-peptide residue contact table,
  3. attach the per-contact TCRen statistical-potential energy (``tcren().value(tcr_aa, pep_aa)``),
  4. render a PyMOL panel (CDR loops + NLV peptide coloured by Kyte-Doolittle hydropathy, showing the
     CDR3beta-peptide contacts) -> figures/struct_nlv_<id>.{png,pdf},
  5. emit the combined contact/energy LaTeX table -> appendix/_tab_nlv_contacts.tex.

This module runs INSIDE the ``tcren-nb`` conda env (needs ``tcren`` + ``pymol``):
    conda run -n tcren-nb python bench/structures_nlv.py

Only public PDB-derived sequences appear in outputs.
"""
from __future__ import annotations

import subprocess
import sys
import urllib.request
from pathlib import Path

import polars as pl

REPO = Path(__file__).resolve().parent.parent
WORK = REPO / "tmp" / "nlv_struct"
MS = Path("/Users/mikesh/vcs/manuscripts/2026-vdjmatch")
FIG = MS / "figures"
FIGDATA = FIG / "data"
TEX = MS / "appendix" / "_tab_nlv_contacts.tex"
PYMOL = "/opt/homebrew/bin/pymol"

# (pdb id, TCR-beta chain, peptide chain) -- one representative complex per crystal form.
COMPLEXES = [
    ("3GSN", "B", "P"),
    ("5D2L", "F", "R"),
    ("5D2N", "E", "I"),
]
PEPTIDE = "NLVPMVATV"
CUTOFF = 5.0

# Kyte-Doolittle hydropathy (J. Mol. Biol. 1982; 157:105-132); used for PyMOL b-factor colouring.
KD = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8,
    "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2,
    "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}


def fetch_pdb(pdb_id: str) -> Path:
    """Download <pdb_id>.pdb from RCSB into WORK (cached)."""
    WORK.mkdir(parents=True, exist_ok=True)
    out = WORK / f"{pdb_id}.pdb"
    if not out.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, out)
    return out


def run_contacts(pdb: Path) -> Path:
    """Run ``tcren contacts`` (TCR<->peptide interface) and return the CSV path."""
    out = WORK / f"{pdb.stem}_contacts.csv"
    subprocess.run(
        ["tcren", "contacts", "-s", str(pdb), "-o", str(out),
         "--interface", "tcr_peptide", "--cutoff", str(CUTOFF)],
        check=True, capture_output=True, text=True,
    )
    return out


def cdr3b_peptide_energies(csv: Path, beta_chain: str, pep_chain: str) -> pl.DataFrame:
    """CDR3beta<->peptide contacts for one complex, with TCRen per-pair energy.

    Collapses atom-level contacts to one row per (TCR residue, peptide residue) at the minimum
    distance, attaches the TCRen statistical-potential value (TCR aa -> peptide aa)."""
    from tcren.potential import tcren
    pot = tcren()

    df = pl.read_csv(csv).filter(
        (pl.col("region.type.from") == "CDR3")
        & (pl.col("chain.type.from") == "TRB")
        & (pl.col("chain.id.from") == beta_chain)
        & (pl.col("chain.id.to") == pep_chain)
    )
    if df.height == 0:
        return df
    # residue-level: minimum atom distance per residue pair
    res = (df.group_by("residue.index.from", "residue.aa.from", "pos.from",
                       "residue.index.to", "residue.aa.to", "pos.to")
             .agg(pl.col("dist").min().alias("dist"))
             .sort("residue.index.from", "residue.index.to"))

    def energy(a, b):
        try:
            return round(pot.value(a, b), 3)
        except KeyError:
            return None

    return res.with_columns(
        pl.struct(["residue.aa.from", "residue.aa.to"])
        .map_elements(lambda s: energy(s["residue.aa.from"], s["residue.aa.to"]),
                      return_dtype=pl.Float64)
        .alias("tcren")
    )


def cdr3b_sequence(pdb: Path, beta_chain: str) -> str:
    """Full CDR3beta sequence of the chosen beta chain (via ``tcren annotate``)."""
    out = WORK / f"{pdb.stem}_annot.csv"
    if not out.exists():
        subprocess.run(["tcren", "annotate", "-s", str(pdb), "-o", str(out), "--regions", "tcr"],
                       check=True, capture_output=True, text=True)
    d = pl.read_csv(out).filter(
        (pl.col("region.type") == "CDR3") & (pl.col("chain.type") == "TRB")
        & (pl.col("chain.id") == beta_chain)
    ).sort("residue.index")
    return "".join(d["residue.aa"].to_list())


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def write_pymol_script(pdb: Path, beta_chain: str, pep_chain: str,
                       contacts: pl.DataFrame, png: Path) -> Path:
    """Generate a batch .pml that colours TCR CDR loops + NLV peptide by KD hydropathy and draws
    the CDR3beta-peptide contacts; ray-traced PNG at >=300 dpi."""
    # set b-factor = KD hydropathy for every residue so spectrum colouring is by hydropathy
    kd_lines = "\n".join(
        f'alter (resn {r}), b={v}' for r, v in KD.items()
    )
    # contact distance objects between specific CDR3b residue and peptide residue
    dist_cmds = []
    for i, row in enumerate(contacts.iter_rows(named=True)):
        ti = row["residue.index.from"]
        pi = row["residue.index.to"]
        dist_cmds.append(
            f'distance ct{i}, (chain {beta_chain} and resi {ti} and not name N+C+O), '
            f'(chain {pep_chain} and resi {pi})'
        )
    dist_block = "\n".join(dist_cmds)
    # alpha chain id: pick the partner TRA chain near this beta (annot has chain types)
    lo = contacts["residue.index.from"].min()
    hi = contacts["residue.index.from"].max()
    script = f"""# auto-generated: {pdb.stem} CDR3b-NLV hydropathy panel
load {pdb}, cplx
hide everything, cplx
bg_color white
set ray_opaque_background, 0
set orthoscopic, 1
set cartoon_side_chain_helper, 1
set surface_quality, 1

# one complex -- TCR beta {beta_chain} plus peptide {pep_chain}, the CDR3b-peptide interface
create scene_obj, cplx and (chain {beta_chain} or chain {pep_chain})
delete cplx

# Kyte-Doolittle hydropathy -> b-factor, then spectrum-colour by it (blue hydrophilic .. red hydrophobic)
{kd_lines}
rebuild

# NLV peptide: a hydropathy-coloured surface + sticks (the ridge the TCR reads)
show surface, scene_obj and chain {pep_chain}
set transparency, 0.25, scene_obj and chain {pep_chain}
show sticks, scene_obj and chain {pep_chain}
set stick_radius, 0.22, scene_obj and chain {pep_chain}

# CDR3 beta: thin cartoon of the loop + sticks for the contacting apex residues
show cartoon, scene_obj and chain {beta_chain} and resi {lo - 2}-{hi + 2}
set cartoon_transparency, 0.45, scene_obj and chain {beta_chain}
show sticks, scene_obj and chain {beta_chain} and resi {lo}-{hi}
set stick_radius, 0.18, scene_obj and chain {beta_chain}

spectrum b, blue_white_red, scene_obj and (chain {pep_chain} or chain {beta_chain}), minimum=-4.5, maximum=4.5
color grey70, scene_obj and chain {beta_chain} and not (resi {lo}-{hi})

# CDR3b-peptide contacts as thin dashes
{dist_block}
hide labels
set dash_color, grey50
set dash_gap, 0.35
set dash_radius, 0.04
set dash_width, 1.5

set ray_shadows, 0
set antialias, 2
orient scene_obj and chain {pep_chain}
turn x, -15
zoom scene_obj and (chain {pep_chain} or (chain {beta_chain} and resi {lo}-{hi})), 4
ray 1800, 1400
png {png}, dpi=300
"""
    pml = WORK / f"{pdb.stem}.pml"
    pml.write_text(script)
    return pml


def render(pml: Path) -> None:
    subprocess.run([PYMOL, "-cq", str(pml)], check=True, capture_output=True, text=True)


def write_tex(all_rows: list[tuple[str, str, pl.DataFrame]]) -> None:
    """Combined CDR3beta-peptide contact/energy table (booktabs, plain \\caption)."""
    lines = [
        r"\begin{table}[h]\centering\small",
        (r"\caption{TCRen CDR3$\beta$--peptide contacts for three native NLV (\texttt{NLVPMVATV}) "
         r"/ HLA-A*02:01 complexes (contact cutoff 5\,\AA). Each row is a CDR3$\beta$ residue in "
         r"contact with an NLV peptide residue, with the TCRen statistical-potential energy "
         r"(more negative = more favourable). Three CDR3$\beta$ loops of distinct length and contact "
         r"footprint engage the same hydrophobic peptide ridge.}"),
        r"\begin{tabular}{@{}llllr@{}}\toprule",
        r"Complex & CDR3$\beta$ & CDR3$\beta$ res. & NLV res. & TCRen \\\midrule",
    ]
    for pid, seq, df in all_rows:
        first = True
        for row in df.iter_rows(named=True):
            tcr = f"{THREE_TO_ONE.get(row['residue.aa.from'], row['residue.aa.from'])}{row['residue.index.from']}"
            pep = f"{THREE_TO_ONE.get(row['residue.aa.to'], row['residue.aa.to'])}{row['residue.index.to']}"
            e = row["tcren"]
            estr = f"{e:+.3f}" if e is not None else "--"
            label = f"{pid} & \\texttt{{{seq}}}" if first else " & "
            lines.append(f"{label} & {tcr} & {pep} & {estr} \\\\")
            first = False
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    lines += [r"\end{tabular}\end{table}", ""]
    TEX.write_text("\n".join(lines))


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    FIGDATA.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for pid, bchain, pchain in COMPLEXES:
        pdb = fetch_pdb(pid)
        csv = run_contacts(pdb)
        contacts = cdr3b_peptide_energies(csv, bchain, pchain)
        seq = cdr3b_sequence(pdb, bchain)
        # persist the per-complex contact/energy table to figures/data
        contacts.write_csv(FIGDATA / f"struct_nlv_{pid.lower()}_contacts.csv")
        n = contacts.height
        emean = contacts["tcren"].drop_nulls().mean()
        print(f"{pid}: CDR3b={seq} len={len(seq)} chain {bchain}/{pchain} "
              f"contacts={n} mean_tcren={emean:+.3f}")
        png = FIG / f"struct_nlv_{pid.lower()}.png"
        pml = write_pymol_script(pdb, bchain, pchain, contacts, png)
        render(pml)
        # PNG -> PDF (ray-traced, 300 dpi) via pymol's png is raster; convert with sips/pillow
        pdf = FIG / f"struct_nlv_{pid.lower()}.pdf"
        png_to_pdf(png, pdf)
        print(f"  wrote {png.name}, {pdf.name}")
        all_rows.append((pid, seq, contacts))
    write_tex(all_rows)
    print(f"wrote {TEX}")


def png_to_pdf(png: Path, pdf: Path) -> None:
    """Wrap the ray-traced PNG into a PDF at 300 dpi (Pillow)."""
    from PIL import Image
    im = Image.open(png).convert("RGB")
    im.save(pdf, "PDF", resolution=300.0)


if __name__ == "__main__":
    sys.exit(main())
