# Bundled resources

Provenance of the scoring resources carried over from the legacy Java vdjmatch
(`legacy-java` branch, tag `1.3.1`). These are the v1 ("vdjam-2018") defaults; the
manuscript re-derives segment/chain-specific matrices ("vdjam-2026") that drop in here
without code changes.

| File | Origin | Format | Use |
|------|--------|--------|-----|
| `vdjam.txt` | legacy `src/main/resources/vdjam.txt` | TSV `aa.1 aa.2 score` (similarity) | TCR CDR3 substitution similarity → `seqtree.SubstitutionMatrix.from_similarity` |
| `score_coef.txt` | legacy `src/main/resources/score_coef.txt` | TSV per (species, gene, fun) | cloglog GLM coefficients over S(V),S(CDR1),S(CDR2),S(J),S(CDR3),indels |
| `segm_score.txt` | legacy `src/main/resources/segm_score.txt` | TSV V/J + CDR1/CDR2 pairwise scores | segment similarity terms S(V), S(J) |

`vdjam.txt` similarity → penalty follows seqtree's convention `pen(a,b)=s_aa+s_bb−2·s_ab`,
which equals the legacy "subtract self-score" rule.
