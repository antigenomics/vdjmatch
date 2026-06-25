"""Standalone TCR-BERT embedder (runs in the cmp-tcrbert conda env; torch/transformers only).

    conda run -n cmp-tcrbert python bench/_tcrbert_embed.py <cdr3_list.txt> <out.npy>
Reads one CDR3 per line, writes an (n, 768) mean-pooled last-hidden-state embedding (MLM-only model,
works for TRA and TRB).
"""
import sys

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer

torch.manual_seed(0)                                                 # reproducible embeddings
np.random.seed(0)
torch.use_deterministic_algorithms(True, warn_only=True)
inp, outp = sys.argv[1], sys.argv[2]
seqs = [ln.strip() for ln in open(inp) if ln.strip()]
tok = AutoTokenizer.from_pretrained("wukevin/tcr-bert-mlm-only")
mdl = AutoModel.from_pretrained("wukevin/tcr-bert-mlm-only").eval()


def embed(seqs, bs=256):
    out = []
    for i in range(0, len(seqs), bs):
        enc = tok([" ".join(list(s)) for s in seqs[i:i + bs]], return_tensors="pt", padding=True)
        with torch.no_grad():
            h = mdl(**enc).last_hidden_state
        m = enc["attention_mask"].unsqueeze(-1)
        out.append(((h * m).sum(1) / m.sum(1)).numpy())
    return np.vstack(out).astype(np.float32)


np.save(outp, embed(seqs))
print(f"[tcrbert] embedded {len(seqs)} -> {outp}", file=sys.stderr)
