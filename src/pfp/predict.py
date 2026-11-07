import argparse
import pandas as pd
import torch

from data import read_fasta, encode
from models import DummyGOModel

BATCH_SIZE = 32
ONTOLOGIES = ["BP", "MF", "CC"]

# Command-line arguments: input FASTA and output CSV
parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

torch.manual_seed(42)   # same random scores every run

# 1. Load GO vocabularies and create one model per ontology
vocabs = {}
models = {}
for name in ONTOLOGIES:
    vocabs[name] = pd.read_csv(f"data/vocab/vocab_{name}.tsv")["go_id"].tolist()
    models[name] = DummyGOModel(n_terms=len(vocabs[name]))

# 2. Read proteins
proteins = read_fasta(args.fasta)
print("Proteins:", len(proteins))

# 3. Predict in batches and write CSV
with open(args.out, "w") as f:
    f.write("ProteinID,GOterm,Score\n")

    for start in range(0, len(proteins), BATCH_SIZE):
        batch = proteins[start:start + BATCH_SIZE]
        ids = [p[0] for p in batch]
        sequences = [p[1] for p in batch]

        x, mask = encode(sequences)                 # (B, N, D), (B, D)

        for name in ONTOLOGIES:
            scores = models[name](x, mask)          # (B, P)
            for i, protein_id in enumerate(ids):
                for k, go_term in enumerate(vocabs[name]):
                    f.write(f"{protein_id},{go_term},{scores[i, k]:.4f}\n")

print("Saved predictions to", args.out)
