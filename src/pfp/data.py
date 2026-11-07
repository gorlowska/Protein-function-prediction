import skbio
import torch

# 20 standard amino acids + X (unknown) + PAD (padding) = 22
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
X_INDEX = 20
PAD_INDEX = 21
N = 22        # amino-acid vocabulary size
D = 1024      # maximum sequence length


def read_fasta(path):
    """Return a list of (protein ID, sequence) pairs."""
    proteins = []
    for record in skbio.io.read(path, format="fasta"):
        protein_id = record.metadata["id"]
        sequence = str(record)
        proteins.append((protein_id, sequence))
    return proteins


def encode(sequences):
    """Encode a list of sequences as a one-hot tensor (B, N, D) and a mask (B, D)."""
    B = len(sequences)
    x = torch.zeros(B, N, D)
    mask = torch.zeros(B, D, dtype=torch.bool)

    for i, seq in enumerate(sequences):
        seq = seq[:D]                          # truncate to D residues
        for j in range(D):
            if j < len(seq):
                aa = seq[j]
                if aa in AMINO_ACIDS:
                    k = AMINO_ACIDS.index(aa)
                else:
                    k = X_INDEX                # non-standard residue
                x[i, k, j] = 1
                mask[i, j] = True              # valid residue
            else:
                x[i, PAD_INDEX, j] = 1         # padding
    return x, mask


# Quick test when this file is run directly
if __name__ == "__main__":
    proteins = read_fasta("data/processed/splits/test.fasta")
    print("Number of proteins:", len(proteins))
    print("First protein:", proteins[0][0], proteins[0][1][:30], "...")

    x, mask = encode(["ACDXU", "M" * 2000])
    print("x shape:", x.shape)
    print("mask shape:", mask.shape)
    print("Valid residues per sequence:", mask.sum(dim=1))
