import torch
from torch import nn


class DummyGOModel(nn.Module):
    """Dummy model for one GO ontology.

    Input:  x    (B, N, D) - one-hot encoded sequences
            mask (B, D)    - True for valid residues, False for padding
    Output: scores (B, P)  - random scores in [0, 1], P = number of GO terms
    """

    def __init__(self, n_terms):
        super().__init__()
        self.n_terms = n_terms

    def forward(self, x, mask):
        B = x.shape[0]
        return torch.rand(B, self.n_terms)


# Quick test when this file is run directly
if __name__ == "__main__":
    model = DummyGOModel(n_terms=2505)
    x = torch.zeros(32, 22, 1024)
    mask = torch.ones(32, 1024, dtype=torch.bool)
    scores = model(x, mask)
    print("Output shape:", scores.shape)
    print("Min / max score:", scores.min().item(), scores.max().item())
