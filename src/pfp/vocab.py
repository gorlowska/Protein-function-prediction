import pandas as pd

# 1. Load data
vocab = pd.read_csv("data/processed/train_label_vocabulary.tsv", sep="\t")
evidence = pd.read_csv("data/processed/label_evidence.tsv", sep="\t")

# 2. Map each GO term to its ontology letter (P, F or C)
go_to_aspect = dict(zip(evidence["go_id"], evidence["aspect"]))

# 3. Add the ontology letter to every training GO term
vocab["aspect"] = vocab["go_id"].map(go_to_aspect)

# 4. Save one vocabulary file per ontology
ontologies = {"BP": "P", "MF": "F", "CC": "C"}

for name, letter in ontologies.items():
    terms = vocab[vocab["aspect"] == letter]["go_id"]
    terms.to_csv(f"data/vocab/vocab_{name}.tsv", index=False)
    print(name, len(terms))
