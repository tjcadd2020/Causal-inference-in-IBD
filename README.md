# Causal-inference-in-IBD
Gao, Sheng, et al. "IBD Subtype-Regulators IFNG and GBP5 Identified by Causal Inference Drive More Intense Innate Immunity and Inflammatory Responses in CD Than Those in UC." Frontiers in Pharmacology 13 (2022).

# Usage
```python
python dowhy_IBD.py
```

# Input
1. A prior graph edges file: it can be downloaded such as .kgml or .gml files in KEGG website, then convert it into a dataframe (see in example files).
2. A dataframe with gene expressions: row as sample, col as geneID (the last col as condition)

# Output
CI_result file with causal inference estimation value of each gene.

