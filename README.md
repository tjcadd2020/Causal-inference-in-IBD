# Causal-inference-in-IBD
Gao, Sheng, et al. "IBD Subtype-Regulators IFNG and GBP5 Identified by Causal Inference Drive More Intense Innate Immunity and Inflammatory Responses in CD Than Those in UC." Frontiers in Pharmacology 13 (2022).

# Usage
python dowhy_IBD.py

# Input
1. prior graph edges, which can be downloaded such as .kgml or .gml files in KEGG website, then convert the format into source/target in each row.
2. dataframe with gene expressions: row as sample, col as geneID (the last col as condition)

