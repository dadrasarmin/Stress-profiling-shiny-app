# Introduction

Co-expression network analysis with RNA-Seq data provides insights into the complex regulatory mechanisms that control gene expression. We can find functional relationships, identify important regulatory genes, and discover biological pathways by analyzing the correlation patterns between gene expression profiles across diverse situations. This method enables the discovery of co-regulated genes that may be implicated in shared biological processes or pathways, providing a comprehensive understanding of gene expression dynamics. Furthermore, co-expression networks provide a systematic framework for selecting candidate genes for additional experimental validation, making it possible to uncover novel biomarkers in a variety of biological situations.

## Weighted Gene Co-Expression Network Analysis (WGCNA)

WGCNA is a widely used method in bioinformatics for creating and analyzing co-expression networks. WGCNA identifies modules of highly connected genes, allowing to discover biologically significant gene clusters associated with specific phenotypes or conditions. By weighting gene-gene correlations based on their significance, WGCNA provides a strong framework for finding co-expression patterns and choosing genes with high relevance to the biological context under study. For more information about this tool, see [this paper](https://doi.org/10.1186/1471-2105-9-559). For network construction, we used the following settings: `merge threshold = 0.20, network type = signed, TOM type = signed, Min. module size = 30, Max. P. outliers = 0.05`. We used `20` as soft threshold for *M. endlicherianum* and *P. patens* and `13` for *Z. circumcarinatum* based on our screening for scale-free network properties. 

Please use the following acronym table for measuring metabolites:

| Acronym | Full name |
|:---------------|:------|
| 11_cis | 11-cis-beta-Carotene |
| 6MHO | 6MHO |
| 9_cis | 9-cis-beta-Carotene |
| 9_cis_neo | 9-cis-Neoxanthin |
| A_Z_per_V_A_Z | (Antheraxanthin+Zeaxanthin) / (Violaxanthin+Antheraxanthin+Zeaxanthin) |
| alpha_car | alpha-Carotene |
| anthera | Antheraxanthin |
| beta_car | beta-Carotene |
| beta_car_per_9_cis_neox | beta-Carotene / 9-cis-beta-Carotene |
| beta_Car_per_beta_CC | beta-Carotene / beta-Cyclocitral |
| beta_Car_per_beta_Io | beta-Carotene / beta-Ionone |
| beta_Car_per_DHA | beta-Carotene / DHA |
| beta_Car_per_beta_CC_beta_Io_DHA | beta-Carotene / (beta-Cyclocitral+beta-Ionone+DHA) |
| beta_CC | beta-Cyclocitral |
| beta_Io | beta-Ionone |
| Chla | Chlorophyll a |
| Chla_per_b | Chlorophyll a / Chlorophyll b |
| Chlb | Chlorophyll b |
| DHA | DHA |
| DHA_per_beta_Io | DHA/beta_Ionone |
| Lut | Lutein |
| viol | Violaxanthin |
| V_A_Z_per_Chla_b | (Violaxanthin+Antheraxanthin+Zeaxanthin) / (Chlorophyll b) |
| Zeax | Zeaxanthin |

**Please be patient. Loading the results in this section take a few seconds.**
