# Introduction

Gene regulatory networks are the intricate webs of interactions between genes and regulatory elements that control biological activities. Understanding these networks is essential for determining the underlying mechanisms that control biological systems. Time-series RNA-Seq data provide a unique chance to dynamically infer gene regulatory networks by capturing the temporal dynamics of gene expression changes in response to different stimuli or perturbations. By studying gene expression profiles across time, we can discover regulatory relationships, identify important regulatory genes, and figure out the dynamic interactions that shape gene expression programs. Time-series RNA-Seq data are thus an effective tool for deciphering the temporal characteristics of gene regulation networks and understanding how they influence cellular behaviours and responses.

# Sliding Window Inference for Network Generation (SWING)

Despite the availability of time-resolved, high-throughput data, many algorithms ignore the temporal delays inherent in regulatory systems, resulting in unreliable network inferences. SWING is used to address this issue by only taking temporal information into account when identifying time-delayed edges. SWING's tolerance to user-defined parameters allows for the successful identification of regulatory mechanisms from time-series gene expression data. SWING uses multivariate Granger causality to capture the regulatory relationships between genes throughout time. SWING, unlike traditional Granger approaches, uses a sliding window approach to evaluate numerous upstream regulators at the same time over a range of time delays. If you would like to investigate more about this tool, please see [this paper](https://doi.org/10.1073/pnas.1710936115).

The algorithm is O(N<sup>2</sup>), indicating that the number of input genes has a significant impact on the computing time. It was not possible to use all expressed genes as input for this tool. We filtered for transcription factors as well as stress response gene homologs of *A. thaliana* (identified using sequence similarity search). For approximately 3,000 genes and 5 to 9 time points, the calculation took 4 to 16 days to complete. As a result, these databases do not include all genes or their interactions. We did two analyses. First, we include all time points where we measured RNA-Seq and metabolite levels for each species under treatment. Second, we only considered RNA-Seq data for each species subjected to a treatment. For *M. endlicherianum* in both instances and transcript-only data sets of *P. patens* and *Z. circumcarinatum-1b*, the following settings were used as input: `k_min = 0, k_max = 1, w = 4, method = RandomForest, trees = 500, and lag_method='mean_mean'`. For *P. patens* and *Z. circumcarinatum-1b*, we utilized the following settings for the data set of RNA-Seq and metabolites: `k_min = 0, k_max = 1, w = 2, method = RandomForest, trees = 500, lag_method="mean_mean"`.

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