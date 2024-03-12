# Methodology

We used six different methods to find interesting genes in our study: [eggNOG-mapper](https://doi.org/10.1093/molbev/msab293), [Orthofinder](https://doi.org/10.1186/s13059-019-1832-y), [TapScan](https://doi.org/10.1093/gbe/evx258), [Gene Ontology](https://geneontology.org/), best [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) hit against [Araport11](https://www.arabidopsis.org/) and [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/Introduction.html). Each of them approach a similar question with a different method and we believe the collective look can help us to target interesting genes better than individual approaches. Please read the corresponding paper for each tool if you need a comprehensive explanation of the tools.

In the context of the analyses described in this section, **it is vital to emphasize that only protein sequences from the indicated species were used.** Consequently, this collection does not include non-coding RNA sequences. Alternative methodology or customized criteria-based studies may be necessary for investigating non-coding RNA or identifying specific genes of interest, potentially necessitating independent study or the use of distinct analytical approaches.

OrthoFinder is a software for identifying and categorizing orthologous genes across species or genomes. It employs an algorithm that considers both sequence similarity and gene evolutionary relationships. By creating a similarity graph based on pairwise comparisons of protein sequences, OrthoFinder clusters genes into orthogroups, which represent genes descended from a common ancestral gene. This approach ensures accurate grouping, even for genes with complex evolutionary histories.To obtain the best results from Orthofinder, we must include a phylogenetically diverse sample of species for orthogroup inference. As a result, this section contains more species.

EggNOG-mapper is a tool for assigning functional annotations to genes by mapping them to orthologous groups and functional categories in the EggNOG database. Users submit gene sequences, and EggNOG-mapper detects orthologous groups and provides functional insights based on known functions in the database. This technology simplifies gene annotation and helps us understand gene functions and evolutionary links across species.

With the increased availability of genomic data, TapScan is a tool designed specifically for plant genome analysis to identify transcription factors with high precision.

Gene ontology (GO) is a database and standardized approach for classifying and categorizing genes and gene products based on their molecular functions, biological processes, and cellular components. It establishes a common vocabulary and framework for identifying gene features and behaviors across various animals, making it easier to evaluate Genomic data and better understand gene function in biological systems. GO organizes genes into hierarchical categories and links, allowing researchers to perform comparative analysis, identify functional similarities, and gain biological insights across species. We integrated the InterProScan and eggNOG-mapper results and filtered the "obsolute" GO IDs to generate comprehensive GO-Gene-Term sets.

InterProScan helps also with the functional annotation of protein sequences. By entering protein sequences into InterProScan, we may detect conserved domains, motifs, and functional signatures. This information contributes to our understanding of these proteins' potential roles and activities in biological processes. Furthermore, InterProScan can predict protein families and connect sequences to known biological pathways, which is critical for understanding the molecular mechanisms that underpin many biological occurrences.

**Unfortunately, the multiple sequence alignment viewer is not compatible with the Firefox browser. It's a bug in the package we use for MSA visualization, and we can't do anything about it. We checked, and MSA viewer works fine with Chrome, Brave, and Safari. Other browsers have not been checked.**

The parameters that were used for each program are listed below.

## Orthofinder

We used Orthofinder with two different settings, and used the results of the second run since we believe it is more accurate. The species tree that we used for the second run is shown below.

```{r, engine='bash'}
# 1st run
orthofinder.py -S diamond -M msa -A mafft -T fasttree -t 50 -a 6 -y -n run_1
# 2nd run
orthofinder.py -t 50 -a 6 -y -n run_2 -ft Results_run_1 -s SpeciesTree_input.txt
```

## eggNOG-mapper

```{r, engine='bash'}
emapper.py -m diamond --itype proteins --data_dir eggnog-mapper-data/ --dmnd_iterate yes --dbmem --cpu 0 --evalue 1e-10 --sensmode ultra-sensitive --tax_scope 33090 --dmnd_db eggnog-mapper-data/eggnog_proteins_default_viridiplantae.dmnd
```

## InterProScan

```{r, engine='bash'}
interproscan.sh -cpu 150 -pa -goterms  &> iprsc.log 
```