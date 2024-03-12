################################################################################
############### Set up
################################################################################
library(tidyverse)
library(ggforce)
library(ggdist)
library(gghalves)
library(treeio)
library(ggtree)
library(DT)
library(shiny)
library(shinydashboard)
library(htmltools)
library(msaR)
library(ape)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(visNetwork)
library(igraph)
################################################################################
############### Perform one time calculations here (optimization)
################################################################################
# find the genes of interests
interproscan_mesotaenium <- read_tsv("assets/interproscan/M_endlicherianum.fa.tsv", col_names = F)
colnames(interproscan_mesotaenium) <- c("Gene", "Protein accession", "Sequence MD5 digest", "Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description", "GO annotations", "Pathways annotations")
interproscan_mesotaenium <- interproscan_mesotaenium[, -c(3, 4, 11, 12)]
interproscan_physcomitrium <- read_tsv("assets/interproscan/P_patens.fa.tsv", col_names = F)
colnames(interproscan_physcomitrium) <- c("Gene", "Protein accession", "Sequence MD5 digest", "Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description", "GO annotations", "Pathways annotations")
interproscan_physcomitrium <- interproscan_physcomitrium[, -c(3, 4, 11, 12)]
interproscan_zygnema <- read_tsv("assets/interproscan/Z_circumcarinatum_SAG698_1b.fa.tsv", col_names = F)
colnames(interproscan_zygnema) <- c("Gene", "Protein accession", "Sequence MD5 digest", "Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description", "GO annotations", "Pathways annotations")
interproscan_zygnema <- interproscan_zygnema[, -c(3, 4, 11, 12)]
eggnog_mesotaenium <- read_tsv("assets/eggNOG/mesotaenium_endlicherianum.emapper.annotations.tsv", skip = 4)
eggnog_physcomitrium <- read_tsv("assets/eggNOG/physcomitrium_patens.emapper.annotations.tsv", skip = 4)
eggnog_zygnema <- read_tsv("assets/eggNOG/zygnema_1b.emapper.annotations.tsv", skip = 4)
eggnog_mesotaenium <- eggnog_mesotaenium[, c("#query", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category")]
eggnog_physcomitrium <- eggnog_physcomitrium[, c("#query", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category")]
eggnog_zygnema <- eggnog_zygnema[, c("#query", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category")]
orthofinder_N0 <- read_tsv("assets/orthofinder/N0.tsv")
tapscan_mesotaenium <- read_csv("assets/tapscan/MESEN.fa.output.1.csv")
tapscan_zygnema <- read_csv("assets/tapscan/ZYGCI.fa.output.1.csv")
tapscan_physcomitrium <- read_csv("assets/tapscan/Ppatens_318_v3.3.protein_primaryTranscriptOnly_shortID.fa-2.csv")
GO_mesotaenium <- read_tsv("assets/GO_gene_term/mesotaenium_GO.tsv")
GO_zygnema <- read_tsv("assets/GO_gene_term/zygnema_1b_GO.tsv")
GO_physcomitrium <- read_tsv("assets/GO_gene_term/physcomitrium_GO.tsv")
blast_mesotaenium <- read_tsv("assets/best_ath_blast_hit/mesotaenium.tsv")
blast_zygnema <- read_tsv("assets/best_ath_blast_hit/zygnema.tsv")
blast_physcomitrium <- read_tsv("assets/best_ath_blast_hit/physcomitrium.tsv")
# Gene expression viewer
# Define a scale function for z-score transformation
scale_z_score <- function(x){ (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }
# There is an error and it can be solved via:
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
mesotaenium_expr <- read_tsv("assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth_processed.tsv")
mesotaenium_expr$condition <- factor(mesotaenium_expr$condition, levels = unique(mesotaenium_expr$condition))
mesotaenium_expr$treatment <- factor(mesotaenium_expr$treatment, levels = unique(mesotaenium_expr$treatment))
physcomitrium_expr <- read_tsv("assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth_processed.tsv")
physcomitrium_expr$condition <- factor(physcomitrium_expr$condition, levels = unique(physcomitrium_expr$condition))
physcomitrium_expr$treatment <- factor(physcomitrium_expr$treatment, levels = unique(physcomitrium_expr$treatment))
zygnema_expr <- read_tsv("assets/gene_profiler/zygnema_count_voom_transformed_qsmooth_processed.tsv")
zygnema_expr$condition <- factor(zygnema_expr$condition, levels = unique(zygnema_expr$condition))
zygnema_expr$treatment <- factor(zygnema_expr$treatment, levels = unique(zygnema_expr$treatment))
# Prepare data for PCA, MDS and hierarchical clustering
## PCA
mesotaenium_study_design <- read_csv("assets/eda/mesotaenium_study_design.csv") %>%
    filter(condition != "constant_light_24")
mesotaenium_sampleLabels <- mesotaenium_study_design$sample_name
# create variables of important columns from study design file
mesotaenium_condition <- mesotaenium_study_design$condition
mesotaenium_condition <- factor(mesotaenium_condition, levels = unique(mesotaenium_study_design$condition))
mesotaenium_replicate <- mesotaenium_study_design$replicate
mesotaenium_replicate <- factor(mesotaenium_replicate, levels = unique(mesotaenium_study_design$replicate))
mesotaenium_treatment <- mesotaenium_study_design$treatment
mesotaenium_treatment <- factor(mesotaenium_treatment, levels = unique(mesotaenium_study_design$treatment))
mesotaenium_time_after_t0 <- mesotaenium_study_design$time
mesotaenium_temperature <- mesotaenium_study_design$temperature
mesotaenium_light_intensity <- mesotaenium_study_design$light_intensity
pca_data_mesotaenium <- mesotaenium_expr %>%
    select(-c(1, dim(mesotaenium_expr)[2]-2, dim(mesotaenium_expr)[2]-1, dim(mesotaenium_expr)[2]))
rownames(pca_data_mesotaenium) <- mesotaenium_expr$sample_name
pca.res_mesotaenium <- prcomp(pca_data_mesotaenium, scale.=F, retx=T)
pc.var_mesotaenium<-pca.res_mesotaenium$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var_mesotaenium/sum(pc.var_mesotaenium)*100, 1)
pca.res_mesotaenium.df <- as_tibble(pca.res_mesotaenium$x)
pca_summary_mesotaenium <- summary(pca.res_mesotaenium) # variance summary for all principal components.
pca_summary_mesotaenium.df <- as_tibble(pca_summary_mesotaenium$importance, rownames = "NA")
# Define a palette
my_colors_mesotaenium <- c("standard_growth_0" = "#e2e6e6", "standard_growth_0.5" = "#c6cdce", "standard_growth_1" = "#a8b4b5", "standard_growth_2" = "#8b9b9d", "standard_growth_4" = "#6e8284", "standard_growth_6" = "#51696c", "standard_growth_24" = "#345053", "constant_light_24" = "#FFFFFF", "cold_0.5" = "#91e0ef", "cold_1" = "#48cae4", "cold_2" = "#00b4d8", "cold_4" = "#0296c7", "cold_6" = "#0077b6", "cold_24" = "#023e8a", "heat_0.5" = "#fdc3c3", "heat_1" = "#feaeae", "heat_2" = "#fe9a9a", "heat_4" = "#fe8585", "heat_6" = "#fe5d5d", "heat_24" = "#fe4848", "high_light_s_0.25" = "#ffea01", "high_light_s_0.5" = "#ffdd02", "high_light_s_2" = "#ffd002", "high_light_s_4" = "#ffb700", "high_light_s_6" = "#ffa201", "high_light_r_0.25" = "#dee602", "high_light_r_0.5" = "#c3d205", "high_light_r_1" = "#a7be07", "high_light_r_2" = "#8daa0b", "high_light_r_4" = "#70960e")
## MDS
mds_mesotaenium_data_mesotaenium <- read_tsv("assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv")
mds_mesotaenium <- plotMDS(mds_mesotaenium_data_mesotaenium[,-1], plot = FALSE)
df_mesotaenium <- data.frame(Sample=colnames(mds_mesotaenium$distance.matrix.squared), mds_mesotaenium$x, mds_mesotaenium$y)
## Clustering
hierarchical_clustering_mesotaenium <- mds_mesotaenium_data_mesotaenium[,-1]
colnames(hierarchical_clustering_mesotaenium) <- mesotaenium_sampleLabels
distance_mesotaenium <- dist(t(hierarchical_clustering_mesotaenium), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters_mesotaenium <- hclust(distance_mesotaenium, method = "ward.D") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
cls_mesotaenium <- split(mesotaenium_study_design$sample_name, mesotaenium_study_design$condition)
clusters_mesotaenium_OTU <- groupOTU(as.phylo(clusters_mesotaenium), cls_mesotaenium)
## PCA
physcomitrium_study_design <- read_csv("assets/eda/physcomitrium_study_design.csv") %>%
    filter(condition != "constant_light_24")
physcomitrium_sampleLabels <- physcomitrium_study_design$sample_name
# create variables of important columns from study design file
physcomitrium_condition <- physcomitrium_study_design$condition
physcomitrium_condition <- factor(physcomitrium_condition, levels = unique(physcomitrium_study_design$condition))
physcomitrium_replicate <- physcomitrium_study_design$replicate
physcomitrium_replicate <- factor(physcomitrium_replicate, levels = unique(physcomitrium_study_design$replicate))
physcomitrium_treatment <- physcomitrium_study_design$treatment
physcomitrium_treatment <- factor(physcomitrium_treatment, levels = unique(physcomitrium_study_design$treatment))
physcomitrium_time_after_t0 <- physcomitrium_study_design$time
physcomitrium_temperature <- physcomitrium_study_design$temperature
physcomitrium_light_intensity <- physcomitrium_study_design$light_intensity
pca_data_physcomitrium <- physcomitrium_expr %>%
    select(-c(1, dim(physcomitrium_expr)[2]-2, dim(physcomitrium_expr)[2]-1, dim(physcomitrium_expr)[2]))
rownames(pca_data_physcomitrium) <- physcomitrium_expr$sample_name
pca.res_physcomitrium <- prcomp(pca_data_physcomitrium, scale.=F, retx=T)
pc.var_physcomitrium<-pca.res_physcomitrium$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var_physcomitrium/sum(pc.var_physcomitrium)*100, 1)
pca.res_physcomitrium.df <- as_tibble(pca.res_physcomitrium$x)
pca_summary_physcomitrium <- summary(pca.res_physcomitrium) # variance summary for all principal components.
pca_summary_physcomitrium.df <- as_tibble(pca_summary_physcomitrium$importance, rownames = "NA")
# Define a palette
my_colors_physcomitrium <- c("standard_growth_0" = "#e2e6e6",  "standard_growth_0.5" = "#c6cdce",  "standard_growth_1" = "#a8b4b5",  "standard_growth_2" = "#8b9b9d",  "standard_growth_4" = "#6e8284",  "standard_growth_6" = "#51696c",  "standard_growth_24" = "#345053",  "constant_light_24" = "#FFFFFF",  "cold_0.5" = "#91e0ef",  "cold_1" = "#48cae4",  "cold_2" = "#00b4d8",  "cold_4" = "#0296c7",  "cold_6" = "#0077b6",  "cold_24" = "#023e8a",  "heat_0.5" = "#fdc3c3",  "heat_1" = "#feaeae",  "heat_2" = "#fe9a9a",  "heat_4" = "#fe8585",  "heat_6" = "#fe5d5d",  "heat_24" = "#fe4848",  "high_light_s_0.25" = "#ffea01",  "high_light_s_0.5" = "#ffdd02",  "high_light_s_2" = "#ffd002",  "high_light_s_4" = "#ffb700",  "high_light_s_6" = "#ffa201",  "high_light_r_0.25" = "#dee602",  "high_light_r_0.5" = "#c3d205",  "high_light_r_1" = "#a7be07",  "high_light_r_2" = "#8daa0b",  "high_light_r_4" = "#70960e")
## MDS
mds_physcomitrium_data_physcomitrium <- read_tsv("assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv")
mds_physcomitrium <- plotMDS(mds_physcomitrium_data_physcomitrium[,-1], plot = FALSE)
df_physcomitrium <- data.frame(Sample=colnames(mds_physcomitrium$distance.matrix.squared), mds_physcomitrium$x, mds_physcomitrium$y)
## Clustering
hierarchical_clustering_physcomitrium <- mds_physcomitrium_data_physcomitrium[,-1]
colnames(hierarchical_clustering_physcomitrium) <- physcomitrium_sampleLabels
distance_physcomitrium <- dist(t(hierarchical_clustering_physcomitrium), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters_physcomitrium <- hclust(distance_physcomitrium, method = "ward.D") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
cls_physcomitrium <- split(physcomitrium_study_design$sample_name, physcomitrium_study_design$condition)
clusters_physcomitrium_OTU <- groupOTU(as.phylo(clusters_physcomitrium), cls_physcomitrium)
## PCA
zygnema_study_design <- read_csv("assets/eda/zygnema_study_design.csv") %>%
    filter(condition != "constant_light_24")
zygnema_sampleLabels <- zygnema_study_design$sample_name
# create variables of important columns from study design file
zygnema_condition <- zygnema_study_design$condition
zygnema_condition <- factor(zygnema_condition, levels = unique(zygnema_study_design$condition))
zygnema_replicate <- zygnema_study_design$replicate
zygnema_replicate <- factor(zygnema_replicate, levels = unique(zygnema_study_design$replicate))
zygnema_treatment <- zygnema_study_design$treatment
zygnema_treatment <- factor(zygnema_treatment, levels = unique(zygnema_study_design$treatment))
zygnema_time_after_t0 <- zygnema_study_design$time
zygnema_temperature <- zygnema_study_design$temperature
zygnema_light_intensity <- zygnema_study_design$light_intensity
pca_data_zygnema <- zygnema_expr %>%
    select(-c(1, dim(zygnema_expr)[2]-2, dim(zygnema_expr)[2]-1, dim(zygnema_expr)[2]))
rownames(pca_data_zygnema) <- zygnema_expr$sample_name
pca.res_zygnema <- prcomp(pca_data_zygnema, scale.=F, retx=T)
pc.var_zygnema<-pca.res_zygnema$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var_zygnema/sum(pc.var_zygnema)*100, 1)
pca.res_zygnema.df <- as_tibble(pca.res_zygnema$x)
pca_summary_zygnema <- summary(pca.res_zygnema) # variance summary for all principal components.
pca_summary_zygnema.df <- as_tibble(pca_summary_zygnema$importance, rownames = "NA")
# Define a palette
my_colors_zygnema <- c("standard_growth_0" = "#e2e6e6",  "standard_growth_0.5" = "#c6cdce",  "standard_growth_1" = "#a8b4b5",  "standard_growth_2" = "#8b9b9d",  "standard_growth_4" = "#6e8284",  "standard_growth_6" = "#51696c",  "standard_growth_24" = "#345053",  "constant_light_24" = "#FFFFFF",  "cold_0.5" = "#91e0ef",  "cold_1" = "#48cae4",  "cold_2" = "#00b4d8",  "cold_4" = "#0296c7",  "cold_6" = "#0077b6",  "cold_24" = "#023e8a",  "heat_0.5" = "#fdc3c3",  "heat_1" = "#feaeae",  "heat_2" = "#fe9a9a",  "heat_4" = "#fe8585",  "heat_6" = "#fe5d5d",  "heat_24" = "#fe4848",  "high_light_s_0.25" = "#ffea01",  "high_light_s_0.5" = "#ffdd02",  "high_light_s_2" = "#ffd002",  "high_light_s_4" = "#ffb700",  "high_light_s_6" = "#ffa201",  "high_light_r_0.25" = "#dee602",  "high_light_r_0.5" = "#c3d205",  "high_light_r_1" = "#a7be07",  "high_light_r_2" = "#8daa0b",  "high_light_r_4" = "#70960e"  )
## MDS
mds_zygnema_data_zygnema <- read_tsv("assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv")
mds_zygnema <- plotMDS(mds_zygnema_data_zygnema[,-1], plot = FALSE)
df_zygnema <- data.frame(Sample=colnames(mds_zygnema$distance.matrix.squared), mds_zygnema$x, mds_zygnema$y)
## Clustering
hierarchical_clustering_zygnema <- mds_zygnema_data_zygnema[,-1]
colnames(hierarchical_clustering_zygnema) <- zygnema_sampleLabels
distance_zygnema <- dist(t(hierarchical_clustering_zygnema), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters_zygnema <- hclust(distance_zygnema, method = "ward.D") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
cls_zygnema <- split(zygnema_study_design$sample_name, zygnema_study_design$condition)
clusters_zygnema_OTU <- groupOTU(as.phylo(clusters_zygnema), cls_zygnema)
myheatcolors <- c("#033270", "#1368aa", "#4091c9", "#9dcee2", "#FFFFFF", "#f29479", "#f26a4f", "#ef3c2d", "#cb1b16")
# Read the tree for the phylogeny tree
tree <- read.tree("assets/orthofinder/SpeciesTree.txt")
species_tree <- ggtree(tree, ladderize=FALSE, branch.length="none") %>%
    ggtree::rotate(26) %>%
    ggtree::rotate(33) %>%
    ggtree::rotate(34) %>%
    ggtree::rotate(39) +
    ggtree::geom_tiplab(as_ylab=TRUE,
                size =12)
################################################################################
############### UI
################################################################################
# This is my header
ui_header <- dashboardHeader(
    disable = FALSE,
    title = "Stress profiling",
    titleWidth = 320
)
# This is my sidebar
ui_sidebar <- dashboardSidebar(
    width = 320,
    sidebarMenu(
        menuItem("Home", tabName = "index"),
        menuItem("How to use this Shiny app?", tabName = "tutorial"),
        menuItem("Data", tabName = "data"),
        menuItem("Methods", tabName = "methods"),
        menuItem("Find your favourite gene(s) here", tabName = "functional_annotation"),
        menuItem("Exploratory data analysis", tabName = "eda"),
        menuItem("Gene expression visualization", tabName = "gene_profiler"),
        menuItem("Differential gene expression analysis", tabName = "dgea"),
        menuItem("Co-expression network analysis (DPGP)", tabName = "co_expression_dpgp"),
        menuItem("Co-expression network analysis (WGCNA)", tabName = "co_expression_wgcna"),
        menuItem("Predicted regulator-regulated pairs", tabName = "grn")
    )
)
# This is my body
ui_body <- dashboardBody(
    # Set the background overflow so the body background does not disappear after some scrolling
    tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }'))),
    # Arrange items
    tabItems(
        tabItem(tabName = "index",
                htmltools::includeMarkdown("markdowns/index.md")
                ),
        tabItem(tabName = "tutorial",
                htmltools::includeMarkdown("markdowns/tutorial.md"),
                HTML('<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/bZrEmW7mzuA?si=VjEScTG1dDuJi_FQ" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>')
        ),
        tabItem(tabName = "functional_annotation",
                box(
                    width = 12,
                    title = "Gene Discovery Toolbox",
                    htmltools::includeMarkdown("markdowns/functional_annotation.md"),
                    plotOutput("species_tree")
                ),
                box(
                    width = 12,
                    title = "Orthofinder",
                    tabBox(
                        width = 12,
                        title = "",
                        tabPanel("Phylogenetic Hierarchical Orthogroups",
                                 DT::dataTableOutput("orthofinder_n0")),
                        tabPanel("Multiple Sequence alignment viewer",
                                 textInput(inputId = "msa_viewer", label = "Write JUST 1 Orthogroup ID", value = "OG0003200"),
                                 msaROutput("msa_viewer_output", width="100%")
                                 )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "eggNOG-mapper",
                    tabBox(
                        width = 12,
                        title = "Species",
                        tabPanel("M. endlicharianum",
                                 DT::dataTableOutput("eggnog_mesotaenium")),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("eggnog_zygnema")),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("eggnog_physcomitrium"))
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "TapScan",
                    tabBox(
                        width = 12,
                        title = "Species",
                        tabPanel("M. endlicharianum",
                                 DT::dataTableOutput("tapscan_mesotaenium")),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("tapscan_zygnema")),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("tapscan_physcomitrium"))
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Gene Ontology",
                    tabBox(
                        width = 12,
                        title = "Species",
                        tabPanel("M. endlicharianum",
                                 DT::dataTableOutput("GO_mesotaenium")),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("GO_zygnema")),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("GO_physcomitrium"))
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Best blast hit against Araport11",
                    tabBox(
                        width = 12,
                        title = "Species",
                        tabPanel("M. endlicharianum",
                                 DT::dataTableOutput("blast_mesotaenium")),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("blast_zygnema")),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("blast_physcomitrium"))
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "InterProScan",
                    tabBox(
                        width = 12,
                        title = "Species",
                        tabPanel("M. endlicharianum",
                                 DT::dataTableOutput("interproscan_mesotaenium")),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("interproscan_zygnema")),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("interproscan_physcomitrium"))
                    ),
                    style = "overflow-x: scroll;"
                )
                ),
        tabItem(
            tabName = "eda",
            box(
                width = 12,
                title = "Exploratory data analysis",
                htmltools::includeMarkdown("markdowns/eda.md")
            ),
            box(
                width = 12,
                title = "Principal component analysis (PCA) plot",
                tabBox(
                    height = "650px",
                    width = "100%",
                    title = "Species",
                    tabPanel("M. endlicherianum",
                             plotOutput("mesotaenium_plot_PCA")),
                    tabPanel("Z. circumcarinatum",
                             plotOutput("zygnema_plot_PCA")),
                    tabPanel("P. patens",
                             plotOutput("physcomitrium_plot_PCA"))
                ),
                style = "overflow-x: scroll;"
            ),
            box(
                width = 12,
                title = "Multidimensional scaling (MDS) plot",
                tabBox(
                    height = "650px",
                    width = "100%",
                    title = "Species",
                    tabPanel("M. endlicherianum",
                             plotOutput("mesotaenium_plotMDS")),
                    tabPanel("Z. circumcarinatum",
                             plotOutput("zygnema_plotMDS")),
                    tabPanel("P. patens",
                             plotOutput("physcomitrium_plotMDS"))
                ),
                style = "overflow-x: scroll;"
            ),
            box(
                width = 12,
                title = "Hierarchical clustering",
                tabBox(
                    height = "1100px",
                    width = "100%",
                    title = "Species",
                    tabPanel("M. endlicherianum",
                             plotOutput("mesotaenium_dendrogram")),
                    tabPanel("Z. circumcarinatum",
                             plotOutput("zygnema_dendrogram")),
                    tabPanel("P. patens",
                             plotOutput("physcomitrium_dendrogram"))
                ),
                style = "overflow-x: scroll;"
            )
            ),
        tabItem(tabName = "gene_profiler",
                box(
                    width = 12,
                    title = "Gene expression visualization",
                    htmltools::includeMarkdown("markdowns/gene_visualization.md")
                ),
                box(
                    width = 12,
                    title = "Plots",
                    radioButtons(inputId = "gene_visualization_type",
                                 label = "Please choose a plot type",
                                 choices = c("Boxplot_conditions", "Boxplot_treatment", "Dotplot", "Heatmap"),
                                 selected = "Boxplot_conditions"),
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 textInput(inputId = "mesotaenium_gene_id_visualization",
                                           label = "Write JUST 1 gene ID",
                                           value = "Me1_v2_0054690"),
                                 plotOutput("mesotaenium_gene_visualization_plot")),
                        tabPanel("Z. circumcarinatum", 
                                 textInput(inputId = "zygnema_gene_id_visualization",
                                           label = "Write JUST 1 gene ID",
                                           value = "Zci_00086"),
                                 plotOutput("zygnema_gene_visualization_plot")),
                        tabPanel("P. patens", 
                                 textInput(inputId = "physcomitrium_gene_id_visualization",
                                           label = "Write JUST 1 gene ID",
                                           value = "Pp3c1_12790V3"),
                                 plotOutput("physcomitrium_gene_visualization_plot"))
                    ),
                    style = "overflow-x: scroll;"
                )
                ),
        tabItem(tabName = "dgea",
                box(
                    width = 8,
                    title = "Differential Gene Expression Analysis",
                    htmltools::includeMarkdown("markdowns/dgea.md")
                ),
                box(
                    width = 4,
                    title = "Parameters for visualization and filtering",
                    selectInput(inputId = "dgea_comparison", label = h4("Select a comparison"),
                                choices = c("mesotaenium_cold_0.5_vs_standard_growth_0.5", "mesotaenium_cold_1_vs_standard_growth_1", "mesotaenium_cold_24_vs_standard_growth_24", "mesotaenium_cold_2_vs_standard_growth_2", "mesotaenium_cold_4_vs_standard_growth_4", "mesotaenium_cold_6_vs_standard_growth_6", "mesotaenium_heat_0.5_vs_standard_growth_0.5", "mesotaenium_heat_24_vs_standard_growth_24", "mesotaenium_heat_2_vs_standard_growth_2", "mesotaenium_heat_4_vs_standard_growth_4", "mesotaenium_heat_6_vs_standard_growth_6", "mesotaenium_high_light_r_0.25_vs_high_light_s_6", "mesotaenium_high_light_r_0.25_vs_standard_growth_0", "mesotaenium_high_light_r_0.5_vs_high_light_s_6", "mesotaenium_high_light_r_0.5_vs_standard_growth_0", "mesotaenium_high_light_r_1_vs_high_light_s_6", "mesotaenium_high_light_r_1_vs_standard_growth_0", "mesotaenium_high_light_r_2_vs_high_light_s_6", "mesotaenium_high_light_r_2_vs_standard_growth_0", "mesotaenium_high_light_r_4_vs_high_light_s_6", "mesotaenium_high_light_r_4_vs_standard_growth_0", "mesotaenium_high_light_s_0.25_vs_standard_growth_0", "mesotaenium_high_light_s_0.5_vs_standard_growth_0.5", "mesotaenium_high_light_s_0.5_vs_standard_growth_0", "mesotaenium_high_light_s_2_vs_standard_growth_0", "mesotaenium_high_light_s_2_vs_standard_growth_2", "mesotaenium_high_light_s_4_vs_standard_growth_0", "mesotaenium_high_light_s_4_vs_standard_growth_4", "mesotaenium_high_light_s_6_vs_standard_growth_0", "mesotaenium_high_light_s_6_vs_standard_growth_6",
                                            "zygnema_cold_0.5_vs_standard_growth_0.5", "zygnema_cold_1_vs_standard_growth_1", "zygnema_cold_24_vs_standard_growth_24", "zygnema_cold_2_vs_standard_growth_2", "zygnema_cold_4_vs_standard_growth_4", "zygnema_cold_6_vs_standard_growth_6", "zygnema_heat_0.5_vs_standard_growth_0.5", "zygnema_heat_24_vs_standard_growth_24", "zygnema_heat_2_vs_standard_growth_2", "zygnema_heat_4_vs_standard_growth_4", "zygnema_heat_6_vs_standard_growth_6", "zygnema_high_light_r_0.25_vs_high_light_s_6", "zygnema_high_light_r_0.25_vs_standard_growth_0", "zygnema_high_light_r_0.5_vs_high_light_s_6", "zygnema_high_light_r_0.5_vs_standard_growth_0", "zygnema_high_light_r_1_vs_high_light_s_6", "zygnema_high_light_r_1_vs_standard_growth_0", "zygnema_high_light_r_2_vs_high_light_s_6", "zygnema_high_light_r_2_vs_standard_growth_0", "zygnema_high_light_r_4_vs_high_light_s_6", "zygnema_high_light_r_4_vs_standard_growth_0", "zygnema_high_light_s_0.25_vs_standard_growth_0", "zygnema_high_light_s_0.5_vs_standard_growth_0.5", "zygnema_high_light_s_0.5_vs_standard_growth_0", "zygnema_high_light_s_2_vs_standard_growth_0", "zygnema_high_light_s_2_vs_standard_growth_2", "zygnema_high_light_s_4_vs_standard_growth_0", "zygnema_high_light_s_4_vs_standard_growth_4", "zygnema_high_light_s_6_vs_standard_growth_0", "zygnema_high_light_s_6_vs_standard_growth_6",
                                            "physcomitrium_cold_0.5_vs_standard_growth_0.5", "physcomitrium_cold_1_vs_standard_growth_1", "physcomitrium_cold_24_vs_standard_growth_24", "physcomitrium_cold_2_vs_standard_growth_2", "physcomitrium_cold_4_vs_standard_growth_4", "physcomitrium_cold_6_vs_standard_growth_6", "physcomitrium_heat_0.5_vs_standard_growth_0.5", "physcomitrium_heat_24_vs_standard_growth_24", "physcomitrium_heat_2_vs_standard_growth_2", "physcomitrium_heat_4_vs_standard_growth_4", "physcomitrium_heat_6_vs_standard_growth_6", "physcomitrium_high_light_r_0.25_vs_high_light_s_6", "physcomitrium_high_light_r_0.25_vs_standard_growth_0", "physcomitrium_high_light_r_0.5_vs_high_light_s_6", "physcomitrium_high_light_r_0.5_vs_standard_growth_0", "physcomitrium_high_light_r_1_vs_high_light_s_6", "physcomitrium_high_light_r_1_vs_standard_growth_0", "physcomitrium_high_light_r_2_vs_high_light_s_6", "physcomitrium_high_light_r_2_vs_standard_growth_0", "physcomitrium_high_light_r_4_vs_high_light_s_6", "physcomitrium_high_light_r_4_vs_standard_growth_0", "physcomitrium_high_light_s_0.25_vs_standard_growth_0", "physcomitrium_high_light_s_0.5_vs_standard_growth_0.5", "physcomitrium_high_light_s_0.5_vs_standard_growth_0", "physcomitrium_high_light_s_2_vs_standard_growth_0", "physcomitrium_high_light_s_2_vs_standard_growth_2", "physcomitrium_high_light_s_4_vs_standard_growth_0", "physcomitrium_high_light_s_4_vs_standard_growth_4", "physcomitrium_high_light_s_6_vs_standard_growth_0", "physcomitrium_high_light_s_6_vs_standard_growth_6"),
                                selected = "physcomitrium_heat_2_vs_standard_growth_2"),
                    numericInput(inputId = "fold_change_threshold", label = h4("Fold change cutoff"), value = 2, min = 1, max = 20, step = 0.1),
                    numericInput(inputId = "adjusted_p_value", label = h4("Adjusted p-value cutoff"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                    selectInput(inputId = "go_domain", label = h4("Select a GO domain for ORA"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process"),
                    numericInput(inputId = "go_show_category", label = h4("The number of categories to show in ORA plots"), value = 20, min = 2, max = 100, step = 1),
                    numericInput(inputId = "go_pvalueCutoff", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                    numericInput(inputId = "go_qvalueCutoff", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                ),
                box(
                    width = 12,
                    title = "GO ORA results",
                    tabBox(
                        height = "1100px",
                        width = "100%",
                        title = "modules",
                        tabPanel("Plot of GO enrichment analysis (ORA): Up-regulated",
                                 plotOutput("go_ora_up_regulated")
                        ),
                        tabPanel("Plot of GO enrichment analysis (ORA): Down-regulated",
                                 plotOutput("go_ora_down_regulated")
                        ),
                        tabPanel("Table of GO enrichment analysis (ORA): Up-regulated",
                                 DT::dataTableOutput("go_enrichment_up_regulated_table")
                        ),
                        tabPanel("Table of GO enrichment analysis (ORA): Down-regulated",
                                 DT::dataTableOutput("go_enrichment_down_regulated_table")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Table of Log2(Fold change) and adjusted P-value compared to the control",
                    DT::dataTableOutput("logfc")
                ),
                box(
                    width = 12,
                    title = "Heatmaps of Z-score transformed of qsmooth normalized and voom transformed counts including only DEGs",
                    tabBox(
                        height = "900px",
                        width = "100%",
                        title = "Gene expression profiles",
                        tabPanel("Heatmap of up-regulated",
                                 plotOutput("expression_plots_up_regulated")
                        ),
                        tabPanel("Heatmap of down-regulated",
                                 plotOutput("expression_plots_down_regulated")
                        ),
                        tabPanel("Heatmap of DEGs (Z-score transformed)",
                                   plotOutput("expression_plots_degs")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Tables of qsmooth normalized and voom transformed counts including only DEGs",
                    tabBox(
                        width = "100%",
                        title = "Gene expression profiles",
                        tabPanel("DEGs (Z-score transformed)",
                                 DT::dataTableOutput("expression_tables_degs")
                        ),
                        tabPanel("Up-regulated DEGs (Z-score transformed)",
                                 DT::dataTableOutput("expression_tables_up_regulated")
                        ),
                        tabPanel("Down-regulated DEGs (Z-score transformed)",
                                 DT::dataTableOutput("expression_tables_down_regulated")
                        )
                    ),
                    style = "overflow-x: scroll;"
                )
                ),
        tabItem(tabName = "co_expression_dpgp",
                box(
                    width = 6,
                    title = "Co-expression network analysis using DPGP",
                    htmltools::includeMarkdown("markdowns/co_expression_dpgp.md")
                ),
                box(
                    width = 6,
                    title = "Pick parameters for analysis and visualization of DPGP results",
                    tabBox(
                        width = "100%",
                        tabPanel("M. endlicherianum",
                                 selectInput(inputId = "dpgp_input_mesotaenium", label = h4("Select a network"),
                                             choices = c("mesotaenium_cold_optimal_clustering", "mesotaenium_heat_optimal_clustering", "mesotaenium_high_light_no_recovery_optimal_clustering"),
                                             selected = "mesotaenium_high_light_no_recovery_optimal_clustering"),
                                 numericInput(inputId = "dpgp_probability_filter_threshold_mesotaenium", label = h4("Threshold of probability inclusion (If you want to have more 'conserved' clusters increase number towards 1)"), value = 0, min = 0, max = 1, step = 0.1),
                                 uiOutput("dpgp_mesotaenium_dynamicInput"),
                                 selectInput(inputId = "go_domain_dpgp_mesotaenium", label = h4("Select a GO domain"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process"),
                                 numericInput(inputId = "go_show_category_dpgp_mesotaenium", label = h4("The number of categories to show"), value = 20, min = 2, max = 100, step = 1),
                                 numericInput(inputId = "go_pvalueCutoff_dpgp_mesotaenium", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                                 numericInput(inputId = "go_qvalueCutoff_dpgp_mesotaenium", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                                 
                        ),
                        tabPanel("Z. circumcarinatum", 
                                 selectInput(inputId = "dpgp_input_zygnema", label = h4("Select a network"),
                                             choices = c("zygnema_cold_optimal_clustering", "zygnema_heat_optimal_clustering", "zygnema_high_light_optimal_clustering"),
                                             selected = "zygnema_high_light_optimal_clustering"),
                                 numericInput(inputId = "dpgp_probability_filter_threshold_zygnema", label = h4("Threshold of probability inclusion (If you want to have more 'conserved' clusters increase number towards 1)"), value = 0, min = 0, max = 1, step = 0.1),
                                 uiOutput("dpgp_zygnema_dynamicInput"),
                                 selectInput(inputId = "go_domain_dpgp_zygnema", label = h4("Select a GO domain"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process"),
                                 numericInput(inputId = "go_show_category_dpgp_zygnema", label = h4("The number of categories to show"), value = 20, min = 2, max = 100, step = 1),
                                 numericInput(inputId = "go_pvalueCutoff_dpgp_zygnema", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                                 numericInput(inputId = "go_qvalueCutoff_dpgp_zygnema", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                        ),
                        tabPanel("P. patens", 
                                 selectInput(inputId = "dpgp_input_physcomitrium", label = h4("Select a network"),
                                             choices = c("physcomitrium_cold_optimal_clustering", "physcomitrium_heat_optimal_clustering", "physcomitrium_high_light_optimal_clustering"),
                                             selected = "physcomitrium_high_light_optimal_clustering"),
                                 numericInput(inputId = "dpgp_probability_filter_threshold_physcomitrium", label = h4("Threshold of probability inclusion (If you want to have more 'conserved' clusters increase number towards 1)"), value = 0, min = 0, max = 1, step = 0.1),
                                 uiOutput("dpgp_physcomitrium_dynamicInput"),
                                 selectInput(inputId = "go_domain_dpgp_physcomitrium", label = h4("Select a GO domain"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process"),
                                 numericInput(inputId = "go_show_category_dpgp_physcomitrium", label = h4("The number of categories to show"), value = 20, min = 2, max = 100, step = 1),
                                 numericInput(inputId = "go_pvalueCutoff_dpgp_physcomitrium", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                                 numericInput(inputId = "go_qvalueCutoff_dpgp_physcomitrium", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                        )
                    )
                ),
                box(
                    width = 6,
                    title = "DPGP cluster expression visualization",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 plotOutput("mesotaenium_dpgp_expression_visualization")
                        ),
                        tabPanel("Z. circumcarinatum", 
                                 plotOutput("zygnema_dpgp_expression_visualization")
                        ),
                        tabPanel("P. patens", 
                                 plotOutput("physcomitrium_dpgp_expression_visualization")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 6,
                    title = "DPGP cluster GO ORA visualization",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 plotOutput("mesotaenium_dpgp_ora_visualization")
                        ),
                        tabPanel("Z. circumcarinatum", 
                                 plotOutput("zygnema_dpgp_ora_visualization")
                        ),
                        tabPanel("P. patens", 
                                 plotOutput("physcomitrium_dpgp_ora_visualization")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "DPGP cluster GO ORA tables",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 DT::dataTableOutput("mesotaenium_dpgp_ora_table")
                        ),
                        tabPanel("Z. circumcarinatum",
                                 uiOutput("grn_zygnema_dynamicInput"),
                                 DT::dataTableOutput("zygnema_dpgp_ora_table")
                        ),
                        tabPanel("P. patens", 
                                 uiOutput("grn_physcomitrium_dynamicInput"),
                                 DT::dataTableOutput("physcomitrium_dpgp_ora_table")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Table of genes clustered via DPGP (all clusters)",
                    tabBox(
                        height = "600px",
                        width = "100%",
                        tabPanel("M. endlicherianum",
                                 DT::dataTableOutput("mesotaenium_dpgp_table")),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("zygnema_dpgp_table")),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("physcomitrium_dpgp_table"))
                        ),
                    style = "overflow-x: scroll;"
                )
                ),
        tabItem(tabName = "co_expression_wgcna",
                box(
                    width = 6,
                    title = "Co-expression network analysis WGCNA",
                    htmltools::includeMarkdown("markdowns/co_expression_wgcna.md")
                ),
                box(
                    width = 6,
                    title = "Pick your parameters for WGCNA visualization",
                    tabBox(
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 selectInput(inputId = "cluster_wgcna_mesotaenium", label = h4("Select a cluster and treatment for ORA analysis and visualization"), choices = c("black", "blue", "brown", "cyan", "darkgreen", "darkgrey", "darkorange", "darkred", "darkturquoise", "green", "greenyellow", "grey", "grey60", "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange", "pink", "purple", "red", "royalblue", "salmon", "tan", "turquoise", "yellow"), selected = "turquoise" ) , 
                                 selectInput(inputId = "go_domain_wgcna_mesotaenium", label = h4("Select a GO domain"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process") ,
                                 numericInput(inputId = "go_show_category_wgcna_mesotaenium", label = h4("The number of categories to show"), value = 20, min = 2, max = 100, step = 1),
                                 numericInput(inputId = "go_pvalueCutoff_wgcna_mesotaenium", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                                 numericInput(inputId = "go_qvalueCutoff_wgcna_mesotaenium", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                                 ),
                        tabPanel("Z. circumcarinatum",
                                 selectInput(inputId = "cluster_wgcna_zygnema", label = h4("Select a cluster and treatment for ORA analysis and visualization"), choices = c("blue", "black", "brown", "cyan", "darkgreen", "darkgrey", "darkorange", "darkred", "darkturquoise", "green", "greenyellow", "grey", "grey60", "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange", "pink", "purple", "red", "royalblue", "salmon", "skyblue", "tan", "turquoise", "white", "yellow"), selected =  "turquoise" ) , 
                                 selectInput(inputId = "go_domain_wgcna_zygnema", label = h4("Select a GO domain"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process") ,
                                 numericInput(inputId = "go_show_category_wgcna_zygnema", label = h4("The number of categories to show"), value = 20, min = 2, max = 100, step = 1),
                                 numericInput(inputId = "go_pvalueCutoff_wgcna_zygnema", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                                 numericInput(inputId = "go_qvalueCutoff_wgcna_zygnema", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                                 ),
                        tabPanel("P. patens",
                                 selectInput(inputId = "cluster_wgcna_physcomitrium", label = h4("Select a cluster and treatment for ORA analysis and visualization"), choices = c("cyan", "black", "blue", "brown", "darkgreen", "darkgrey", "darkolivegreen", "darkorange", "darkred", "darkturquoise", "green", "greenyellow", "grey", "grey60", "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange", "paleturquoise", "pink", "purple", "red", "royalblue", "saddlebrown", "salmon", "skyblue", "steelblue", "tan", "turquoise", "violet", "white", "yellow"), selected = "brown" ) , 
                                 selectInput(inputId = "go_domain_wgcna_physcomitrium", label = h4("Select a GO domain"), choices = list("Biological Process" = "biological_process", "Molecular function" = "molecular_function", "Cellular component" = "cellular_component"), selected = "biological_process") ,
                                 numericInput(inputId = "go_show_category_wgcna_physcomitrium", label = h4("The number of categories to show"), value = 20, min = 2, max = 100, step = 1),
                                 numericInput(inputId = "go_pvalueCutoff_wgcna_physcomitrium", label = h4("p-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01),
                                 numericInput(inputId = "go_qvalueCutoff_wgcna_physcomitrium", label = h4("q-value cutoff for ORA"), value = 0.05, min = 0.01, max = 0.3, step = 0.01)
                                 )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Summary table of WGCNA analysis",
                    tabBox(
                        height = "600px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 DT::dataTableOutput("mesotaenium_wgcna_table")
                                 ),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("zygnema_wgcna_table")
                                 ),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("physcomitrium_wgcna_table")
                                 )
                        ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "WGCNA Gene Significance visualization",
                    selectInput(inputId = "wgcna_gene_significance_type", label = h4("Gene Significance or absolute value of value of Gene Significance?"), choices = c("absolute value", "value"), selected = "value"),
                    selectInput(inputId = "wgcna_gene_significance_measure", label = h4("Pick a measurement for Gene Significance"), choices = c("beta_car_per_9_cis_neox", "A_Z_per_V_A_Z", "alpha_car_micromol_per_mgDW", "anthera_micromol_per_mgDW", "beta_car_micromol_per_mgDW", "beta_Car_per_beta_CC_beta_Io_DHA", "beta_Car_per_beta_CC", "beta_Car_per_beta_Io", "beta_Car_per_DHA", "beta_CC_nmol_per_mgDW", "beta_Io_nmol_per_mgDW", "Chla_micromol_per_mgDW", "Chla_per_b", "Chlb_micromol_per_mgDW", "DHA_nmol_per_mgDW", "DHA_per_beta_Io", "FvFm", "light_intensity", "Lut_micromol_per_mgDW", "replicate", "temperature", "time", "V_A_Z_per_Chla_b", "viol_micromol_per_mgDW", "X6MHO_nmol_per_mgDW", "X9_cis_micromol_per_mgDW", "X9_cis_neo_micromol_per_mgDW", "X11_cis_micromol_per_mgDW", "Zeax_micromol_per_mgDW"), selected = "light_intensity"),
                    tabBox(
                        height = "1050px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 plotOutput("mesotaenium_wgcna_GS_visualization")
                        ),
                        tabPanel("Z. circumcarinatum", 
                                 plotOutput("zygnema_wgcna_GS_visualization")
                        ),
                        tabPanel("P. patens", 
                                 plotOutput("physcomitrium_wgcna_GS_visualization")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "WGCNA cluster expression visualization",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 imageOutput("mesotaenium_wgcna_expression_visualization")
                        ),
                        tabPanel("Z. circumcarinatum", 
                                 imageOutput("zygnema_wgcna_expression_visualization")
                        ),
                        tabPanel("P. patens", 
                                 imageOutput("physcomitrium_wgcna_expression_visualization")
                        )
                    ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 6,
                    title = "Table of top 20 hubs for each module of WGCNA analysis",
                    tabBox(
                        height = "950px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 DT::dataTableOutput("mesotaenium_wgcna_top_hubs")
                                 ),
                        tabPanel("Z. circumcarinatum",
                                 DT::dataTableOutput("zygnema_wgcna_top_hubs")
                                 ),
                        tabPanel("P. patens",
                                 DT::dataTableOutput("physcomitrium_wgcna_top_hubs")
                                 )
                        ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 6,
                    title = "WGCNA cluster GO ORA visualization",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 plotOutput("mesotaenium_wgcna_ora_visualization")
                                 ),
                        tabPanel("Z. circumcarinatum", 
                                 plotOutput("zygnema_wgcna_ora_visualization")
                                 ),
                        tabPanel("P. patens", 
                                 plotOutput("physcomitrium_wgcna_ora_visualization")
                                 )
                        ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "WGCNA cluster GO ORA tables",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 DT::dataTableOutput("mesotaenium_wgcna_ora_table")
                                 ),
                        tabPanel("Z. circumcarinatum", 
                                 DT::dataTableOutput("zygnema_wgcna_ora_table")
                                 ),
                        tabPanel("P. patens", 
                                 DT::dataTableOutput("physcomitrium_wgcna_ora_table")
                                 )
                        ),
                    style = "overflow-x: scroll;"
                )
                ),
        tabItem(tabName = "grn",
                box(
                    width = 12,
                    title = "Gene regulatory networks",
                    htmltools::includeMarkdown("markdowns/gene_regulatory_network.md")
                ),
                box(
                    width = 12,
                    title = "Gene regulatory network tables",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 selectInput(inputId = "grn_input_mesotaenium", label = h4("Select a network"),
                                             choices = c("mesotaenium_with_metabolites_highlight", "mesotaenium_cold_only_transcripts", "mesotaenium_heat_only_transcripts", "mesotaenium_HL_only_transcripts", "mesotaenium_with_metabolites_cold", "mesotaenium_with_metabolites_heat"),
                                             selected = "mesotaenium_with_metabolites_highlight"),
                                 downloadButton("download_mesotaenium_grn_table", "Download complete table"),
                                 sliderInput(inputId = "grn_top_n_mesotaenium", label = h4("Show top N interactions (the maximum number on the website is 2000):"), value = 700, min = 1, max = 2000, step = 1),
                                 textInput(inputId = "grn_gene_filter_mesotaenium", label = "Write JUST 1 gene ID. Or, leave it empty.", value = ""),
                                 selectInput(inputId = "grn_filter_gene_category_mesotaenium", label = h4("Filter for this gene in:"), choices = list("Regulators" = "regulator", "Regulated" = "target", "Both" = "both"), selected = "both"),
                                 DT::dataTableOutput("mesotaenium_grn_table")
                                 ),
                        tabPanel("Z. circumcarinatum", 
                                 selectInput(inputId = "grn_input_zygnema", label = h4("Select a network"),
                                             choices = c("zygnema_only_transcripts_cold", "zygnema_only_transcripts_heat", "zygnema_only_transcripts_highlight", "zygnema_with_metabolites_cold", "zygnema_with_metabolites_heat", "zygnema_with_metabolites_highlight"),
                                             selected = "zygnema_only_transcripts_highlight"),
                                 downloadButton("download_zygnema_grn_table", "Download complete table"),
                                 sliderInput(inputId = "grn_top_n_zygnema", label = h4("Show top N interactions (the maximum number on the website is 2000)"), value = 700, min = 1, max = 2000, step = 1),
                                 textInput(inputId = "grn_gene_filter_zygnema", label = "Write JUST 1 gene ID. Or, leave it empty.", value = ""),
                                 selectInput(inputId = "grn_filter_gene_category_zygnema", label = h4("Filter for this gene in:"), choices = list("Regulators" = "regulator", "Regulated" = "target", "Both" = "both"), selected = "both"),
                                 DT::dataTableOutput("zygnema_grn_table")
                                 ),
                        tabPanel("P. patens", 
                                 selectInput(inputId = "grn_input_physcomitrium", label = h4("Select a network"),
                                             choices = c("physcomitrium_only_transcripts_cold", "physcomitrium_only_transcripts_heat", "physcomitrium_only_transcripts_highlight", "physcomitrium_with_metabolites_cold", "physcomitrium_with_metabolites_heat", "physcomitrium_with_metabolites_highlight"),
                                             selected = "physcomitrium_only_transcripts_highlight"),
                                 downloadButton("download_physcomitrium_grn_table", "Download complete table"),
                                 sliderInput(inputId = "grn_top_n_physcomitrium", label = h4("Show top N interactions (the maximum number on the website is 2000)"), value = 700, min = 1, max = 2000, step = 1),
                                 textInput(inputId = "grn_gene_filter_physcomitrium", label = "Write JUST 1 gene ID. Or, leave it empty.", value = ""),
                                 selectInput(inputId = "grn_filter_gene_category_physcomitrium", label = h4("Filter for this gene in:"), choices = list("Regulators" = "regulator", "Regulated" = "target", "Both" = "both"), selected = "both"),
                                 DT::dataTableOutput("physcomitrium_grn_table")
                                 )
                        ),
                    style = "overflow-x: scroll;"
                ),
                box(
                    width = 12,
                    title = "Gene regulatory network visualization",
                    tabBox(
                        height = "1000px",
                        width = "100%",
                        tabPanel("M. endlicherianum ",
                                 visNetworkOutput("mesotaenium_grn_visualization", height = "900px")
                                 ),
                        tabPanel("Z. circumcarinatum", 
                                 visNetworkOutput("zygnema_grn_visualization", height = "900px")
                                 ),
                        tabPanel("P. patens", 
                                 visNetworkOutput("physcomitrium_grn_visualization", height = "900px")
                                 )
                        ),
                    style = "overflow-x: scroll;"
                )
                ),
        tabItem(tabName = "data",
                width = 12,
                htmltools::includeMarkdown("markdowns/data.md")
                ),
        tabItem(tabName = "methods",
                width = 12,
                htmltools::includeMarkdown("markdowns/method.md")
                )
        )
)
# Bring everything together for the UI
ui <- fluidPage(
    # Include Roboto font from Google Fonts and specify the font size
    tags$head(
        tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Roboto"),
        tags$style(
            HTML("
      /* Specify Roboto as the font family for all text elements */
      body, input, select, textarea {
        font-family: 'Roboto', sans-serif !important;
        font-size: 16px !important; /* Body text size */
      }
      
      /* Specify font sizes for different elements */
      h1 { font-size: 30px !important; } /* Heading 1 */
      h2 { font-size: 24px !important; } /* Heading 2 */
      h3 { font-size: 20px !important; } /* Heading 3 */
      .sidebar .control-label { font-size: 16px !important; } /* Sidebar labels */
      .main-panel .control-label { font-size: 16px !important; } /* Main panel labels */
      ")
      )
    ),
    dashboardPage(header = ui_header,
                  sidebar = ui_sidebar,
                  body = ui_body,
                  title = "Stress profiling",
                  skin = "black"
    )
)
################################################################################
############### Server
################################################################################
server <- function(input, output, session) {
    # Functional annotation part
    ## Species tree
    output$species_tree <- renderPlot({species_tree}, width = 400, height = 400)
    ## InterProScan tables
    output$interproscan_mesotaenium <- DT::renderDataTable({datatable(interproscan_mesotaenium,
                                                                  options = list(paging = TRUE,    ## paginate the output
                                                                                 pageLength = 5,  ## number of rows to output for each page
                                                                                 lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                   c("5", "10", "25", "50", "100", 'All')),
                                                                                 scrollX = TRUE,   ## enable scrolling on X axis
                                                                                 scrollY = "500px",   ## enable scrolling on Y axis
                                                                                 rowCallback = JS(
                                                                                     "function(row, data, index) {",
                                                                                     "  $('td', row).css('vertical-align', 'top');",
                                                                                     "}"
                                                                                 ),
                                                                                 autoWidth = FALSE, ## use smart column width handling
                                                                                 server = TRUE,   ## use client-side processing
                                                                                 processing = TRUE,
                                                                                 dom = 'lfrtBip',
                                                                                 buttons = c('csv', 'excel')
                                                                  ),
                                                                  extensions = 'Buttons',
                                                                  selection = 'single', ## enable selection of a single row
                                                                  filter = 'top',              ## include column filters at the top
                                                                  rownames = FALSE                ## don't show row numbers/names
                                                                  )
        })
    output$interproscan_physcomitrium <- DT::renderDataTable({datatable(interproscan_physcomitrium,
                                                                        options = list(paging = TRUE,    ## paginate the output
                                                                                       pageLength = 5,  ## number of rows to output for each page
                                                                                       lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                         c("5", "10", "25", "50", "100", 'All')),
                                                                                       scrollX = TRUE,   ## enable scrolling on X axis
                                                                                       scrollY = "500px",   ## enable scrolling on Y axis
                                                                                       rowCallback = JS(
                                                                                           "function(row, data, index) {",
                                                                                           "  $('td', row).css('vertical-align', 'top');",
                                                                                           "}"
                                                                                       ),
                                                                                       autoWidth = FALSE, ## use smart column width handling
                                                                                       server = TRUE,   ## use client-side processing
                                                                                       processing = TRUE,
                                                                                       dom = 'lfrtBip',
                                                                                       buttons = c('csv', 'excel')
                                                                        ),
                                                                        extensions = 'Buttons',
                                                                        selection = 'single', ## enable selection of a single row
                                                                        filter = 'top',              ## include column filters at the top
                                                                        rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$interproscan_zygnema <- DT::renderDataTable({datatable(interproscan_zygnema,
                                                                  options = list(paging = TRUE,    ## paginate the output
                                                                                 pageLength = 5,  ## number of rows to output for each page
                                                                                 lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                   c("5", "10", "25", "50", "100", 'All')),
                                                                                 scrollX = TRUE,   ## enable scrolling on X axis
                                                                                 scrollY = "500px",   ## enable scrolling on Y axis
                                                                                 rowCallback = JS(
                                                                                     "function(row, data, index) {",
                                                                                     "  $('td', row).css('vertical-align', 'top');",
                                                                                     "}"
                                                                                 ),
                                                                                 autoWidth = FALSE, ## use smart column width handling
                                                                                 server = TRUE,   ## use client-side processing
                                                                                 processing = TRUE,
                                                                                 dom = 'lfrtBip',
                                                                                 buttons = c('csv', 'excel')
                                                                  ),
                                                                  extensions = 'Buttons',
                                                                  selection = 'single', ## enable selection of a single row
                                                                  filter = 'top',              ## include column filters at the top
                                                                  rownames = FALSE                ## don't show row numbers/names
    )
    })
    ## eggNOG-mapper tables
    output$eggnog_mesotaenium <- DT::renderDataTable({datatable(eggnog_mesotaenium,
                                                                options = list(paging = TRUE,    ## paginate the output
                                                                               pageLength = 25,  ## number of rows to output for each page
                                                                               lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                 c("5", "10", "25", "50", "100", 'All')),
                                                                               scrollX = TRUE,   ## enable scrolling on X axis
                                                                               scrollY = "500px",   ## enable scrolling on Y axis
                                                                               rowCallback = JS(
                                                                                   "function(row, data, index) {",
                                                                                   "  $('td', row).css('vertical-align', 'top');",
                                                                                   "}"
                                                                               ),
                                                                               autoWidth = FALSE, ## use smart column width handling
                                                                               server = TRUE,   ## use client-side processing
                                                                               processing = TRUE,
                                                                               dom = 'lfrtBip',
                                                                               buttons = c('csv', 'excel')
                                                                ),
                                                                extensions = 'Buttons',
                                                                selection = 'single', ## enable selection of a single row
                                                                filter = 'top',              ## include column filters at the top
                                                                rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$eggnog_physcomitrium <- DT::renderDataTable({datatable(eggnog_physcomitrium,
                                                                  options = list(paging = TRUE,    ## paginate the output
                                                                                 pageLength = 25,  ## number of rows to output for each page
                                                                                 lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                   c("5", "10", "25", "50", "100", 'All')),
                                                                                 scrollX = TRUE,   ## enable scrolling on X axis
                                                                                 scrollY = "500px",   ## enable scrolling on Y axis
                                                                                 rowCallback = JS(
                                                                                     "function(row, data, index) {",
                                                                                     "  $('td', row).css('vertical-align', 'top');",
                                                                                     "}"
                                                                                 ),
                                                                                 autoWidth = FALSE, ## use smart column width handling
                                                                                 server = TRUE,   ## use client-side processing
                                                                                 processing = TRUE,
                                                                                 dom = 'lfrtBip',
                                                                                 buttons = c('csv', 'excel')
                                                                  ),
                                                                  extensions = 'Buttons',
                                                                  selection = 'single', ## enable selection of a single row
                                                                  filter = 'top',              ## include column filters at the top
                                                                  rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$eggnog_zygnema <- DT::renderDataTable({datatable(eggnog_zygnema,
                                                            options = list(paging = TRUE,    ## paginate the output
                                                                           pageLength = 25,  ## number of rows to output for each page
                                                                           lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                             c("5", "10", "25", "50", "100", 'All')),
                                                                           scrollX = TRUE,   ## enable scrolling on X axis
                                                                           scrollY = "500px",   ## enable scrolling on Y axis
                                                                           rowCallback = JS(
                                                                               "function(row, data, index) {",
                                                                               "  $('td', row).css('vertical-align', 'top');",
                                                                               "}"
                                                                           ),
                                                                           autoWidth = FALSE, ## use smart column width handling
                                                                           server = TRUE,   ## use client-side processing
                                                                           processing = TRUE,
                                                                           dom = 'lfrtBip',
                                                                           buttons = c('csv', 'excel')
                                                            ),
                                                            extensions = 'Buttons',
                                                            selection = 'single', ## enable selection of a single row
                                                            filter = 'top',              ## include column filters at the top
                                                            rownames = FALSE                ## don't show row numbers/names
    )
    })
    ## tapscan tables
    output$tapscan_mesotaenium <- DT::renderDataTable({datatable(tapscan_mesotaenium,
                                                                 options = list(paging = TRUE,    ## paginate the output
                                                                                pageLength = 25,  ## number of rows to output for each page
                                                                                lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                  c("5", "10", "25", "50", "100", 'All')),
                                                                                scrollX = TRUE,   ## enable scrolling on X axis
                                                                                scrollY = "500px",   ## enable scrolling on Y axis
                                                                                rowCallback = JS(
                                                                                    "function(row, data, index) {",
                                                                                    "  $('td', row).css('vertical-align', 'top');",
                                                                                    "}"
                                                                                ),
                                                                                autoWidth = FALSE, ## use smart column width handling
                                                                                server = TRUE,   ## use client-side processing
                                                                                processing = TRUE,
                                                                                dom = 'lfrtBip',
                                                                                buttons = c('csv', 'excel')
                                                                 ),
                                                                 extensions = 'Buttons',
                                                                 selection = 'single', ## enable selection of a single row
                                                                 filter = 'top',              ## include column filters at the top
                                                                 rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$tapscan_physcomitrium <- DT::renderDataTable({datatable(tapscan_physcomitrium,
                                                                   options = list(paging = TRUE,    ## paginate the output
                                                                                  pageLength = 25,  ## number of rows to output for each page
                                                                                  lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                    c("5", "10", "25", "50", "100", 'All')),
                                                                                  scrollX = TRUE,   ## enable scrolling on X axis
                                                                                  scrollY = "500px",   ## enable scrolling on Y axis
                                                                                  rowCallback = JS(
                                                                                      "function(row, data, index) {",
                                                                                      "  $('td', row).css('vertical-align', 'top');",
                                                                                      "}"
                                                                                  ),
                                                                                  autoWidth = FALSE, ## use smart column width handling
                                                                                  server = TRUE,   ## use client-side processing
                                                                                  processing = TRUE,
                                                                                  dom = 'lfrtBip',
                                                                                  buttons = c('csv', 'excel')
                                                                   ),
                                                                   extensions = 'Buttons',
                                                                   selection = 'single', ## enable selection of a single row
                                                                   filter = 'top',              ## include column filters at the top
                                                                   rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$tapscan_zygnema <- DT::renderDataTable({datatable(tapscan_zygnema,
                                                             options = list(paging = TRUE,    ## paginate the output
                                                                            pageLength = 25,  ## number of rows to output for each page
                                                                            lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                              c("5", "10", "25", "50", "100", 'All')),
                                                                            scrollX = TRUE,   ## enable scrolling on X axis
                                                                            scrollY = "500px",   ## enable scrolling on Y axis
                                                                            rowCallback = JS(
                                                                                "function(row, data, index) {",
                                                                                "  $('td', row).css('vertical-align', 'top');",
                                                                                "}"
                                                                            ),
                                                                            autoWidth = FALSE, ## use smart column width handling
                                                                            server = TRUE,   ## use client-side processing
                                                                            processing = TRUE,
                                                                            dom = 'lfrtBip',
                                                                            buttons = c('csv', 'excel')
                                                             ),
                                                             extensions = 'Buttons',
                                                             selection = 'single', ## enable selection of a single row
                                                             filter = 'top',              ## include column filters at the top
                                                             rownames = FALSE                ## don't show row numbers/names
    )
    })
    ## GO tables
    output$GO_mesotaenium <- DT::renderDataTable({datatable(GO_mesotaenium,
                                                            options = list(paging = TRUE,    ## paginate the output
                                                                           pageLength = 25,  ## number of rows to output for each page
                                                                           lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                             c("5", "10", "25", "50", "100", 'All')),
                                                                           scrollX = TRUE,   ## enable scrolling on X axis
                                                                           scrollY = "500px",   ## enable scrolling on Y axis
                                                                           rowCallback = JS(
                                                                               "function(row, data, index) {",
                                                                               "  $('td', row).css('vertical-align', 'top');",
                                                                               "}"
                                                                           ),
                                                                           autoWidth = FALSE, ## use smart column width handling
                                                                           server = TRUE,   ## use client-side processing
                                                                           processing = TRUE,
                                                                           dom = 'lfrtBip',
                                                                           buttons = c('csv', 'excel')
                                                            ),
                                                            extensions = 'Buttons',
                                                            selection = 'single', ## enable selection of a single row
                                                            filter = 'top',              ## include column filters at the top
                                                            rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$GO_physcomitrium <- DT::renderDataTable({datatable(GO_physcomitrium,
                                                              options = list(paging = TRUE,    ## paginate the output
                                                                             pageLength = 25,  ## number of rows to output for each page
                                                                             lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                               c("5", "10", "25", "50", "100", 'All')),
                                                                             scrollX = TRUE,   ## enable scrolling on X axis
                                                                             scrollY = "500px",   ## enable scrolling on Y axis
                                                                             rowCallback = JS(
                                                                                 "function(row, data, index) {",
                                                                                 "  $('td', row).css('vertical-align', 'top');",
                                                                                 "}"
                                                                             ),
                                                                             autoWidth = FALSE, ## use smart column width handling
                                                                             server = TRUE,   ## use client-side processing
                                                                             processing = TRUE,
                                                                             dom = 'lfrtBip',
                                                                             buttons = c('csv', 'excel')
                                                              ),
                                                              extensions = 'Buttons',
                                                              selection = 'single', ## enable selection of a single row
                                                              filter = 'top',              ## include column filters at the top
                                                              rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$GO_zygnema <- DT::renderDataTable({datatable(GO_zygnema,
                                                        options = list(paging = TRUE,    ## paginate the output
                                                                       pageLength = 25,  ## number of rows to output for each page
                                                                       lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                         c("5", "10", "25", "50", "100", 'All')),
                                                                       scrollX = TRUE,   ## enable scrolling on X axis
                                                                       scrollY = "500px",   ## enable scrolling on Y axis
                                                                       rowCallback = JS(
                                                                           "function(row, data, index) {",
                                                                           "  $('td', row).css('vertical-align', 'top');",
                                                                           "}"
                                                                       ),
                                                                       autoWidth = FALSE, ## use smart column width handling
                                                                       server = TRUE,   ## use client-side processing
                                                                       processing = TRUE,
                                                                       dom = 'lfrtBip',
                                                                       buttons = c('csv', 'excel')
                                                        ),
                                                        extensions = 'Buttons',
                                                        selection = 'single', ## enable selection of a single row
                                                        filter = 'top',              ## include column filters at the top
                                                        rownames = FALSE                ## don't show row numbers/names
    )
    })
    ## blast tables
    output$blast_mesotaenium <- DT::renderDataTable({datatable(blast_mesotaenium,
                                                               options = list(paging = TRUE,    ## paginate the output
                                                                              pageLength = 25,  ## number of rows to output for each page
                                                                              lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                c("5", "10", "25", "50", "100", 'All')),
                                                                              scrollX = TRUE,   ## enable scrolling on X axis
                                                                              scrollY = "500px",   ## enable scrolling on Y axis
                                                                              rowCallback = JS(
                                                                                  "function(row, data, index) {",
                                                                                  "  $('td', row).css('vertical-align', 'top');",
                                                                                  "}"
                                                                              ),
                                                                              autoWidth = FALSE, ## use smart column width handling
                                                                              server = TRUE,   ## use client-side processing
                                                                              processing = TRUE,
                                                                              dom = 'lfrtBip',
                                                                              buttons = c('csv', 'excel')
                                                               ),
                                                               extensions = 'Buttons',
                                                               selection = 'single', ## enable selection of a single row
                                                               filter = 'top',              ## include column filters at the top
                                                               rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$blast_physcomitrium <- DT::renderDataTable({datatable(blast_physcomitrium,
                                                                 options = list(paging = TRUE,    ## paginate the output
                                                                                pageLength = 25,  ## number of rows to output for each page
                                                                                lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                                  c("5", "10", "25", "50", "100", 'All')),
                                                                                scrollX = TRUE,   ## enable scrolling on X axis
                                                                                scrollY = "500px",   ## enable scrolling on Y axis
                                                                                rowCallback = JS(
                                                                                    "function(row, data, index) {",
                                                                                    "  $('td', row).css('vertical-align', 'top');",
                                                                                    "}"
                                                                                ),
                                                                                autoWidth = FALSE, ## use smart column width handling
                                                                                server = TRUE,   ## use client-side processing
                                                                                processing = TRUE,
                                                                                dom = 'lfrtBip',
                                                                                buttons = c('csv', 'excel')
                                                                 ),
                                                                 extensions = 'Buttons',
                                                                 selection = 'single', ## enable selection of a single row
                                                                 filter = 'top',              ## include column filters at the top
                                                                 rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$blast_zygnema <- DT::renderDataTable({datatable(blast_zygnema,
                                                           options = list(paging = TRUE,    ## paginate the output
                                                                          pageLength = 25,  ## number of rows to output for each page
                                                                          lengthMenu = list(c(5, 10, 25, 50, 100, -1), 
                                                                                            c("5", "10", "25", "50", "100", 'All')),
                                                                          scrollX = TRUE,   ## enable scrolling on X axis
                                                                          scrollY = "500px",   ## enable scrolling on Y axis
                                                                          rowCallback = JS(
                                                                              "function(row, data, index) {",
                                                                              "  $('td', row).css('vertical-align', 'top');",
                                                                              "}"
                                                                          ),
                                                                          autoWidth = FALSE, ## use smart column width handling
                                                                          server = TRUE,   ## use client-side processing
                                                                          processing = TRUE,
                                                                          dom = 'lfrtBip',
                                                                          buttons = c('csv', 'excel')
                                                           ),
                                                           extensions = 'Buttons',
                                                           selection = 'single', ## enable selection of a single row
                                                           filter = 'top',              ## include column filters at the top
                                                           rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$orthofinder_n0 <- DT::renderDataTable({datatable(orthofinder_N0,
                                                            options = list(paging = TRUE,    ## paginate the output
                                                                           pageLength = 1,  ## number of rows to output for each page
                                                                           lengthMenu = list(c(1, 5, 10, 25, 50, 100, -1), 
                                                                                             c("1", "5", "10", "25", "50", "100", 'All')),
                                                                           scrollX = TRUE,   ## enable scrolling on X axis
                                                                           scrollY = "500px",   ## enable scrolling on Y axis
                                                                           rowCallback = JS(
                                                                               "function(row, data, index) {",
                                                                               "  $('td', row).css('vertical-align', 'top');",
                                                                               "}"
                                                                           ),
                                                                           autoWidth = FALSE, ## use smart column width handling
                                                                           server = TRUE,   ## use client-side processing
                                                                           processing = TRUE,
                                                                           dom = 'lfrtBip',
                                                                           buttons = c('csv', 'excel')
                                                            ),
                                                            extensions = 'Buttons',
                                                            selection = 'single', ## enable selection of a single row
                                                            filter = 'top',              ## include column filters at the top
                                                            rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$msa_viewer_output <- renderMsaR(
        msaR(ape::read.FASTA(paste0("./assets/orthofinder/MultipleSequenceAlignments/", input$msa_viewer, ".fa"), type="AA"),
             menu=FALSE, overviewbox = FALSE,  colorscheme = "clustal2",
             width = NULL, height = NULL, rowheight = 20, alignmentHeight = 500,
             seqlogo = TRUE, conservation = TRUE, markers = TRUE, metacell = FALSE,
             leftheader = TRUE, labels = TRUE, labelname = TRUE, labelid = FALSE,
             labelNameLength = 300, overviewboxWidth = "auto", overviewboxHeight = "fixed")
    )
    # Gene expression profiler
    output$mesotaenium_gene_visualization_plot <- renderPlot({
        if(input$gene_visualization_type == "Boxplot_conditions"){
            # Visualize gene expression for all conditions
            # Define a palette
            my_colors <- c("standard_growth_0" = "#e2e6e6",
                           "standard_growth_0.5" = "#c6cdce",
                           "standard_growth_1" = "#a8b4b5",
                           "standard_growth_2" = "#8b9b9d",
                           "standard_growth_4" = "#6e8284",
                           "standard_growth_6" = "#51696c",
                           "standard_growth_24" = "#345053",
                           "cold_0.5" = "#91e0ef",
                           "cold_1" = "#48cae4",
                           "cold_2" = "#00b4d8",
                           "cold_4" = "#0296c7",
                           "cold_6" = "#0077b6",
                           "cold_24" = "#023e8a",
                           "heat_0.5" = "#fdc3c3",
                           "heat_1" = "#feaeae",
                           "heat_2" = "#fe9a9a",
                           "heat_4" = "#fe8585",
                           "heat_6" = "#fe5d5d",
                           "heat_24" = "#fe4848",
                           "high_light_s_0.25" = "#ffea01",
                           "high_light_s_0.5" = "#ffdd02",
                           "high_light_s_2" = "#ffd002",
                           "high_light_s_4" = "#ffb700",
                           "high_light_s_6" = "#ffa201",
                           "high_light_r_0.25" = "#dee602",
                           "high_light_r_0.5" = "#c3d205",
                           "high_light_r_1" = "#a7be07",
                           "high_light_r_2" = "#8daa0b",
                           "high_light_r_4" = "#70960e"
            )
            p <- ggplot(mesotaenium_expr, aes(condition, scale_z_score(get(input$mesotaenium_gene_id_visualization)))) +
                gghalves::geom_half_point(side = "l", range_scale = .5, alpha = .5, color="black", fill= my_colors, shape=21, size=1 ) +
                geom_boxplot(width = .3, outlier.shape = NA, fill= my_colors,  alpha = 0.3) +
                ggtitle(paste0("Box plot of Z-score transformed log2(CPM) normalized and voom transformed of ", input$mesotaenium_gene_id_visualization," - M. endlicherianum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Condition") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        } else if(input$gene_visualization_type == "Boxplot_treatment"){
            # Visualize gene expression for all treatments
            my_colors <- c("standard_growth" = "#345053",
                           "cold" = "#023e8a",
                           "heat" = "#fe4848",
                           "high_light_s" = "#ffa201",
                           "high_light_r" = "#70960e"
            )
            p <- ggplot(mesotaenium_expr, aes(treatment, scale_z_score(get(input$mesotaenium_gene_id_visualization)))) +
                ggdist::stat_halfeye(adjust = .5, width = .5, .width = 0, justification = -.4) +
                gghalves::geom_half_point(side = "l", range_scale = .5, alpha = .5, color="black", fill= my_colors, shape=21, size=1 ) +
                geom_boxplot(width = .3, outlier.shape = NA, fill= my_colors,  alpha = 0.3) +
                ggtitle(paste0("Box plot of Z-score transformed Log2(CPM) normalized and voom transformed of ", input$mesotaenium_gene_id_visualization," - M. endlicherianum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Treatment") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        } else if(input$gene_visualization_type == "Dotplot"){
            my_colors <- c("standard_growth" = "#345053",
                           "cold" = "#023e8a",
                           "heat" = "#fe4848",
                           "high_light_s" = "#ffa201",
                           "high_light_r" = "#70960e"
            )
            p <- ggplot(mesotaenium_expr, aes(x = time, y = scale_z_score(get(input$mesotaenium_gene_id_visualization)), color=treatment)) +
                geom_point(aes(fill=treatment), shape=21, size=2, alpha = .3) +
                stat_smooth(method = loess, formula = y ~ x,  linewidth = 0.5, alpha = .2, se = FALSE, show.legend = FALSE) +
                ggtitle(paste0("Dot plot of Z-score transformed Log2(CPM) normalized and voom transformed of ", input$mesotaenium_gene_id_visualization," - M. endlicherianum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Time (hour)") +
                scale_fill_manual(values = my_colors) +
                scale_color_manual(values = my_colors) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      ) +
                scale_x_continuous(breaks = unique(mesotaenium_expr$time))
        } else if (input$gene_visualization_type == "Heatmap"){
            p <- ggplot(mesotaenium_expr,  aes(x = sample_name, y = 1, fill = scale_z_score(get(input$mesotaenium_gene_id_visualization)))) +
                geom_tile(color = "white",
                          lwd = 0.5,
                          linetype = 1) +
                scale_fill_gradient2(low = "#033270",
                                     mid = "#FFFFFF",
                                     high = "#cb1b16", name = "Z-Score") +
                ggtitle(paste0("Heatmap of z-score transformed of Log2(CPM) normalized and voom transformed of ", input$mesotaenium_gene_id_visualization," - M. endlicherianum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Condition") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        }
        p
    }, width = 1000, height = 800)
    output$zygnema_gene_visualization_plot <- renderPlot({
        if(input$gene_visualization_type == "Boxplot_conditions"){
            # Visualize gene expression for all conditions
            # Define a palette
            my_colors <- c("standard_growth_0" = "#e2e6e6",
                           "standard_growth_0.5" = "#c6cdce",
                           "standard_growth_1" = "#a8b4b5",
                           "standard_growth_2" = "#8b9b9d",
                           "standard_growth_4" = "#6e8284",
                           "standard_growth_6" = "#51696c",
                           "standard_growth_24" = "#345053",
                           "cold_0.5" = "#91e0ef",
                           "cold_1" = "#48cae4",
                           "cold_2" = "#00b4d8",
                           "cold_4" = "#0296c7",
                           "cold_6" = "#0077b6",
                           "cold_24" = "#023e8a",
                           "heat_0.5" = "#fdc3c3",
                           "heat_1" = "#feaeae",
                           "heat_2" = "#fe9a9a",
                           "heat_4" = "#fe8585",
                           "heat_6" = "#fe5d5d",
                           "heat_24" = "#fe4848",
                           "high_light_s_0.25" = "#ffea01",
                           "high_light_s_0.5" = "#ffdd02",
                           "high_light_s_2" = "#ffd002",
                           "high_light_s_4" = "#ffb700",
                           "high_light_s_6" = "#ffa201",
                           "high_light_r_0.25" = "#dee602",
                           "high_light_r_0.5" = "#c3d205",
                           "high_light_r_1" = "#a7be07",
                           "high_light_r_2" = "#8daa0b",
                           "high_light_r_4" = "#70960e"
            )
            p <- ggplot(zygnema_expr, aes(condition, scale_z_score(get(input$zygnema_gene_id_visualization)))) +
                gghalves::geom_half_point(side = "l", range_scale = .5, alpha = .5, color="black", fill= my_colors, shape=21, size=1 ) +
                geom_boxplot(width = .3, outlier.shape = NA, fill= my_colors,  alpha = 0.3) +
                ggtitle(paste0("Box plot of Z-score transformed log2(CPM) normalized and voom transformed of ", input$zygnema_gene_id_visualization," - Z. circumcarinatum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Condition") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        } else if(input$gene_visualization_type == "Boxplot_treatment"){
            # Visualize gene expression for all treatments
            my_colors <- c("standard_growth" = "#345053",
                           "cold" = "#023e8a",
                           "heat" = "#fe4848",
                           "high_light_s" = "#ffa201",
                           "high_light_r" = "#70960e"
            )
            p <- ggplot(zygnema_expr, aes(treatment, scale_z_score(get(input$zygnema_gene_id_visualization)))) +
                ggdist::stat_halfeye(adjust = .5, width = .5, .width = 0, justification = -.4) +
                gghalves::geom_half_point(side = "l", range_scale = .5, alpha = .5, color="black", fill= my_colors, shape=21, size=1 ) +
                geom_boxplot(width = .3, outlier.shape = NA, fill= my_colors,  alpha = 0.3) +
                ggtitle(paste0("Box plot of Z-score transformed Log2(CPM) normalized and voom transformed of ", input$zygnema_gene_id_visualization," - Z. circumcarinatum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Treatment") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        } else if(input$gene_visualization_type == "Dotplot"){
            my_colors <- c("standard_growth" = "#345053",
                           "cold" = "#023e8a",
                           "heat" = "#fe4848",
                           "high_light_s" = "#ffa201",
                           "high_light_r" = "#70960e"
            )
            p <- ggplot(zygnema_expr, aes(x = time, y = scale_z_score(get(input$zygnema_gene_id_visualization)), color=treatment)) +
                geom_point(aes(fill=treatment), shape=21, size=2, alpha = .3) +
                stat_smooth(method = loess, formula = y ~ x,  linewidth = 0.5, alpha = .2, se = FALSE, show.legend = FALSE) +
                ggtitle(paste0("Dot plot of Z-score transformed Log2(CPM) normalized and voom transformed of ", input$zygnema_gene_id_visualization," - Z. circumcarinatum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Time (hour)") +
                scale_fill_manual(values = my_colors) +
                scale_color_manual(values = my_colors) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      ) +
                scale_x_continuous(breaks = unique(zygnema_expr$time))
        } else if (input$gene_visualization_type == "Heatmap"){
            p <- ggplot(zygnema_expr,  aes(x = sample_name, y = 1, fill = scale_z_score(get(input$zygnema_gene_id_visualization)))) +
                geom_tile(color = "white",
                          lwd = 0.5,
                          linetype = 1) +
                scale_fill_gradient2(low = "#033270",
                                     mid = "#FFFFFF",
                                     high = "#cb1b16", name = "Z-Score") +
                ggtitle(paste0("Heatmap of z-score transformed of Log2(CPM) normalized and voom transformed of ", input$zygnema_gene_id_visualization," - Z. circumcarinatum")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Condition") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        }
        p
    }, width = 1000, height = 800)
    output$physcomitrium_gene_visualization_plot <- renderPlot({
        if(input$gene_visualization_type == "Boxplot_conditions"){
            # Visualize gene expression for all conditions
            # Define a palette
            my_colors <- c("standard_growth_0" = "#e2e6e6",
                           "standard_growth_0.5" = "#c6cdce",
                           "standard_growth_1" = "#a8b4b5",
                           "standard_growth_2" = "#8b9b9d",
                           "standard_growth_4" = "#6e8284",
                           "standard_growth_6" = "#51696c",
                           "standard_growth_24" = "#345053",
                           "cold_0.5" = "#91e0ef",
                           "cold_1" = "#48cae4",
                           "cold_2" = "#00b4d8",
                           "cold_4" = "#0296c7",
                           "cold_6" = "#0077b6",
                           "cold_24" = "#023e8a",
                           "heat_0.5" = "#fdc3c3",
                           "heat_1" = "#feaeae",
                           "heat_2" = "#fe9a9a",
                           "heat_4" = "#fe8585",
                           "heat_6" = "#fe5d5d",
                           "heat_24" = "#fe4848",
                           "high_light_s_0.25" = "#ffea01",
                           "high_light_s_0.5" = "#ffdd02",
                           "high_light_s_2" = "#ffd002",
                           "high_light_s_4" = "#ffb700",
                           "high_light_s_6" = "#ffa201",
                           "high_light_r_0.25" = "#dee602",
                           "high_light_r_0.5" = "#c3d205",
                           "high_light_r_1" = "#a7be07",
                           "high_light_r_2" = "#8daa0b",
                           "high_light_r_4" = "#70960e"
            )
            p <- ggplot(physcomitrium_expr, aes(condition, scale_z_score(get(input$physcomitrium_gene_id_visualization)))) +
                gghalves::geom_half_point(side = "l", range_scale = .5, alpha = .5, color="black", fill= my_colors, shape=21, size=1 ) +
                geom_boxplot(width = .3, outlier.shape = NA, fill= my_colors,  alpha = 0.3) +
                ggtitle(paste0("Box plot of Z-score transformed log2(CPM) normalized and voom transformed of ", input$physcomitrium_gene_id_visualization," - P. patens")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Condition") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        } else if(input$gene_visualization_type == "Boxplot_treatment"){
            # Visualize gene expression for all treatments
            my_colors <- c("standard_growth" = "#345053",
                           "cold" = "#023e8a",
                           "heat" = "#fe4848",
                           "high_light_s" = "#ffa201",
                           "high_light_r" = "#70960e"
            )
            p <- ggplot(physcomitrium_expr, aes(treatment, scale_z_score(get(input$physcomitrium_gene_id_visualization)))) +
                ggdist::stat_halfeye(adjust = .5, width = .5, .width = 0, justification = -.4) +
                gghalves::geom_half_point(side = "l", range_scale = .5, alpha = .5, color="black", fill= my_colors, shape=21, size=1 ) +
                geom_boxplot(width = .3, outlier.shape = NA, fill= my_colors,  alpha = 0.3) +
                ggtitle(paste0("Box plot of Z-score transformed Log2(CPM) normalized and voom transformed of ", input$physcomitrium_gene_id_visualization," - P. patens")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Treatment") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        } else if(input$gene_visualization_type == "Dotplot"){
            my_colors <- c("standard_growth" = "#345053",
                           "cold" = "#023e8a",
                           "heat" = "#fe4848",
                           "high_light_s" = "#ffa201",
                           "high_light_r" = "#70960e"
            )
            p <- ggplot(physcomitrium_expr, aes(x = time, y = scale_z_score(get(input$physcomitrium_gene_id_visualization)), color=treatment)) +
                geom_point(aes(fill=treatment), shape=21, size=2, alpha = .3) +
                stat_smooth(method = loess, formula = y ~ x,  linewidth = 0.5, alpha = .2, se = FALSE, show.legend = FALSE) +
                ggtitle(paste0("Dot plot of Z-score transformed Log2(CPM) normalized and voom transformed of ", input$physcomitrium_gene_id_visualization," - P. patens")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Time (hour)") +
                scale_fill_manual(values = my_colors) +
                scale_color_manual(values = my_colors) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      ) +
                scale_x_continuous(breaks = unique(physcomitrium_expr$time))
        } else if (input$gene_visualization_type == "Heatmap"){
            p <- ggplot(physcomitrium_expr,  aes(x = sample_name, y = 1, fill = scale_z_score(get(input$physcomitrium_gene_id_visualization)))) +
                geom_tile(color = "white",
                          lwd = 0.5,
                          linetype = 1) +
                scale_fill_gradient2(low = "#033270",
                                     mid = "#FFFFFF",
                                     high = "#cb1b16", name = "Z-Score") +
                ggtitle(paste0("Heatmap of z-score transformed of Log2(CPM) normalized and voom transformed of ", input$physcomitrium_gene_id_visualization," - P. patens")) +
                labs(y = "Z-score transformed normalized Log2(CPM)", x = "Condition") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                      )
        }
        p
    }, width = 1000, height = 800)
    # Co-expression networks
    # WGCNA
    output$mesotaenium_wgcna_table <- DT::renderDataTable({
        read_csv(normalizePath(file.path("./assets/WGCNA/tables/mesotaenium/WGCNA_Results.csv"))) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$zygnema_wgcna_table <- DT::renderDataTable({
        read_csv(normalizePath(file.path("./assets/WGCNA/tables/zygnema/WGCNA_Results.csv"))) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$physcomitrium_wgcna_table <- DT::renderDataTable({
        read_csv(normalizePath(file.path("./assets/WGCNA/tables/physcomitrium/WGCNA_Results.csv"))) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$mesotaenium_wgcna_top_hubs <- DT::renderDataTable({
        hubs_mesotaenium <- read_csv(normalizePath(file.path("./assets/WGCNA/tables/mesotaenium/hubs/", paste0(input$cluster_wgcna_mesotaenium, ".csv"))))
        colnames(hubs_mesotaenium) <- c("Rank", "Gene ID")
        hubs_mesotaenium %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 20,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 20, -1),
                                                       c("5", "10", "20", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$zygnema_wgcna_top_hubs <- DT::renderDataTable({
        hubs_zygnema <- read_csv(normalizePath(file.path("./assets/WGCNA/tables/zygnema/hubs/", paste0(input$cluster_wgcna_zygnema, ".csv")))) 
        colnames(hubs_zygnema) <- c("Rank", "Gene ID")
        hubs_zygnema %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 20,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 20, -1),
                                                       c("5", "10", "20", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$physcomitrium_wgcna_top_hubs <- DT::renderDataTable({
        hubs_physcomitrium <- read_csv(normalizePath(file.path("./assets/WGCNA/tables/physcomitrium/hubs/", paste0(input$cluster_wgcna_physcomitrium, ".csv")))) 
        colnames(hubs_physcomitrium) <- c("Rank", "Gene ID")
        hubs_physcomitrium%>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 20,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 20, -1),
                                                       c("5", "10", "20", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$mesotaenium_wgcna_ora_table <- DT::renderDataTable({
        cluster_genes_mesotaenium_wgcna_table <- read_csv(normalizePath(file.path("./assets/wgcna/tables/mesotaenium/WGCNA_Results.csv"))) %>%
            filter(moduleColors == input$cluster_wgcna_mesotaenium) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        universe_mesotaenium_wgcna_table <- read_csv(normalizePath(file.path("./assets/wgcna/tables/mesotaenium/WGCNA_Results.csv"))) %>%
            left_join(GO_mesotaenium, by = join_by("GeneID" == "gene")) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        enrichment_mesotaenium_wgcna_table <- enricher(
            cluster_genes_mesotaenium_wgcna_table,
            TERM2GENE = GO_mesotaenium %>%
                filter(domain == input$go_domain_wgcna_mesotaenium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_mesotaenium %>%
                filter(domain == input$go_domain_wgcna_mesotaenium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_wgcna_mesotaenium,
            universe = universe_mesotaenium_wgcna_table,
            qvalueCutoff = input$go_qvalueCutoff_wgcna_mesotaenium
        )
        if (!is.null(enrichment_mesotaenium_wgcna_table@result)){
            enrichment_mesotaenium_wgcna_table@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$zygnema_wgcna_ora_table <- DT::renderDataTable({
        cluster_genes_zygnema_wgcna_table <- read_csv(normalizePath(file.path("./assets/wgcna/tables/zygnema/WGCNA_Results.csv"))) %>%
            filter(moduleColors == input$cluster_wgcna_zygnema) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        universe_zygnema_wgcna_table <- read_csv(normalizePath(file.path("./assets/wgcna/tables/zygnema/WGCNA_Results.csv"))) %>%
            left_join(GO_zygnema, by = join_by("GeneID" == "gene")) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        enrichment_zygnema_wgcna_table <- enricher(
            cluster_genes_zygnema_wgcna_table,
            TERM2GENE = GO_zygnema %>%
                filter(domain == input$go_domain_wgcna_zygnema) %>%
                select(c(term, gene)),
            TERM2NAME = GO_zygnema %>%
                filter(domain == input$go_domain_wgcna_zygnema) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_wgcna_zygnema,
            universe = universe_zygnema_wgcna_table,
            qvalueCutoff = input$go_qvalueCutoff_wgcna_zygnema
        )
        if (!is.null(enrichment_zygnema_wgcna_table@result)){
            enrichment_zygnema_wgcna_table@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$physcomitrium_wgcna_ora_table <- DT::renderDataTable({
        cluster_genes_physcomitrium_wgcna_table <- read_csv(normalizePath(file.path("./assets/wgcna/tables/physcomitrium/WGCNA_Results.csv"))) %>%
            filter(moduleColors == input$cluster_wgcna_physcomitrium) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        universe_physcomitrium_wgcna_table <- read_csv(normalizePath(file.path("./assets/wgcna/tables/physcomitrium/WGCNA_Results.csv"))) %>%
            left_join(GO_physcomitrium, by = join_by("GeneID" == "gene")) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        enrichment_physcomitrium_wgcna_table <- enricher(
            cluster_genes_physcomitrium_wgcna_table,
            TERM2GENE = GO_physcomitrium %>%
                filter(domain == input$go_domain_wgcna_physcomitrium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_physcomitrium %>%
                filter(domain == input$go_domain_wgcna_physcomitrium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_wgcna_physcomitrium,
            universe = universe_physcomitrium_wgcna_table,
            qvalueCutoff = input$go_qvalueCutoff_wgcna_physcomitrium
        )
        if (!is.null(enrichment_physcomitrium_wgcna_table@result)){
            enrichment_physcomitrium_wgcna_table@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$mesotaenium_wgcna_ora_visualization <- renderPlot({
        cluster_genes_mesotaenium_wgcna <- read_csv(normalizePath(file.path("./assets/wgcna/tables/mesotaenium/WGCNA_Results.csv"))) %>%
            filter(moduleColors == input$cluster_wgcna_mesotaenium) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        universe_mesotaenium_wgcna <- read_csv(normalizePath(file.path("./assets/wgcna/tables/mesotaenium/WGCNA_Results.csv"))) %>%
            left_join(GO_mesotaenium, by = join_by("GeneID" == "gene")) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        enrichment_mesotaenium_wgcna <- enricher(
            cluster_genes_mesotaenium_wgcna,
            TERM2GENE = GO_mesotaenium %>%
                filter(domain == input$go_domain_wgcna_mesotaenium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_mesotaenium %>%
                filter(domain == input$go_domain_wgcna_mesotaenium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_wgcna_mesotaenium,
            universe = universe_mesotaenium_wgcna,
            qvalueCutoff = input$go_qvalueCutoff_wgcna_mesotaenium
        )
        if (!is.null(enrichment_mesotaenium_wgcna@result)){
            p <- dotplot(enrichment_mesotaenium_wgcna,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category_wgcna_mesotaenium,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 600, height = 900)
    output$zygnema_wgcna_ora_visualization <- renderPlot({
        cluster_genes_zygnema_wgcna <- read_csv(normalizePath(file.path("./assets/wgcna/tables/zygnema/WGCNA_Results.csv"))) %>%
            filter(moduleColors == input$cluster_wgcna_zygnema) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        universe_zygnema_wgcna <- read_csv(normalizePath(file.path("./assets/wgcna/tables/zygnema/WGCNA_Results.csv"))) %>%
            left_join(GO_zygnema, by = join_by("GeneID" == "gene")) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        enrichment_zygnema_wgcna <- enricher(
            cluster_genes_zygnema_wgcna,
            TERM2GENE = GO_zygnema %>%
                filter(domain == input$go_domain_wgcna_zygnema) %>%
                select(c(term, gene)),
            TERM2NAME = GO_zygnema %>%
                filter(domain == input$go_domain_wgcna_zygnema) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_wgcna_zygnema,
            universe = universe_zygnema_wgcna,
            qvalueCutoff = input$go_qvalueCutoff_wgcna_zygnema
        )
        if (!is.null(enrichment_zygnema_wgcna@result)){
            p <- dotplot(enrichment_zygnema_wgcna,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category_wgcna_zygnema,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 600, height = 900)
    output$physcomitrium_wgcna_ora_visualization <- renderPlot({
        cluster_genes_physcomitrium_wgcna <- read_csv(normalizePath(file.path("./assets/wgcna/tables/physcomitrium/WGCNA_Results.csv"))) %>%
            filter(moduleColors == input$cluster_wgcna_physcomitrium) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        universe_physcomitrium_wgcna <- read_csv(normalizePath(file.path("./assets/wgcna/tables/physcomitrium/WGCNA_Results.csv"))) %>%
            left_join(GO_physcomitrium, by = join_by("GeneID" == "gene")) %>%
            select(GeneID) %>%
            unlist() %>%
            as.vector()
        enrichment_physcomitrium_wgcna <- enricher(
            cluster_genes_physcomitrium_wgcna,
            TERM2GENE = GO_physcomitrium %>%
                filter(domain == input$go_domain_wgcna_physcomitrium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_physcomitrium %>%
                filter(domain == input$go_domain_wgcna_physcomitrium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_wgcna_physcomitrium,
            universe = universe_physcomitrium_wgcna,
            qvalueCutoff = input$go_qvalueCutoff_wgcna_physcomitrium
        )
        if (!is.null(enrichment_physcomitrium_wgcna@result)){
            p <- dotplot(enrichment_physcomitrium_wgcna,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category_wgcna_physcomitrium,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 600, height = 900)
    output$mesotaenium_wgcna_expression_visualization <- renderImage({
        filename <- normalizePath(file.path("./assets/WGCNA/plots/mesotaenium/expression_profiles", paste0(input$cluster_wgcna_mesotaenium ,".pdf.png")))
        list(src = filename,
             alt = paste("Image number"),
             width = "1142px",
             height = "826"
        )
    }, deleteFile = FALSE)
    output$zygnema_wgcna_expression_visualization <- renderImage({
        filename <- normalizePath(file.path("./assets/WGCNA/plots/zygnema/expression_profiles", paste0(input$cluster_wgcna_zygnema ,".pdf.png")))
        list(src = filename,
             alt = paste("Image number"),
             width = "1142px",
             height = "826"
        )
    }, deleteFile = FALSE)
    output$physcomitrium_wgcna_expression_visualization <- renderImage({
        filename <- normalizePath(file.path("./assets/WGCNA/plots/physcomitrium/expression_profiles", paste0(input$cluster_wgcna_physcomitrium ,".pdf.png")))
        list(src = filename,
             alt = paste("Image number"),
             width = "1142px",
             height = "826"
        )
    }, deleteFile = FALSE)
    output$mesotaenium_wgcna_GS_visualization <- renderImage({
        if (input$wgcna_gene_significance_type == "absolute value"){
            filename <- normalizePath(file.path("./assets/WGCNA/plots/mesotaenium/Gene_Significance/absolute_value", paste0(input$wgcna_gene_significance_measure ,".pdf.png")))
            list(src = filename,
                 alt = paste("Image number"),
                 width = "1800px",
                 height = "1000"
            )
        } else if (input$wgcna_gene_significance_type == "value"){
            filename <- normalizePath(file.path("./assets/WGCNA/plots/mesotaenium/Gene_Significance/value", paste0(input$wgcna_gene_significance_measure ,".pdf.png")))
            list(src = filename,
                 alt = paste("Image number"),
                 width = "1800px",
                 height = "1000"
            )
        }
    }, deleteFile = FALSE)
    output$zygnema_wgcna_GS_visualization <- renderImage({
        if (input$wgcna_gene_significance_type == "absolute value"){
            filename <- normalizePath(file.path("./assets/WGCNA/plots/zygnema/Gene_Significance/absolute_value", paste0(input$wgcna_gene_significance_measure ,".pdf.png")))
            list(src = filename,
                 alt = paste("Image number"),
                 width = "1800px",
                 height = "1000"
            )
        } else if (input$wgcna_gene_significance_type == "value"){
            filename <- normalizePath(file.path("./assets/WGCNA/plots/zygnema/Gene_Significance/value", paste0(input$wgcna_gene_significance_measure ,".pdf.png")))
            list(src = filename,
                 alt = paste("Image number"),
                 width = "1800px",
                 height = "1000"
            )
        }
    }, deleteFile = FALSE)
    output$physcomitrium_wgcna_GS_visualization <- renderImage({
        if (input$wgcna_gene_significance_type == "absolute value"){
            filename <- normalizePath(file.path("./assets/WGCNA/plots/physcomitrium/Gene_Significance/absolute_value", paste0(input$wgcna_gene_significance_measure ,".pdf.png")))
            list(src = filename,
                 alt = paste("Image number"),
                 width = "1800px",
                 height = "1000"
            )
        } else if (input$wgcna_gene_significance_type == "value"){
            filename <- normalizePath(file.path("./assets/WGCNA/plots/physcomitrium/Gene_Significance/value", paste0(input$wgcna_gene_significance_measure ,".pdf.png")))
            list(src = filename,
                 alt = paste("Image number"),
                 width = "1800px",
                 height = "1000"
            )
        }
    }, deleteFile = FALSE)
    # DPGP
    output$dpgp_mesotaenium_dynamicInput <- renderUI({
        # Render different input based on the selection
        if (input$dpgp_input_mesotaenium == "mesotaenium_cold_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_mesotaenium", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("cold 1", "cold 2", "cold 3", "cold 4", "cold 5", "cold 6", "cold 7", "cold 8", "cold 9", "cold 10", "cold 11", "cold 12", "cold 13", "cold 14", "cold 15", "cold 16"), selected = "cold 1")
        } else if (input$dpgp_input_mesotaenium == "mesotaenium_heat_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_mesotaenium", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("heat 1", "heat 2", "heat 3", "heat 4", "heat 5", "heat 6", "heat 7", "heat 8", "heat 9", "heat 10", "heat 11", "heat 12"), selected = "heat 1")
        } else if(input$dpgp_input_mesotaenium == "mesotaenium_high_light_no_recovery_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_mesotaenium", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("highlight 1", "highlight 2", "highlight 3", "highlight 4", "highlight 5", "highlight 6", "highlight 7", "highlight 8", "highlight 9", "highlight 10", "highlight 11", "highlight 12", "highlight 13"), selected = "highlight 6")
        }
    })
    output$dpgp_zygnema_dynamicInput <- renderUI({
        # Render different input based on the selection
        if (input$dpgp_input_zygnema == "zygnema_cold_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_zygnema", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("cold 1", "cold 2", "cold 3", "cold 4", "cold 5", "cold 6", "cold 7", "cold 8", "cold 9", "cold 10", "cold 11", "cold 12", "cold 13", "cold 14", "cold 15", "cold 16"), selected = "cold 1")
        } else if (input$dpgp_input_zygnema == "zygnema_heat_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_zygnema", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("heat 1", "heat 2", "heat 3", "heat 4", "heat 5", "heat 6", "heat 7", "heat 8", "heat 9", "heat 10", "heat 11", "heat 12"), selected = "heat 1")
        } else if(input$dpgp_input_zygnema == "zygnema_high_light_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_zygnema", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("highlight 1", "highlight 2", "highlight 3", "highlight 4", "highlight 5", "highlight 6", "highlight 7", "highlight 8", "highlight 9", "highlight 10", "highlight 11", "highlight 12", "highlight 13"), selected = "highlight 6")
        }
    })
    output$dpgp_physcomitrium_dynamicInput <- renderUI({
        # Render different input based on the selection
        if (input$dpgp_input_physcomitrium == "physcomitrium_cold_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_physcomitrium", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("cold 1", "cold 2", "cold 3", "cold 4", "cold 5", "cold 6", "cold 7", "cold 8", "cold 9", "cold 10", "cold 11", "cold 12", "cold 13", "cold 14", "cold 15", "cold 16"), selected = "cold 1")
        } else if (input$dpgp_input_physcomitrium == "physcomitrium_heat_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_physcomitrium", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("heat 1", "heat 2", "heat 3", "heat 4", "heat 5", "heat 6", "heat 7", "heat 8", "heat 9", "heat 10", "heat 11", "heat 12"), selected = "heat 1")
        } else if(input$dpgp_input_physcomitrium == "physcomitrium_high_light_optimal_clustering") {
            selectInput(inputId = "cluster_dpgp_physcomitrium", label = h4("Select a cluster of this treatment for ORA analysis and visualization"), choices = c("highlight 1", "highlight 2", "highlight 3", "highlight 4", "highlight 5", "highlight 6", "highlight 7", "highlight 8", "highlight 9", "highlight 10", "highlight 11", "highlight 12", "highlight 13"), selected = "highlight 6")
        }
    })
    output$mesotaenium_dpgp_table <- DT::renderDataTable({
        read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_mesotaenium) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                         ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$zygnema_dpgp_table <- DT::renderDataTable({
        read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_zygnema) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$physcomitrium_dpgp_table <- DT::renderDataTable({
        read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_physcomitrium) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$mesotaenium_dpgp_ora_table <- DT::renderDataTable({
        cluster_genes_mesotaenium_table <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_mesotaenium) %>%
            filter(cluster == str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][2]) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        universe_mesotaenium_table <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
            left_join(GO_mesotaenium, by = join_by("gene" == "gene")) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        enrichment_mesotaenium_table <- enricher(
            cluster_genes_mesotaenium_table,
            TERM2GENE = GO_mesotaenium %>%
                filter(domain == input$go_domain_dpgp_mesotaenium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_mesotaenium %>%
                filter(domain == input$go_domain_dpgp_mesotaenium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_dpgp_mesotaenium,
            universe = universe_mesotaenium_table,
            qvalueCutoff = input$go_qvalueCutoff_dpgp_mesotaenium
        )
        if (!is.null(enrichment_mesotaenium_table@result)){
            enrichment_mesotaenium_table@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$zygnema_dpgp_ora_table <- DT::renderDataTable({
        cluster_genes_zygnema_table <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_zygnema) %>%
            filter(cluster == str_split(input$cluster_dpgp_zygnema, " ")[[1]][2]) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        universe_zygnema_table <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
            left_join(GO_zygnema, by = join_by("gene" == "gene")) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        enrichment_zygnema_table <- enricher(
            cluster_genes_zygnema_table,
            TERM2GENE = GO_zygnema %>%
                filter(domain == input$go_domain_dpgp_zygnema) %>%
                select(c(term, gene)),
            TERM2NAME = GO_zygnema %>%
                filter(domain == input$go_domain_dpgp_zygnema) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_dpgp_zygnema,
            universe = universe_zygnema_table,
            qvalueCutoff = input$go_qvalueCutoff_dpgp_zygnema
        )
        if (!is.null(enrichment_zygnema_table@result)){
            enrichment_zygnema_table@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$physcomitrium_dpgp_ora_table <- DT::renderDataTable({
        cluster_genes_physcomitrium_table <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_physcomitrium) %>%
            filter(cluster == str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][2]) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        universe_physcomitrium_table <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
            left_join(GO_physcomitrium, by = join_by("gene" == "gene")) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        enrichment_physcomitrium_table <- enricher(
            cluster_genes_physcomitrium_table,
            TERM2GENE = GO_physcomitrium %>%
                filter(domain == input$go_domain_dpgp_physcomitrium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_physcomitrium %>%
                filter(domain == input$go_domain_dpgp_physcomitrium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_dpgp_physcomitrium,
            universe = universe_physcomitrium_table,
            qvalueCutoff = input$go_qvalueCutoff_dpgp_physcomitrium
        )
        if (!is.null(enrichment_physcomitrium_table@result)){
            enrichment_physcomitrium_table@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$mesotaenium_dpgp_ora_visualization <- renderPlot({
        cluster_genes_mesotaenium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_mesotaenium) %>%
            filter(cluster == str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][2]) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        universe_mesotaenium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
            left_join(GO_mesotaenium, by = join_by("gene" == "gene")) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        enrichment_mesotaenium <- enricher(
            cluster_genes_mesotaenium,
            TERM2GENE = GO_mesotaenium %>%
                filter(domain == input$go_domain_dpgp_mesotaenium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_mesotaenium %>%
                filter(domain == input$go_domain_dpgp_mesotaenium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_dpgp_mesotaenium,
            universe = universe_mesotaenium,
            qvalueCutoff = input$go_qvalueCutoff_dpgp_mesotaenium
        )
        if (!is.null(enrichment_mesotaenium@result)){
            p <- dotplot(enrichment_mesotaenium,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category_dpgp_mesotaenium,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 600, height = 900)
    output$physcomitrium_dpgp_ora_visualization <- renderPlot({
        cluster_genes_physcomitrium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_physcomitrium) %>%
            filter(cluster == str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][2]) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        universe_physcomitrium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
            left_join(GO_physcomitrium, by = join_by("gene" == "gene")) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        enrichment_physcomitrium <- enricher(
            cluster_genes_physcomitrium,
            TERM2GENE = GO_physcomitrium %>%
                filter(domain == input$go_domain_dpgp_physcomitrium) %>%
                select(c(term, gene)),
            TERM2NAME = GO_physcomitrium %>%
                filter(domain == input$go_domain_dpgp_physcomitrium) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_dpgp_physcomitrium,
            universe = universe_physcomitrium,
            qvalueCutoff = input$go_qvalueCutoff_dpgp_physcomitrium
        )
        if (!is.null(enrichment_physcomitrium@result)){
            p <- dotplot(enrichment_physcomitrium,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category_dpgp_physcomitrium,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 600, height = 900)
    output$zygnema_dpgp_ora_visualization <- renderPlot({
        cluster_genes_zygnema <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
            filter(probability >= input$dpgp_probability_filter_threshold_zygnema) %>%
            filter(cluster == str_split(input$cluster_dpgp_zygnema, " ")[[1]][2]) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        universe_zygnema <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
            left_join(GO_zygnema, by = join_by("gene" == "gene")) %>%
            select(gene) %>%
            unlist() %>%
            as.vector()
        enrichment_zygnema <- enricher(
            cluster_genes_zygnema,
            TERM2GENE = GO_zygnema %>%
                filter(domain == input$go_domain_dpgp_zygnema) %>%
                select(c(term, gene)),
            TERM2NAME = GO_zygnema %>%
                filter(domain == input$go_domain_dpgp_zygnema) %>%
                select(c(term, name)),
            pvalueCutoff = input$go_pvalueCutoff_dpgp_zygnema,
            universe = universe_zygnema,
            qvalueCutoff = input$go_qvalueCutoff_dpgp_zygnema
        )
        if (!is.null(enrichment_zygnema@result)){
            p <- dotplot(enrichment_zygnema,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category_dpgp_zygnema,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 600, height = 900)
    output$mesotaenium_dpgp_expression_visualization <- renderPlot({
        if (str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][1] == "cold"){
            dpgp_cluster_mesotaenium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_mesotaenium) %>%
                filter(cluster == str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][2])
            cold_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/mesotaenium_cold.tsv"))) %>%
                right_join(dpgp_cluster_mesotaenium, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(cold_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#91e0ef") +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(0.5, 1)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(1, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(6, 24)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.5, 1, 2, 4, 6, 24)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        } else if (str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][1] == "heat"){
            dpgp_cluster_mesotaenium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_mesotaenium) %>%
                filter(cluster == str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][2])
            heat_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/mesotaenium_heat.tsv"))) %>%
                right_join(dpgp_cluster_mesotaenium, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(heat_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#fdc3c3") +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(0.5, 1)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(1, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(6, 24)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.5, 1, 2, 4, 6, 24)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        } else if (str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][1] == "highlight"){
            dpgp_cluster_mesotaenium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_mesotaenium, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_mesotaenium) %>%
                filter(cluster == str_split(input$cluster_dpgp_mesotaenium, " ")[[1]][2])
            highlight_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/mesoteanium_high_light.tsv"))) %>%
                right_join(dpgp_cluster_mesotaenium, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(highlight_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#ffea01") +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(0.25, 0.5)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(0.5, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c( 0.25, 0.5, 2, 4, 6)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        }
    }, width = 900, height = 900)
    output$zygnema_dpgp_expression_visualization <- renderPlot({
        if (str_split(input$cluster_dpgp_zygnema, " ")[[1]][1] == "cold"){
            dpgp_cluster_zygnema <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_zygnema) %>%
                filter(cluster == str_split(input$cluster_dpgp_zygnema, " ")[[1]][2])
            cold_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/zygnema_cold.tsv"))) %>%
                right_join(dpgp_cluster_zygnema, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(cold_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#91e0ef") +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(0.5, 1)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(1, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(6, 24)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.5, 1, 2, 4, 6, 24)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        } else if (str_split(input$cluster_dpgp_zygnema, " ")[[1]][1] == "heat"){
            dpgp_cluster_zygnema <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_zygnema) %>%
                filter(cluster == str_split(input$cluster_dpgp_zygnema, " ")[[1]][2])
            heat_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/zygnema_heat.tsv"))) %>%
                right_join(dpgp_cluster_zygnema, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(heat_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#fdc3c3") +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(0.5, 1)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(1, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(6, 24)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.5, 1, 2, 4, 6, 24)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        } else if (str_split(input$cluster_dpgp_zygnema, " ")[[1]][1] == "highlight"){
            dpgp_cluster_zygnema <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_zygnema, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_zygnema) %>%
                filter(cluster == str_split(input$cluster_dpgp_zygnema, " ")[[1]][2])
            highlight_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/zygnema_high_light.tsv"))) %>%
                right_join(dpgp_cluster_zygnema, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(highlight_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#ffea01") +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(0.25, 0.5)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(0.5, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(6, 6.25)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(6.25, 6.5)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(6.5, 7)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(7, 8)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(8, 10)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.25, 0.5, 2, 4, 6, 6.25, 6.5, 7, 8, 10)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        }
    }, width = 900, height = 900)
    output$physcomitrium_dpgp_expression_visualization <- renderPlot({
        if (str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][1] == "cold"){
            dpgp_cluster_physcomitrium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_physcomitrium) %>%
                filter(cluster == str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][2])
            cold_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/physcomitrium_cold.tsv"))) %>%
                right_join(dpgp_cluster_physcomitrium, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(cold_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#91e0ef") +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(0.5, 1)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(1, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                stat_smooth(data = subset(cold_FC_t_longer, Time_point %in% c(6, 24)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#023e8a",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.5, 1, 2, 4, 6, 24)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        } else if (str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][1] == "heat"){
            dpgp_cluster_physcomitrium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_physcomitrium) %>%
                filter(cluster == str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][2])
            heat_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/physcomitrium_heat.tsv"))) %>%
                right_join(dpgp_cluster_physcomitrium, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(heat_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#fdc3c3") +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(0.5, 1)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(1, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                stat_smooth(data = subset(heat_FC_t_longer, Time_point %in% c(6, 24)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#fe4848",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.5, 1, 2, 4, 6, 24)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        } else if (str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][1] == "highlight"){
            dpgp_cluster_physcomitrium <- read_tsv(normalizePath(file.path("./assets/DPGP/results/", paste0(input$dpgp_input_physcomitrium, ".txt")))) %>%
                filter(probability >= input$dpgp_probability_filter_threshold_physcomitrium) %>%
                filter(cluster == str_split(input$cluster_dpgp_physcomitrium, " ")[[1]][2])
            highlight_FC_t_longer <- read_tsv(normalizePath(file.path("./assets/DPGP/inputs/physcomitrium_high_light.tsv"))) %>%
                right_join(dpgp_cluster_physcomitrium, by = join_by("gene" == "gene")) %>%
                pivot_longer(!c(gene, cluster, probability), names_to = "Time_point", values_to = "logfc")
            p <- ggplot(highlight_FC_t_longer, aes(as.numeric(Time_point), y = logfc)) +
                geom_line(aes(group = gene), linewidth = 0.2, alpha = 0.3, color = "#ffea01") +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(0.25, 0.5)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(0.5, 2)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(2, 4)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(4, 6)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#ffa201",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(6, 6.25)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(6.25, 6.5)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(6.5, 7)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(7, 8)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                stat_smooth(data = subset(highlight_FC_t_longer, Time_point %in% c(8, 10)),
                            aes(x = as.numeric(Time_point), y = logfc), method = "loess", se = FALSE, linetype = "solid", color = "#70960e",  alpha = 1) +
                scale_y_continuous(limits = c(-2.1, 2.1)) +
                scale_x_continuous(breaks = c(0.25, 0.5, 2, 4, 6, 6.25, 6.5, 7, 8, 10)) +
                labs(y = "Z-score standardized log2 Fold Change (CPM)",
                     x = "Time (hours)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90),
                      axis.title = element_text(size = 14),     # Font size for axis titles
                      axis.text = element_text(size = 12),       # Font size for axis labels
                      plot.title = element_text(size = 16)       # Font size for plot title
                )
            p
        }
    }, width = 900, height = 900)
    # GRN
    output$download_mesotaenium_grn_table <- downloadHandler(
        filename = function() { paste0(input$grn_input_mesotaenium, "_", Sys.Date(), ".tsv") },
        content = function(file) { write_tsv(read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))), file) }
    )
    output$download_physcomitrium_grn_table <- downloadHandler(
        filename = function() { paste0(input$grn_input_physcomitrium, "_", Sys.Date(), ".tsv") },
        content = function(file) { write_tsv(read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))), file) }
    )
    output$download_zygnema_grn_table <- downloadHandler(
        filename = function() { paste0(input$grn_input_zygnema, "_", Sys.Date(), ".tsv") },
        content = function(file) { write_tsv(read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))), file) }
    )
    output$mesotaenium_grn_table <- DT::renderDataTable({
        if (input$grn_gene_filter_mesotaenium != ""){
            if (input$grn_filter_gene_category_mesotaenium == "regulator"){
                grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                    filter(regulator == input$grn_gene_filter_mesotaenium) %>%
                    slice_head(n = input$grn_top_n_mesotaenium)
            } else if (input$grn_filter_gene_category_mesotaenium == "target"){
                grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                    filter(target == input$grn_gene_filter_mesotaenium) %>%
                    slice_head(n = input$grn_top_n_mesotaenium)
            } else if (input$grn_filter_gene_category_mesotaenium == "both"){
                grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                    filter((regulator == input$grn_gene_filter_mesotaenium) | (target == input$grn_gene_filter_mesotaenium)) %>%
                    slice_head(n = input$grn_top_n_mesotaenium)
            }
        } else {
            grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                slice_head(n = input$grn_top_n_mesotaenium)
        }
        grn_mesotaenium %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$physcomitrium_grn_table <- DT::renderDataTable({
        if (input$grn_gene_filter_physcomitrium != ""){
            if (input$grn_filter_gene_category_physcomitrium == "regulator"){
                grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                    filter(regulator == input$grn_gene_filter_physcomitrium) %>%
                    slice_head(n = input$grn_top_n_physcomitrium)
            } else if (input$grn_filter_gene_category_physcomitrium == "target"){
                grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                    filter(target == input$grn_gene_filter_physcomitrium) %>%
                    slice_head(n = input$grn_top_n_physcomitrium)
            } else if (input$grn_filter_gene_category_physcomitrium == "both"){
                grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                    filter((regulator == input$grn_gene_filter_physcomitrium) | (target == input$grn_gene_filter_physcomitrium)) %>%
                    slice_head(n = input$grn_top_n_physcomitrium)
            }
        } else {
            grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                slice_head(n = input$grn_top_n_physcomitrium)
        }
        grn_physcomitrium %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$zygnema_grn_table <- DT::renderDataTable({
        if (input$grn_gene_filter_zygnema != ""){
            if (input$grn_filter_gene_category_zygnema == "regulator"){
                grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                    filter(regulator == input$grn_gene_filter_zygnema) %>%
                    slice_head(n = input$grn_top_n_zygnema)
            } else if (input$grn_filter_gene_category_zygnema == "target"){
                grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                    filter(target == input$grn_gene_filter_zygnema) %>%
                    slice_head(n = input$grn_top_n_zygnema)
            } else if (input$grn_filter_gene_category_zygnema == "both"){
                grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                    filter((regulator == input$grn_gene_filter_zygnema) | (target == input$grn_gene_filter_zygnema)) %>%
                    slice_head(n = input$grn_top_n_zygnema)
            }
        } else {
            grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                slice_head(n = input$grn_top_n_zygnema)
        }
        grn_zygnema %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
            ),
            extensions = 'Buttons',
            selection = 'single', ## enable selection of a single row
            filter = 'top',              ## include column filters at the top
            rownames = FALSE                ## don't show row numbers/names
            )
    })
    output$mesotaenium_grn_visualization <- renderVisNetwork({
        if (input$grn_gene_filter_mesotaenium != ""){
            if (input$grn_filter_gene_category_mesotaenium == "regulator"){
                grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                    filter(regulator == input$grn_gene_filter_mesotaenium) %>%
                    slice_head(n = input$grn_top_n_mesotaenium)
            } else if (input$grn_filter_gene_category_mesotaenium == "target"){
                grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                    filter(target == input$grn_gene_filter_mesotaenium) %>%
                    slice_head(n = input$grn_top_n_mesotaenium)
            } else if (input$grn_filter_gene_category_mesotaenium == "both"){
                grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                    filter((regulator == input$grn_gene_filter_mesotaenium) | (target == input$grn_gene_filter_mesotaenium)) %>%
                    slice_head(n = input$grn_top_n_mesotaenium)
            }
        } else {
            grn_mesotaenium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_mesotaenium, ".tsv")))) %>%
                slice_head(n = input$grn_top_n_mesotaenium)
        }
        mesotaenium_grn_network <- graph_from_data_frame(d=grn_mesotaenium, directed=T)
        mesotaenium_visNetowrk <- toVisNetworkData(mesotaenium_grn_network)
        visNetwork(nodes = mesotaenium_visNetowrk$nodes, edges = mesotaenium_visNetowrk$edges, height = "900px") %>%
            visIgraphLayout(physics = TRUE, layout = "layout_nicely") %>%
            visNodes(shape = "circle", size = rep(1, input$grn_top_n_mesotaenium),
                     color = list(background = "#a8dadc", 
                                  border = "#457b9d",
                                  highlight = "#ffc300"),
                     shadow = list(enabled = FALSE, size = 1)) %>%
            visEdges(shadow = FALSE,
                     arrows =list(to = list(enabled = TRUE, scaleFactor = 1)),
                     color = list(color = "#03045e", highlight = "#cad2c5")) %>%
            visOptions(highlightNearest = list(enabled = T, hover = T),
                       nodesIdSelection = T) %>% 
            visPhysics(solver = "repulsion", repulsion = list(centralGravity = 0, springLength = 200, springConstant = 0.1)) %>%
            visInteraction(navigationButtons = TRUE)
    })
    output$physcomitrium_grn_visualization <- renderVisNetwork({
        if (input$grn_gene_filter_physcomitrium != ""){
            if (input$grn_filter_gene_category_physcomitrium == "regulator"){
                grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                    filter(regulator == input$grn_gene_filter_physcomitrium) %>%
                    slice_head(n = input$grn_top_n_physcomitrium)
            } else if (input$grn_filter_gene_category_physcomitrium == "target"){
                grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                    filter(target == input$grn_gene_filter_physcomitrium) %>%
                    slice_head(n = input$grn_top_n_physcomitrium)
            } else if (input$grn_filter_gene_category_physcomitrium == "both"){
                grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                    filter((regulator == input$grn_gene_filter_physcomitrium) | (target == input$grn_gene_filter_physcomitrium)) %>%
                    slice_head(n = input$grn_top_n_physcomitrium)
            }
        } else {
            grn_physcomitrium <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_physcomitrium, ".tsv")))) %>%
                slice_head(n = input$grn_top_n_physcomitrium)
        }
        physcomitrium_grn_network <- graph_from_data_frame(d=grn_physcomitrium, directed=T)
        physcomitrium_visNetowrk <- toVisNetworkData(physcomitrium_grn_network)
        visNetwork(nodes = physcomitrium_visNetowrk$nodes, edges = physcomitrium_visNetowrk$edges, height = "900px") %>%
            visIgraphLayout(physics = TRUE, layout = "layout_nicely") %>%
            visNodes(shape = "circle", size = rep(1, input$grn_top_n_physcomitrium),
                     color = list(background = "#a8dadc", 
                                  border = "#457b9d",
                                  highlight = "#ffc300"),
                     shadow = list(enabled = FALSE, size = 1)) %>%
            visEdges(shadow = FALSE,
                     arrows =list(to = list(enabled = TRUE, scaleFactor = 1)),
                     color = list(color = "#03045e", highlight = "#cad2c5")) %>%
            visOptions(highlightNearest = list(enabled = T, hover = T),
                       nodesIdSelection = T) %>% 
            visPhysics(solver = "repulsion", repulsion = list(centralGravity = 0)) %>%
            visInteraction(navigationButtons = TRUE)
    })
    output$zygnema_grn_visualization <- renderVisNetwork({
        if (input$grn_gene_filter_zygnema != ""){
            if (input$grn_filter_gene_category_zygnema == "regulator"){
                grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                    filter(regulator == input$grn_gene_filter_zygnema) %>%
                    slice_head(n = input$grn_top_n_zygnema)
            } else if (input$grn_filter_gene_category_zygnema == "target"){
                grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                    filter(target == input$grn_gene_filter_zygnema) %>%
                    slice_head(n = input$grn_top_n_zygnema)
            } else if (input$grn_filter_gene_category_zygnema == "both"){
                grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                    filter((regulator == input$grn_gene_filter_zygnema) | (target == input$grn_gene_filter_zygnema)) %>%
                    slice_head(n = input$grn_top_n_zygnema)
            }
        } else {
            grn_zygnema <- read_tsv(normalizePath(file.path("./assets/SWING/", paste0(input$grn_input_zygnema, ".tsv")))) %>%
                slice_head(n = input$grn_top_n_zygnema)
        }
        zygnema_grn_network <- graph_from_data_frame(d=grn_zygnema, directed=T)
        zygnema_visNetowrk <- toVisNetworkData(zygnema_grn_network)
        visNetwork(nodes = zygnema_visNetowrk$nodes, edges = zygnema_visNetowrk$edges, height = "900px") %>%
            visIgraphLayout(physics = TRUE, layout = "layout_nicely") %>%
            visNodes(shape = "circle", size = rep(1, input$grn_top_n_zygnema),
                     color = list(background = "#a8dadc", 
                                  border = "#457b9d",
                                  highlight = "#ffc300"),
                     shadow = list(enabled = FALSE, size = 1)) %>%
            visEdges(shadow = FALSE,
                     arrows =list(to = list(enabled = TRUE, scaleFactor = 1)),
                     color = list(color = "#03045e", highlight = "#cad2c5")) %>%
            visOptions(highlightNearest = list(enabled = T, hover = T),
                       nodesIdSelection = T) %>% 
            visPhysics(solver = "repulsion", repulsion = list(centralGravity = 0)) %>%
            visInteraction(navigationButtons = TRUE)
    })
    # DEG
    output$logfc <- DT::renderDataTable({read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            datatable(options = list(paging = TRUE,    ## paginate the output
                                     pageLength = 10,  ## number of rows to output for each page
                                     lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                       c("5", "10", "25", "50", "100", 'All')),
                                     scrollX = TRUE,   ## enable scrolling on X axis
                                     scrollY = TRUE,   ## enable scrolling on Y axis
                                     rowCallback = JS(
                                         "function(row, data, index) {",
                                         "  $('td', row).css('vertical-align', 'top');",
                                         "}"
                                     ),
                                     autoWidth = FALSE, ## use smart column width handling
                                     server = TRUE,   ## use client-side processing
                                     processing = TRUE,
                                     dom = 'lfrtBip',
                                     buttons = c('csv', 'excel')
                                     ),
                      extensions = 'Buttons',
                      selection = 'single', ## enable selection of a single row
                      filter = 'top',              ## include column filters at the top
                      rownames = FALSE                ## don't show row numbers/names
    )
    })
    output$expression_tables_degs <- DT::renderDataTable({
        deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((abs(logFC) >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        }
    })
    output$expression_tables_up_regulated <- DT::renderDataTable({
        deg_gene_up_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((logFC >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_up_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_up_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_up_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        }
    })
    output$expression_tables_down_regulated <- DT::renderDataTable({
        deg_gene_down_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((logFC <= (log2(input$fold_change_threshold)*-1)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_down_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_down_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            read_tsv(normalizePath(file.path("./assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_down_IDs) %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 10,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        }
    })
    output$go_enrichment_up_regulated_table <- DT::renderDataTable({
        deg_gene_up_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((logFC >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            universe_up <- mds_mesotaenium_data_mesotaenium %>%
                right_join(GO_mesotaenium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_up <- GO_mesotaenium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            universe_up <- mds_physcomitrium_data_physcomitrium %>%
                right_join(GO_physcomitrium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_up <- GO_physcomitrium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            universe_up <- mds_zygnema_data_zygnema %>%
                right_join(GO_zygnema, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_up <- GO_zygnema
        }
        
        enrichment_up <- enricher(deg_gene_up_IDs,
                               TERM2GENE = GO_up %>%
                                   filter(domain == input$go_domain) %>%
                                   select(c(term, gene)),
                               TERM2NAME = GO_up %>%
                                   filter(domain == input$go_domain) %>%
                                   select(c(term, name)),
                               pvalueCutoff = input$go_pvalueCutoff,
                               universe = universe_up,
                               qvalueCutoff = input$go_qvalueCutoff)
        
        if (!is.null(enrichment_up@result)){
            enrichment_up@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 25,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$go_enrichment_down_regulated_table <- DT::renderDataTable({
        deg_gene_down_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((logFC <= (log2(input$fold_change_threshold)*-1)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            universe_down <- mds_mesotaenium_data_mesotaenium %>%
                right_join(GO_mesotaenium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_down <- GO_mesotaenium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            universe_down <- mds_physcomitrium_data_physcomitrium %>%
                right_join(GO_physcomitrium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_down <- GO_physcomitrium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            universe_down <- mds_zygnema_data_zygnema %>%
                right_join(GO_zygnema, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_down <- GO_zygnema
        }
        
        enrichment_down <- enricher(deg_gene_down_IDs,
                                  TERM2GENE = GO_down %>%
                                      filter(domain == input$go_domain) %>%
                                      select(c(term, gene)),
                                  TERM2NAME = GO_down %>%
                                      filter(domain == input$go_domain) %>%
                                      select(c(term, name)),
                                  pvalueCutoff = input$go_pvalueCutoff,
                                  universe = universe_down,
                                  qvalueCutoff = input$go_qvalueCutoff)
        
        if (!is.null(enrichment_down@result)){
            enrichment_down@result %>%
                datatable(options = list(paging = TRUE,    ## paginate the output
                                         pageLength = 25,  ## number of rows to output for each page
                                         lengthMenu = list(c(5, 10, 25, 50, 100, -1),
                                                           c("5", "10", "25", "50", "100", 'All')),
                                         scrollX = TRUE,   ## enable scrolling on X axis
                                         scrollY = TRUE,   ## enable scrolling on Y axis
                                         rowCallback = JS(
                                             "function(row, data, index) {",
                                             "  $('td', row).css('vertical-align', 'top');",
                                             "}"
                                         ),
                                         autoWidth = FALSE, ## use smart column width handling
                                         server = TRUE,   ## use client-side processing
                                         processing = TRUE,
                                         dom = 'lfrtBip',
                                         buttons = c('csv', 'excel')
                ),
                extensions = 'Buttons',
                selection = 'single', ## enable selection of a single row
                filter = 'top',              ## include column filters at the top
                rownames = FALSE                ## don't show row numbers/names
                )
        } else{
            tibble()
        }
    })
    output$expression_plots_degs <- renderPlot({
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((abs(logFC) >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs <- read_tsv(normalizePath(file.path("./assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((abs(logFC) >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs <- read_tsv(normalizePath(file.path("./assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((abs(logFC) >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs <- read_tsv(normalizePath(file.path("./assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv"))) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        }
    }, width = 1200, height = 800)
    output$expression_plots_up_regulated <- renderPlot({
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((logFC >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs_up <- read_tsv(normalizePath(file.path("./assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv"))) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs_up, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((logFC >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs_up <- read_tsv(normalizePath(file.path("./assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv"))) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs_up, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((logFC >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs_up <- read_tsv(normalizePath(file.path("./assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv"))) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs_up, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        }
    }, width = 1200, height = 800)
    output$expression_plots_down_regulated <- renderPlot({
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((logFC <= (log2(input$fold_change_threshold)*-1)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs_down <- read_tsv(normalizePath(file.path("./assets/gene_profiler/mesotaenium_count_voom_transformed_qsmooth.tsv"))) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs_down, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((logFC <= (log2(input$fold_change_threshold)*-1)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs_down <- read_tsv(normalizePath(file.path("./assets/gene_profiler/physcomitrium_count_voom_transformed_qsmooth.tsv"))) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs_down, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        } else if(str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            deg_gene_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
                filter((logFC <= (log2(input$fold_change_threshold)*-1)) & (adj.P.Val <= input$adjusted_p_value)) %>%
                select(geneID) %>% unlist() %>% as.vector()
            
            degs_down <- read_tsv(normalizePath(file.path("./assets/gene_profiler/zygnema_count_voom_transformed_qsmooth.tsv"))) %>%
                mutate_at(vars(-geneID), scale_z_score) %>%
                filter(geneID %in% deg_gene_IDs) %>%
                pivot_longer(cols = -geneID, names_to = "sample", values_to = "expression")
            ggplot(degs_down, aes(x = sample, y = geneID, fill = expression)) +
                geom_tile() +
                scale_fill_gradientn(colors = myheatcolors, breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), name = "Z-score") +
                labs(x = "Samples", y = "Genes", fill = "Z-Score Normalized expression") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10), axis.text.y = element_blank(), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
                coord_fixed(ratio = 1/100)
        }
    }, width = 1200, height = 800)
    output$go_ora_up_regulated <- renderPlot({
        deg_gene_up_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((logFC >= log2(input$fold_change_threshold)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            universe_up <- mds_mesotaenium_data_mesotaenium %>%
                right_join(GO_mesotaenium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_up <- GO_mesotaenium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            universe_up <- mds_physcomitrium_data_physcomitrium %>%
                right_join(GO_physcomitrium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_up <- GO_physcomitrium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            universe_up <- mds_zygnema_data_zygnema %>%
                right_join(GO_zygnema, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_up <- GO_zygnema
        }
        
        enrichment_up <- enricher(deg_gene_up_IDs,
                                  TERM2GENE = GO_up %>%
                                      filter(domain == input$go_domain) %>%
                                      select(c(term, gene)),
                                  TERM2NAME = GO_up %>%
                                      filter(domain == input$go_domain) %>%
                                      select(c(term, name)),
                                  pvalueCutoff = input$go_pvalueCutoff,
                                  universe = universe_up,
                                  qvalueCutoff = input$go_qvalueCutoff)
        if (!is.null(enrichment_up@result)){
            p <- dotplot(enrichment_up,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
        }, width = 1200, height = 1000)
    output$go_ora_down_regulated <- renderPlot({
        deg_gene_down_IDs <- read_tsv(normalizePath(file.path("./assets/dgea/fold_change_tables/", paste0(input$dgea_comparison, "_logFC_ebFit.tsv")))) %>%
            filter((logFC <= (log2(input$fold_change_threshold)*-1)) & (adj.P.Val <= input$adjusted_p_value)) %>%
            select(geneID) %>% unlist() %>% as.vector()
        
        if (str_split(input$dgea_comparison, "_")[[1]][1] == "mesotaenium"){
            universe_down <- mds_mesotaenium_data_mesotaenium %>%
                right_join(GO_mesotaenium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_down <- GO_mesotaenium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "physcomitrium"){
            universe_down <- mds_physcomitrium_data_physcomitrium %>%
                right_join(GO_physcomitrium, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_down <- GO_physcomitrium
        } else if (str_split(input$dgea_comparison, "_")[[1]][1] == "zygnema"){
            universe_down <- mds_zygnema_data_zygnema %>%
                right_join(GO_zygnema, by = join_by("geneID" == "gene")) %>%
                select(geneID) %>%
                unlist() %>%
                as.vector()
            GO_down <- GO_zygnema
        }
        
        enrichment_down <- enricher(deg_gene_down_IDs,
                                    TERM2GENE = GO_down %>%
                                        filter(domain == input$go_domain) %>%
                                        select(c(term, gene)),
                                    TERM2NAME = GO_down %>%
                                        filter(domain == input$go_domain) %>%
                                        select(c(term, name)),
                                    pvalueCutoff = input$go_pvalueCutoff,
                                    universe = universe_down,
                                    qvalueCutoff = input$go_qvalueCutoff)
        if (!is.null(enrichment_down@result)){
            p <- dotplot(enrichment_down,
                         x= "geneRatio", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         color="p.adjust",
                         orderBy = "x", # Options: GeneRatio, BgRatio, pvalue, p.adjust, qvalue
                         showCategory=input$go_show_category,
                         font.size=12) +
                ggtitle("Dotplot of ORA")
            print(p)
        } else {
            print(ggplot()) # if nothing is enriched!
        }
    }, width = 1200, height = 1000)
    output$mesotaenium_dendrogram <- renderPlot({
        ggtree(clusters_mesotaenium_OTU, aes(color=group)) +
            scale_color_manual(values = my_colors_mesotaenium) +
            geom_tiplab(size=3) + 
            theme(legend.position='none', plot.margin = margin(r = 1, unit = "cm"))
    }, width = 1500, height = 1000)
    output$zygnema_dendrogram <- renderPlot({
        ggtree(clusters_zygnema_OTU, aes(color=group)) +
            scale_color_manual(values = my_colors_zygnema) +
            geom_tiplab(size=3) + 
            theme(legend.position='none', plot.margin = margin(r = 1, unit = "cm"))
    }, width = 1500, height = 1000)
    output$physcomitrium_dendrogram <- renderPlot({
        ggtree(clusters_physcomitrium_OTU, aes(color=group)) +
            scale_color_manual(values = my_colors_physcomitrium) +
            geom_tiplab(size=3) + 
            theme(legend.position='none', plot.margin = margin(r = 1, unit = "cm"))
    }, width = 1500, height = 1000)
    output$mesotaenium_plot_PCA <- renderPlot({
        ggplot(pca.res_mesotaenium.df) +
                aes(x=PC1, y=PC2, label=mesotaenium_sampleLabels, shape = mesotaenium_replicate, color = mesotaenium_condition) +
                geom_point(size=4) +
                xlab(paste0("PC1 (",pc.per[1],"%",")")) +
                ylab(paste0("PC2 (",pc.per[2],"%",")")) +
                scale_color_manual(values = my_colors_mesotaenium) +
                scale_fill_brewer(palette="Set3") +
                labs(title=paste0("PCA plot of filtered and normalized via qsmooth")) +
                coord_fixed() +
                theme_bw() + 
            theme(
                axis.title = element_text(size = 14),     # Font size for axis titles
                axis.text = element_text(size = 12),       # Font size for axis labels
                plot.title = element_text(size = 16),       # Font size for plot title
                legend.text = element_text(size = 12)       # Font size for plot legend
            )
    }, width = 1000, height = 600)
    output$zygnema_plot_PCA <- renderPlot({
        ggplot(pca.res_zygnema.df) +
            aes(x=PC1, y=PC2, label=zygnema_sampleLabels, shape = zygnema_replicate, color = zygnema_condition) +
            geom_point(size=4) +
            xlab(paste0("PC1 (",pc.per[1],"%",")")) +
            ylab(paste0("PC2 (",pc.per[2],"%",")")) +
            scale_color_manual(values = my_colors_zygnema) +
            scale_fill_brewer(palette="Set3") +
            labs(title=paste0("PCA plot of filtered and normalized via qsmooth")) +
            coord_fixed() +
            theme_bw() + 
            theme(
                axis.title = element_text(size = 14),     # Font size for axis titles
                axis.text = element_text(size = 12),       # Font size for axis labels
                plot.title = element_text(size = 16),       # Font size for plot title
                legend.text = element_text(size = 12)       # Font size for plot legend
            )
    }, width = 1000, height = 600)
    output$physcomitrium_plot_PCA <- renderPlot({
        ggplot(pca.res_physcomitrium.df) +
            aes(x=PC1, y=PC2, label=physcomitrium_sampleLabels, shape = physcomitrium_replicate, color = physcomitrium_condition) +
            geom_point(size=4) +
            xlab(paste0("PC1 (",pc.per[1],"%",")")) +
            ylab(paste0("PC2 (",pc.per[2],"%",")")) +
            scale_color_manual(values = my_colors_physcomitrium) +
            scale_fill_brewer(palette="Set3") +
            labs(title=paste0("PCA plot of filtered and normalized via qsmooth")) +
            coord_fixed() +
            theme_bw() + 
            theme(
                axis.title = element_text(size = 14),     # Font size for axis titles
                axis.text = element_text(size = 12),       # Font size for axis labels
                plot.title = element_text(size = 16),       # Font size for plot title
                legend.text = element_text(size = 12)       # Font size for plot legend
            )
    }, width = 1000, height = 600)
    output$mesotaenium_plotMDS <- renderPlot({
        ggplot(df_mesotaenium, aes(x=mds_mesotaenium.x, y=mds_mesotaenium.y, label=Sample, shape = mesotaenium_replicate, color = mesotaenium_condition)) +
            geom_point(size=4) +
            xlab(paste0("Leading logFC dim 1 (", round(mds_mesotaenium$var.explained[1] * 100, 2),"%)")) +
            ylab(paste0("Leading logFC dim 2 (", round(mds_mesotaenium$var.explained[2] * 100, 2),"%)")) +
            scale_color_manual(values = my_colors_mesotaenium) +
            labs(title=paste0("MDS plot of filtered and qsmooth normalized")) +
            coord_fixed() +
            theme_bw() + 
            theme(
                axis.title = element_text(size = 14),     # Font size for axis titles
                axis.text = element_text(size = 12),       # Font size for axis labels
                plot.title = element_text(size = 16),       # Font size for plot title
                legend.text = element_text(size = 12)       # Font size for plot legend
            )
    }, width = 1000, height = 600)
    output$zygnema_plotMDS <- renderPlot({
        ggplot(df_zygnema, aes(x=mds_zygnema.x, y=mds_zygnema.y, label=Sample, shape = zygnema_replicate, color = zygnema_condition)) +
            geom_point(size=4) +
            xlab(paste0("Leading logFC dim 1 (", round(mds_zygnema$var.explained[1] * 100, 2),"%)")) +
            ylab(paste0("Leading logFC dim 2 (", round(mds_zygnema$var.explained[2] * 100, 2),"%)")) +
            scale_color_manual(values = my_colors_zygnema) +
            labs(title=paste0("MDS plot of filtered and qsmooth normalized")) +
            coord_fixed() +
            theme_bw() + 
            theme(
                axis.title = element_text(size = 14),     # Font size for axis titles
                axis.text = element_text(size = 12),       # Font size for axis labels
                plot.title = element_text(size = 16),       # Font size for plot title
                legend.text = element_text(size = 12)       # Font size for plot legend
            )
    }, width = 1000, height = 600)
    output$physcomitrium_plotMDS <- renderPlot({
        ggplot(df_physcomitrium, aes(x=mds_physcomitrium.x, y=mds_physcomitrium.y, label=Sample, shape = physcomitrium_replicate, color = physcomitrium_condition)) +
            geom_point(size=4) +
            xlab(paste0("Leading logFC dim 1 (", round(mds_physcomitrium$var.explained[1] * 100, 2),"%)")) +
            ylab(paste0("Leading logFC dim 2 (", round(mds_physcomitrium$var.explained[2] * 100, 2),"%)")) +
            scale_color_manual(values = my_colors_physcomitrium) +
            labs(title=paste0("MDS plot of filtered and qsmooth normalized")) +
            coord_fixed() +
            theme_bw() + 
            theme(
                axis.title = element_text(size = 14),     # Font size for axis titles
                axis.text = element_text(size = 12),       # Font size for axis labels
                plot.title = element_text(size = 16),       # Font size for plot title
                legend.text = element_text(size = 12)       # Font size for plot legend
            )
    }, width = 1000, height = 600)
    # Add custom CSS for top alignment of images
    output$custom_css <- renderUI({
        tags$head(
            tags$style(HTML("
        #top_aligned_image img {
          vertical-align: top;
        }
      "))
        )
    })
}
################################################################################
############### Start the webapp
################################################################################
shiny::shinyApp(ui, server)