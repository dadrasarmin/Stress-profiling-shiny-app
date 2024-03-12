# Introduction

We utilized principal component analysis (PCA) to extract the core of the data and find the most significant sources of variation between samples. Multidimensional Scaling (MDS) goes one step further, allowing us to see the subtle interactions between samples, discovering underlying structures that would otherwise go undiscovered. Finally, our Hierarchical Clustering reveals natural data groupings, allowing us to identify physiologically significant groups.

# Methods

We employed a range of R packages for our exploratory data analysis, and you can find the comprehensive list in the Methods section. To give an overview, our process involved importing `Kallisto` quantification files using the `txImport` tool and applying a `lengthScaledTPM` transformation. We then filtered the data, retaining reads with a minimum Count-Per-Million (CPM) of at least 10 across a minimum of three samples. Ensuring robustness, we conducted Smooth Quantile Normalization, accounting for varying experimental conditions, and subsequently transformed the data using the `voom` function from the `limma` package. To perform hierarchical clustering, we computed distances using the Euclidean method and opted for the "ward.D"" method for agglomeration.

# A note on abbreviations

Each treatment has a name composed of two components separated by an underscore. The first portion is the treatment, and the second part is the time (in hours) since the start of the experiment. We had two highlight settings: stress ("s") and recovery ("r").