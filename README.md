# Stress profiling website

## Setup

In order to install the required packages please install the following packages using `microconda` or `conda` or any other way that you like. Please do not change the version of any software.

```
micromamba create -n stress_profiling_website -c conda-forge -c bioconda r-tidyverse=2.0.0 r-ggforce=0.4.1 r-ggdist=3.3.1 r-gghalves=0.1.4 bioconductor-treeio=1.26.0 bioconductor-ggtree=3.10.0 r-dt=0.31 r-shiny=1.8.0 r-shinydashboard=0.7.2 r-htmltools=0.5.7 r-msar=0.6.0 r-ape=5.7_1 bioconductor-limma=3.58.1 bioconductor-clusterprofiler=4.10.0 bioconductor-enrichplot=1.22.0 r-visnetwork=2.1.2 r-igraph=1.6.0 r-base=4.3.2 r-markdown=1.12
```

Then, you can activate the environment as follows:

```
micromamba activate stress_profiling_website
```

You need to clone this repository:

```
git clone https://gitlab.gwdg.de/applbionf/stress_profiling_website.git
```

## Download the datasets

You need to download the datasets and put them in the right location in order to get a functioning Shiny app. It is available (here on Zenodo)[https://zenodo.org/records/10808881]. To do that follow these steps.
### Download the compressed file

```
curl https://zenodo.org/api/records/10808881/files-archive -o assets_zenodo.zip
```
### Decompress it

```
unzip assets_zenodo.zip
unzip assets.zip
```

Then, you can move the `assets` folder to the git project folder (if it is not already there).

### Clean up

In order to save space, you can remove the compressed files as follows:

```
rm -rf assets_zenodo.zip assets.zip
``` 

## Run the app

When you are in the repository on your local machine, execute the following command and when it is ready, you will get an IP. Enter the IP into your browser to explore the website.

```
Rscript stress_profiling_website.R
```

I added the exported environment as a yaml file to this repository (`micromamba_environment.yaml`).
