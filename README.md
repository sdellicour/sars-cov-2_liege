This repo gathers all the input files and scripts related to our study entitled "Exploiting genomic surveillance to map the spatio-temporal dispersal of SARS-CoV-2 spike mutations in Belgium across 2020" (Bollen *et al*. *submitted*): R scripts to prepare, analyse, and visualise the phylogeographic analyses performed in this study, Geographic Information System (GIS) files used in those scripts, as well as resulting XML files for running the phylogeographic analyses with the program BEAST 1.10.

At the end of 2020, several new variants of SARS-CoV-2 - designated variants of concern - were detected and quickly suspected to be associated with a higher transmissibility and possible escape of vaccine-induced immunity. In Belgium, this discovery has motivated the initiation of a more ambitious genomic surveillance program, which is drastically increasing the number of SARS-CoV-2 genomes to analyse for monitoring the circulation of viral lineages and variants of concern. In order to efficiently analyse the massive collection of genomic data that are the result of such increased sequencing efforts, streamlined analytical strategies are crucial. In our study, we illustrate how to efficiently map the spatio-temporal dispersal of target mutations at a regional level. As a proof of concept, we focus on the Belgian province of Liège that has been consistently sampled throughout 2020, but was also one of the main epicenters of the second European epidemic wave.

Here, we summarise the four main steps of the analytical pipeline implemented for this study.

## 1. Sequence data collection
The first step consists in gathering fine-grained genomic sequences collected within the geographic area of interest (i.e. the province of Liège in our case), but also reference genomic sequences collected outside the target area and representing the genetic diversity of the virus in this region of the world. We included 2,114 Belgian sequences (including 869 from the province of Liège) as well as 3,543 reference sequences used in the European Nextstrain build on December 1, 2020. Once all the viral genomes are gathered, they have to be aligned (in our study, we used the program [MAFFT](https://mafft.cbrc.jp/alignment/software/source.html)).

How to download reference sequences using [Nextstrain](https://nextstrain.org/) and [GISAID](https://www.gisaid.org/):
- Select the Nextstrain build that provides enough phylogenetic context for the target geographic area. In our study, we selected the European Nextstrain build: https://nextstrain.org/ncov/gisaid/europe?f_region=Europe.
- Scroll to the end of the webpage and click on “DOWNLOAD DATA”, then click on “AUTHOR METADATA (TSV)”.
- The metadata file contains the EPI ISL ID’s for all the sequences. Open the GISAID EpiCoV database. Click on “Select” and copy-paste the list of EPI ISL ID’s, then download the selected sequences.

## 2. Inferring a time-scaled phylogenetic tree
The second step consists in inferring a time-calibrated phylogenetic tree. In our pipeline, we use a maximum likelihood (ML) rather than a Bayesian method, which would be notably time-consuming with such a large data set.

How to infer a time-calibrated ML tree using [IQ-TREE](http://www.iqtree.org/) and [TreeTime](https://github.com/neherlab/treetime):
- Run an IQ-TREE analysis based on the genomic sequence alignment obtained in the previous step. Follow the beginner’s tutorial found [here](http://www.iqtree.org/doc/Tutorial) or simply use the following command in your Terminal: iqtree -s alignment.phy -m MFP
- It is possible that IQ-TREE flagged some sequences as outliers. Remove these sequences from your tree file (named “alignment.phy.treefile” by default) before moving forward.
- Run a TreeTime analysis based on the maximum likelihood phylogenetic tree file obtained from the IQ-TREE analysis, along with the alignment file and a metadata file containing the names and dates of all genomic sequences. A tutorial can be found [here](https://treetime.readthedocs.io/en/latest/tutorials.html), or use the following command line in your Terminal: treetime --tree alignment.phy.treefile --aln alignment.phy --dates dates.tsv
- TreeTime might also have flagged some sequences as outliers, which also have to be removed from the tree file before continuing. 

## 3. Discrete phylogeographic analysis
The second step consists in performing a preliminary discrete phylogeographic analysis. The aim of this step is to provide phylogenetic context for the fine-grained analysis performed in the final step. Specifically, we want to identify the individual introduction events into the area of interest (province of Liège in our case).

How to construct a discrete phylogeography and identify introduction events, using [BEAST 1.10](https://beast.community/) and R:
- In order to run this preliminary discrete phylogeographic analysis, we follow the first two steps (**steps 1 and 2**) outlined in the R script file associated with this GitHub repository (**Phylogenetic_analyses.r**). A tutorial for running a discrete trait phylogeography in BEAST can also be found here: https://beast.community/workshop_discrete_diffusion. Note that in the “traits” tab, only two possible locations are considered: “location of interest” (province of Liège) or “other”. Importantly, for this analysis, we use the tree obtained in the previous step as a fixed empirical tree.
- Check that the MCMC chain has reached adequate mixing and convergence using the program [Tracer](http://tree.bio.ed.ac.uk/software/tracer/) (part of the BEAST utility package). An Effective Sample Size (ESS) value of >200 for all estimated parameters is often taken as the threshold. A tutorial on how to use Tracer can be found [here]( https://beast.community/analysing_beast_output).
- Construct a maximum clade credibility (MCC) tree using the program [TreeAnnotator](https://beast.community/treeannotator) (part of the BEAST utility package), discarding 10% (or more) of the initial states as burn-in. A tutorial on how to use TreeAnnotator can be found [here](https://beast.community/treeannotator).
- Extract all clades whose estimated origin is the location of interest using R. The code for this can be found in the **Phylogenetic_analyses.r** script. Follow **step 3** named “Analysing the outputs of the preliminary discrete phylogeographic analysis” and **step 4** named “Identifying the different clusters (clades following introduction events)”.

## 4. Continuous phylogeographic analysis
Finally, we can run a continuous phylogeographic analysis along the set of clades identified in the previous step and corresponding to distinct introduction events. In our study, we also take advantage of the continuous phylogeographic reconstruction to track specific SARS-CoV-2 spike protein mutations. This part of the analysis relies heavily on the **Phylogenetic_analyses.r** script.

How to construct continuous phylogeographic trees using [BEAST 1.10](https://beast.community/) and R:
- Generate an XML file to perform a continuous phylogeographic reconstruction along each clade identified during the previous step and corresponding to a distinct introduction event into the target geographic area. For this purpose, we can use R code in **Phylogenetic_analyses.r** (**step 5**). We can then run this XML file using BEAST and construct MCC trees (**step 6**).
- Extract the spatiotemporal information embedded in spatially-annotated trees using the “treeExtractions” function of the R package “seraphim” (**step 7**). This information can then be used to generate the map figures reported in our study (**step 8** in the script).
- Track the dispersal of mutation of interest (**step 9** in the script).
