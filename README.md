# MR_analysis
## Install MR packages based on DNAnexus
1. First, create the environment.yml file <br>
2. conda env create -f environment.yml <br> ** Take a long time to install **
3. conda activate mr_analysis <br>
4. Next, open R within the activated conda environment and install the missing packages: <br>
- Install TwoSampleMR
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
- Install MRPRESSO from GitHub
devtools::install_github("rondolab/MR-PRESSO") <br>
** take a long time to install dependency **
## Note
- Coloc answers: "Do gene expression and the trait share the same causal variant(s)?" <br>
- MR answers: "Does gene expression cause a change in the trait? And if so, how large is that effect?" <br>
Often, fine-mapping plus colocalization will yield a single "top" variant (or a single CS in tight LD). In that scenario, we can't do a multi-instrument (IVW) MR becasue we lack multiple independent SNPs.
