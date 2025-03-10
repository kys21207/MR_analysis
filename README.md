# MR_analysis
## install MR packages based on DNAnexus
1. First, create the environment.yml file <br>
2. conda env create -f environment.yml <br>
3. conda activate mr_analysis <br>
4. Next, open R within the activated conda environment and install the missing packages: <br>
- Install TwoSampleMR
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
- Install MRPRESSO from GitHub
devtools::install_github("rondolab/MR-PRESSO")
