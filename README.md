# BRStudio
**An automated NGS data analyzing tool from [Braunstein Lab](http://www.braunstein.team/).**  
Current version: v1.0.2

## Setting up the environment
- This script has been tested under ubuntu 19.10. To set up the environment with terminal:  
```
  sudo apt update && sudo apt install r-base r-cran-rjava -y  
  sudo apt install libxml2-dev libssl-dev libcurl4-openssl-dev libxml2-dev -y  
  sudo -i R  
  install.packages(c("expss", "XML", "vcfR", "foreach", "doParallel", "R.utils", "RSelenium","readr", "dplyr", "tidyr", "filesstrings","readxl", "writexl", "wdman", "stringr", "stringi", "rlang"))  
```
- If you are using it under other systems, please make sure that all packages are properly installed.  

## Running the script
- This script exploits multiple threads for parallel computing.  
- This script automatically retrives variant info from [GnomAD](https://gnomad.broadinstitute.org/), [NCBI SNP](https://www.ncbi.nlm.nih.gov/snp/), [Provean](http://provean.jcvi.org/index.php) and [PolyPhen](http://genetics.bwh.harvard.edu/pph2/bgi.shtml).  
- You may run it with multiple panel files, but all of which should be formatted to 4 columns: "Gene", "Transcript ID", "Protein ID", "Protein ID for PP2".  
- You may run it with multiple pool files, but all of which should be formatted to 2 columns: "Sample", "Diagnosis".    
- This script will only process ".genome.vcf" files. If other vcf files are found, they will not be processed but only be archived.  
- Running no more than 96 samples in each batch is recommended.  
- For each batch, 20 to 30 minutes might be required, please make sure that your computer/server won't go to sleep.  
- Put ".panel.xlsx", ".pool.xlsx", and ".vcf.gz" files in the same folder with this script.  
- Open the script directly from the parent folder with RStudio to start.  
