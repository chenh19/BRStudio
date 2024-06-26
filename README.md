# BRStudio
**An automated NGS data analyzing tool from Braunstein Lab.**  
*Current version: v3.1.0*  
  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10989791.svg)](https://zenodo.org/doi/10.5281/zenodo.6632037)


## Setting up the environment
- [BRStudio](https://github.com/chenh19/BRStudio/blob/master/BRStudio.R) has been tested under [Debian 12.5.0](https://cdimage.debian.org/debian-cd/current-live/amd64/iso-hybrid/). 
- To set up the environment, please refer to [MyWorkspace](https://github.com/chenh19/MyWorkspace). You might need ```google-chrome-stable_114.0.5735.90-1_amd64.deb``` package, which can be downloaded from [Google-1](https://dl.google.com/linux/deb/pool/main/g/google-chrome-stable/google-chrome-stable_114.0.5735.90-1_amd64.deb) / [Google-2](https://dl.google.com/linux/chrome/deb/pool/main/g/google-chrome-stable/google-chrome-stable_114.0.5735.90-1_amd64.deb) / [UChicago mirror](http://mirror.cs.uchicago.edu/google-chrome/pool/main/g/google-chrome-stable/google-chrome-stable_114.0.5735.90-1_amd64.deb).
- If you are using it under other systems, please make sure that all relevant packages are properly installed.  

## Running the script
- This script exploits multiple threads for parallel computing.  
- This script automatically retrives variant info from [GnomAD](https://gnomad.broadinstitute.org/), [NCBI SNP](https://www.ncbi.nlm.nih.gov/snp/), [Provean](http://provean.jcvi.org/index.php) and [PolyPhen](http://genetics.bwh.harvard.edu/pph2/bgi.shtml), please make sure that you are connected to the internet.    
- You may run it with multiple panel files, but all of which should be formatted to 5 columns: "Gene", "Transcript_ID", "Protein_ID", "Protein_ID_for_PP2", "ENSG_ID"  (refer to [examples](https://github.com/chenh19/BRStudio/tree/master/examples)).  
- You may run it with multiple pool files, but all of which should be formatted to 2 columns: "Sample", "Diagnosis" (refer to [examples](https://github.com/chenh19/BRStudio/tree/master/examples)).  
- This script will only process ".genome.vcf" files. If other vcf files are found, they will not be processed but only be archived.  
- Running no more than 200 samples in each batch is recommended.  
- For each batch, 30 tp 60 minutes might be required, please make sure that your computer/server won't go to sleep.  
- Put ".panel.xlsx", ".pool.xlsx", and ".vcf.gz" files in the same folder with this script.  
- Launch RStudio, go to **File** menu and **Open File** > select the script, then go to **Session** menu and **Set Working Directory** > **To Source File Location**.  
- Run the script.
