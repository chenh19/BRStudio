# BRStudio
**An automated NGS data analyzing tool from [Braunstein Lab](https://www.braunstein.team/).**  
*Current version: v3.0.0*

## Setting up the environment
- [BRStudio](https://github.com/chenh19/BRStudio/blob/master/BRStudio.R) has been tested under [Kubuntu 22.04 LTS](https://kubuntu.org/) and linux systems based on Ubuntu 22.04 LTS. 
- To set up the environment, simply execute [setup.sh](https://github.com/chenh19/BRStudio/blob/master/setup.sh) with ```bash ./setup.sh``` in terminal.
- If you are using it under other systems, please make sure that all relevant packages are properly installed.  

## Running the script
- This script exploits multiple threads for parallel computing.  
- This script automatically retrives variant info from [GnomAD](https://gnomad.broadinstitute.org/), [NCBI SNP](https://www.ncbi.nlm.nih.gov/snp/), [Provean](http://provean.jcvi.org/index.php) and [PolyPhen](http://genetics.bwh.harvard.edu/pph2/bgi.shtml), please make sure that you are connected to the internet.    
- You may run it with multiple panel files, but all of which should be formatted to 4 columns: "Gene", "Transcript ID", "Protein ID", "Protein ID for PP2".  
- You may run it with multiple pool files, but all of which should be formatted to 2 columns: "Sample", "Diagnosis".    
- This script will only process ".genome.vcf" files. If other vcf files are found, they will not be processed but only be archived.  
- Running no more than 200 samples in each batch is recommended.  
- For each batch, 30 tp 60 minutes might be required, please make sure that your computer/server won't go to sleep.  
- Put ".panel.xlsx", ".pool.xlsx", and ".vcf.gz" files in the same folder with this script.  
- Launch RStudio, go to **File** menu and **Open File...** > select the script, then go to **Session** menu and **Set Working Directory** > **To Source File Location**.  
- Run the script.
