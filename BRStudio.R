# BRStudio v3.0.2
# This script exploits multiple threads for parallel computing
# This script has been tested under Kubuntu 22.04 LTS. If you are using it under other systems, please make sure that all packages are properly installed
# You may run it with multiple panel files, but all of which should be formatted to 5 columns: Gene, Transcript_ID, Protein_ID, Protein_ID_for_PP2, ENSG_ID
# You may run it with multiple pool files, but all of which should be formatted to 2 columns: Sample, Diagnosis
# This script will only process ".genome.vcf" files. If other vcf files are found, they will not be processed but only be archived
# Running no more than 200 samples in each batch is recommended
# For each batch, 30 to 60 minutes might be required, please make sure that your computer/server won't go to sleep
# Put ".panel.xlsx", ".pool.xlsx", and ".vcf.gz" files in the same folder with this script
# Launch RStudio, go to File menu and Open File > select the script, then go to Session menu and Set Working Directory > To Source File Location
# Run the script


# Load R packages
lapply(c("XML", "vcfR", "parallel", "foreach", "doParallel", "R.utils", "RSelenium", "expss", "rlang", "readr", 
               "dplyr", "tidyr", "filesstrings", "readxl", "writexl", "wdman", "stringr", "stringi"), 
             require, character.only = TRUE)
print("Welcome to BRStudio, the system is loading...")


# Set the folder containing current R script as working directory
setwd(".")
print(paste("Current working dir:", getwd()))
if (dir.exists("original_files")==FALSE){
  dir.create("original_files")
}
if (dir.exists("cache")==FALSE){
  dir.create("cache")
}


# Set cpu cores for parallel computing
numCores <- detectCores(all.tests = FALSE, logical = TRUE)
print(paste("Parallel computing:", numCores, "cores will be used for data processing"))
Sys.sleep(1)


# Load gene panel files, remove duplicates
panelfiles <- list.files(pattern='*.panel.xlsx')
if (is_empty(panelfiles)==FALSE){
  panel <- lapply(panelfiles, read_excel)
  panel <- bind_rows(panel)
  panel <- panel %>% 
    group_by(Gene) %>% 
    summarise_all(funs(trimws(paste(., collapse = ''))))
  panel$"Transcript_ID" <- gsub('NA', '', panel$"Transcript_ID")
  panel$"Protein_ID" <- gsub('NA', '', panel$"Protein_ID")
  panel$"Protein_ID_for_PP2" <- gsub('NA', '', panel$"Protein_ID_for_PP2")
  panel$"ENSG_ID" <- gsub('NA', '', panel$"ENSG_ID")
  panel <- panel[c("Gene", "Transcript_ID", "Protein_ID", "Protein_ID_for_PP2", "ENSG_ID")]
  write.table(panel, file="all.panel.csv", sep=",", row.names = FALSE)
  registerDoParallel(numCores)
  foreach (panelfile = panelfiles) %dopar% {
    file.move(panelfile, "./original_files/a_panel_files", overwrite=TRUE)
  }
  genes <- read.csv("all.panel.csv")[,"Gene"]
  #ensg <- read.csv("all.panel.csv")[,"ENSG_ID"]
  genes <- split(genes, ceiling(seq_along(genes)/5))
  #ensg <- split(genes, ceiling(seq_along(genes)/5))
  
  panel <- read.csv("all.panel.csv")
  file.move("all.panel.csv", "./cache/a_gene_panel", overwrite=TRUE)
  print(paste("All panels loaded:", nrow(panel), "genes in total"))
}
rm("panelfiles")
Sys.sleep(1)


# Load pool files, remove duplicates
poolfiles <- list.files(pattern='*.pool.xlsx')
if (is_empty(poolfiles)==FALSE){
  pool <- lapply(poolfiles, read_excel)
  pool <- bind_rows(pool)
  pool <- pool %>% 
    group_by(Sample) %>% 
    summarise_all(funs(trimws(paste(., collapse = ''))))
  write.table(pool, file="all.pool.csv", sep=",", row.names = FALSE)
  registerDoParallel(numCores)
  foreach (poolfile = poolfiles) %dopar% {
    file.move(poolfile, "./original_files/b_pool_files", overwrite=TRUE)
  }
  pool <- read.csv("all.pool.csv")
  file.move("all.pool.csv", "./cache/b_sample_pool", overwrite=TRUE)
  print(paste("Sample pools loaded:", nrow(pool), "samples in total"))
}
rm("poolfiles")
Sys.sleep(1)


# Retrieve variant info from GnomAD v2.1
print("Downloading ref library for each gene...(this might take a while)")
for (i in genes){
  current_ver <- sapply(binman::list_versions("chromedriver"),"[[",1)
  cDrv <- chrome(version = current_ver, verbose = FALSE)
  eCaps <- list(
    chromeOptions = list(
      args = c('--headless', '--disable-gpu', '--window-size=1280,800'),
      prefs = list(
        "profile.default_content_settings.popups" = 0L,
        "download.prompt_for_download" = FALSE,
        "download.default_directory" = "."
      )
    )
  )
  remDr <- remoteDriver(browserName = "chrome", port = 4567L, extraCapabilities = eCaps)
  remDr$open()
  remDr$setTimeout(type = "implicit", milliseconds = 50000)
  remDr$setTimeout(type = "page load", milliseconds = 50000)
  
  for (gene in i){
    ensg_id <- panel[panel$Gene == gene,]$ENSG_ID
    geneurl <- paste0("https://gnomad.broadinstitute.org/gene/", ensg_id)
    remDr$navigate(geneurl) # remDr$getCurrentUrl()
    csv_link <- remDr$findElement(using="xpath", "//button[contains(.,'Export variants to CSV')]")
    Sys.sleep(2)
    csv_link$clickElement()
    Sys.sleep(2)
    reffile <- list.files(pattern=".csv")
    if (is_empty(reffile)==FALSE){
      filename <- paste0(gene, ".ref")
      file.rename(reffile, filename)
    }
  }
  remDr$navigate("https://gnomad.broadinstitute.org/")
  remDr$close()
  cDrv$stop()
  Sys.sleep(10)
}
reffiles <- list.files(pattern=".ref")
registerDoParallel(numCores)
foreach (reffile = reffiles) %dopar% {
  filename <- paste0(reffile, ".csv")
  file.rename(reffile, filename)
}
rm("cDrv", "eCaps", "remDr", "gene", "genes", "geneurl", "filename", "csv_link", "i","ensg_id", "current_ver")
print("Ref library for each gene is downloaded")
Sys.sleep(10)


# Organize variant ref info
nans <- read.csv("NANS.ref.csv", header = TRUE, na.strings=c("","NA"))
ash_jew <- c("Allele_Count_Ashkenazi_Jewish", "Allele_Number_Ashkenazi_Jewish", "Homozygote_Count_Ashkenazi_Jewish", "Hemizygote_Count_Ashkenazi_Jewish")
nans[ , ash_jew] <- 0
nans[,"Allele_Number_Ashkenazi_Jewish"] <-1
colnames(nans) <- c("Chromosome", "Position", "rsID", "Reference", "Alternate", "Source", "Filters_exomes", "Filters_genomes", "Transcript",
                    "HGVS_Consequence", "Protein_Consequence", "Transcript_Consequence", "VEP_Annotation", "ClinVar_Clinical_Significance", "ClinVar_Variation_ID", "Flags", 
                    "Allele_Count", "Allele_Number", "Allele_Frequency", "Homozygote_Count", "Hemizygote_Count", 
                    "Allele_Count_African_African_American", "Allele_Number_African_African_American", "Homozygote_Count_African_African_American", "Hemizygote_Count_African_African_American", 
                    "Allele_Count_Latino_Admixed_American", "Allele_Number_Latino_Admixed_American", "Homozygote_Count_Latino_Admixed_American", "Hemizygote_Count_Latino_Admixed_American", 
                    "Allele_Count_East_Asian", "Allele_Number_East_Asian", "Homozygote_Count_East_Asian", "Hemizygote_Count_East_Asian", 
                    "Allele_Count_European_Finnish", "Allele_Number_European_Finnish", "Homozygote_Count_European_Finnish", "Hemizygote_Count_European_Finnish", 
                    "Allele_Count_European_non_Finnish", "Allele_Number_European_non_Finnish", "Homozygote_Count_European_non_Finnish", "Hemizygote_Count_European_non_Finnish", 
                    "Allele_Count_Other", "Allele_Number_Other", "Homozygote_Count_Other", "Hemizygote_Count_Other", 
                    "Allele_Count_South_Asian", "Allele_Number_South_Asian", "Homozygote_Count_South_Asian", "Hemizygote_Count_South_Asian",
                    "Allele_Count_Ashkenazi_Jewish", "Allele_Number_Ashkenazi_Jewish", "Homozygote_Count_Ashkenazi_Jewish", "Hemizygote_Count_Ashkenazi_Jewish")
nans <- nans[c("Chromosome", "Position", "rsID", "Reference", "Alternate", "Source", "Filters_exomes", "Filters_genomes", "Transcript",
               "HGVS_Consequence", "Protein_Consequence", "Transcript_Consequence", "VEP_Annotation", "ClinVar_Clinical_Significance", "ClinVar_Variation_ID", "Flags", 
               "Allele_Count", "Allele_Number", "Allele_Frequency", "Homozygote_Count", "Hemizygote_Count", 
               "Allele_Count_African_African_American", "Allele_Number_African_African_American", "Homozygote_Count_African_African_American", "Hemizygote_Count_African_African_American", 
               "Allele_Count_Latino_Admixed_American", "Allele_Number_Latino_Admixed_American", "Homozygote_Count_Latino_Admixed_American", "Hemizygote_Count_Latino_Admixed_American", 
               "Allele_Count_Ashkenazi_Jewish", "Allele_Number_Ashkenazi_Jewish", "Homozygote_Count_Ashkenazi_Jewish", "Hemizygote_Count_Ashkenazi_Jewish", 
               "Allele_Count_East_Asian", "Allele_Number_East_Asian", "Homozygote_Count_East_Asian", "Hemizygote_Count_East_Asian", 
               "Allele_Count_European_Finnish", "Allele_Number_European_Finnish", "Homozygote_Count_European_Finnish", "Hemizygote_Count_European_Finnish", 
               "Allele_Count_European_non_Finnish", "Allele_Number_European_non_Finnish", "Homozygote_Count_European_non_Finnish", "Hemizygote_Count_European_non_Finnish", 
               "Allele_Count_Other", "Allele_Number_Other", "Homozygote_Count_Other", "Hemizygote_Count_Other", 
               "Allele_Count_South_Asian", "Allele_Number_South_Asian", "Homozygote_Count_South_Asian", "Hemizygote_Count_South_Asian")]
write.table(nans, file="NANS.ref.csv", sep=",", row.names = FALSE)


reffiles <- list.files(pattern='.ref.csv')
registerDoParallel(numCores)
foreach (reffile = reffiles) %dopar% {
  csv1 <- read.csv(reffile, header = TRUE, na.strings=c("","NA"))
  colnames(csv1) <- c("Chromosome", "Position", "rsID", "Reference", "Alternate", "Source", "Filters_exomes", "Filters_genomes", "Transcript",
                      "HGVS_Consequence", "Protein_Consequence", "Transcript_Consequence", "VEP_Annotation", "ClinVar_Clinical_Significance", "ClinVar_Variation_ID", "Flags", 
                      "Allele_Count", "Allele_Number", "Allele_Frequency", "Homozygote_Count", "Hemizygote_Count", 
                      "Allele_Count_African_African_American", "Allele_Number_African_African_American", "Homozygote_Count_African_African_American", "Hemizygote_Count_African_African_American", 
                      "Allele_Count_Latino_Admixed_American", "Allele_Number_Latino_Admixed_American", "Homozygote_Count_Latino_Admixed_American", "Hemizygote_Count_Latino_Admixed_American", 
                      "Allele_Count_Ashkenazi_Jewish", "Allele_Number_Ashkenazi_Jewish", "Homozygote_Count_Ashkenazi_Jewish", "Hemizygote_Count_Ashkenazi_Jewish", 
                      "Allele_Count_East_Asian", "Allele_Number_East_Asian", "Homozygote_Count_East_Asian", "Hemizygote_Count_East_Asian", 
                      "Allele_Count_European_Finnish", "Allele_Number_European_Finnish", "Homozygote_Count_European_Finnish", "Hemizygote_Count_European_Finnish", 
                      "Allele_Count_European_non_Finnish", "Allele_Number_European_non_Finnish", "Homozygote_Count_European_non_Finnish", "Hemizygote_Count_European_non_Finnish", 
                      "Allele_Count_Other", "Allele_Number_Other", "Homozygote_Count_Other", "Hemizygote_Count_Other", 
                      "Allele_Count_South_Asian", "Allele_Number_South_Asian", "Homozygote_Count_South_Asian", "Hemizygote_Count_South_Asian")

  csv1$Gene <- rep(reffile,nrow(csv1))
  csv1$Gene <- gsub('.ref.csv', '', csv1$Gene)
  
  csv1$Variant <- paste0(csv1$Chromosome, "-", csv1$Position, "-", csv1$Reference, "-", csv1$Alternate)
  csv1 <- mutate(csv1, Global_AF=Allele_Count/Allele_Number)
  csv1 <- mutate(csv1, African_African_American_AF=Allele_Count_African_African_American/Allele_Number_African_African_American)
  csv1 <- mutate(csv1, Latino_Admixed_American_AF=Allele_Count_Latino_Admixed_American/Allele_Number_Latino_Admixed_American)
  csv1 <- mutate(csv1, Ashkenazi_Jewish_AF=Allele_Count_Ashkenazi_Jewish/Allele_Number_Ashkenazi_Jewish)
  csv1 <- mutate(csv1, East_Asian_AF=Allele_Count_East_Asian/Allele_Number_East_Asian)
  csv1 <- mutate(csv1, European_Finnish_AF=Allele_Count_European_Finnish/Allele_Number_European_Finnish)
  csv1 <- mutate(csv1, European_Non_Finnish_AF=Allele_Count_European_non_Finnish/Allele_Number_European_non_Finnish)
  csv1 <- mutate(csv1, Other_AF=Allele_Count_Other/Allele_Number_Other)
  csv1 <- mutate(csv1, South_Asian_AF=Allele_Count_South_Asian/Allele_Number_South_Asian)
  
  csv1 <- add_columns(csv1, panel, by="Gene")
  
  csv1 <- csv1[c("Chromosome", "rsID", "Source", "Filters_exomes", "Filters_genomes", "Transcript", "HGVS_Consequence", "Protein_Consequence", 
                 "Transcript_Consequence", "VEP_Annotation", "ClinVar_Clinical_Significance", "ClinVar_Variation_ID", "Flags", "Gene", "Variant", 
                 "Global_AF", "African_African_American_AF", "Latino_Admixed_American_AF", "Ashkenazi_Jewish_AF", "East_Asian_AF", "European_Finnish_AF", 
                 "European_Non_Finnish_AF", "Other_AF", "South_Asian_AF", "Transcript_ID", "Protein_ID", "Protein_ID_for_PP2")]
  filename <- paste0(reffile, ".af.csv")
  write.table(csv1, file=filename, sep=",", row.names = FALSE)
  file.move(reffile, "./original_files/d_gene_ref", overwrite=TRUE)
}
Sys.sleep(2)


# split variant ref by chr
reffiles <- list.files(pattern='.af.csv')
reflib <- lapply(reffiles, read.csv)
reflib <- do.call(rbind.data.frame, reflib)
registerDoParallel(numCores)
foreach (reffile = reffiles) %dopar% {
  file.move(reffile, "./cache/c_gene_ref_trim", overwrite=TRUE)
}
chrs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
          "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
registerDoParallel(numCores)
foreach (chr = chrs) %dopar% {
  csv1 <- filter(reflib, Chromosome == chr)
  if(nrow(csv1)>0){
    csv1 <- csv1[c("Variant", "rsID", "Gene", "HGVS_Consequence", "VEP_Annotation", "Protein_ID", "Protein_Consequence", 
                   "Transcript_ID", "Transcript_Consequence", "Filters_exomes", "Filters_genomes", "Global_AF", 
                   "African_African_American_AF", "Latino_Admixed_American_AF", "Ashkenazi_Jewish_AF", "East_Asian_AF", "South_Asian_AF", 
                   "European_Finnish_AF", "European_Non_Finnish_AF", "Other_AF", "Source", "Protein_ID_for_PP2", 
                   "Flags", "ClinVar_Clinical_Significance","ClinVar_Variation_ID")]
    filename <- paste0("chr", chr, ".ref.csv")
    write.table(csv1, file=filename, sep=",", row.names = FALSE)
  }
}
rm("reflib", "chrs", "reffile", "reffiles", "ash_jew", "nans")
Sys.sleep(1)


# Decompress all .genome.vcf.gz files
gzfiles <- list.files(pattern=".genome.vcf.gz")
if (is_empty(gzfiles)==FALSE){
  print("Scanning for all .genome.vcf.gz files...")
  print(gzfiles)
  registerDoParallel(numCores)
  foreach (gzfile = gzfiles) %dopar% {
    gunzip(gzfile, remove=FALSE, overwrite=TRUE)
    file.move(gzfile, "./original_files/c_vcf_gz_files", overwrite=TRUE)
  }
  print("Decompressed all vcf.gz files")
}
Sys.sleep(1)
gzfiles <- list.files(pattern=".gz")
if (is_empty(gzfiles)==FALSE){
  registerDoParallel(numCores)
  foreach (gzfile = gzfiles) %dopar% {
    file.move(gzfile, "./original_files/c_vcf_gz_files/a_non_genome_vcf", overwrite=TRUE)
  }
}
rm("gzfiles")
Sys.sleep(2)


# Convert .vcf files to .csv files
vcffiles <- list.files(pattern=".genome.vcf")
samplesinthisbatch <- vcffiles
if (is_empty(vcffiles)==FALSE){
  print("Extracting variants info from VCF files...(this might take a while)")
  registerDoParallel(numCores)
  foreach (vcffile = vcffiles) %dopar% {
    allvariants <- read.vcfR(vcffile, verbose=FALSE)
    filename <- paste0(vcffile, ".allvariants.csv")
    write.table(cbind(allvariants@fix, allvariants@gt), 
                file=filename, sep=",", row.name=FALSE)
    file.move(vcffile, "./cache/d_unzipped_vcf", overwrite=TRUE)
  }
  print("All variants info is extracetd")
}
rm("vcffiles")
Sys.sleep(2)


# Translate "Genotype" part
rawcsvs <- list.files(pattern=".allvariants.csv")
if (is_empty(rawcsvs)==FALSE){
  print("Translating variants info...(This might take a while)")
  registerDoParallel(numCores)
  foreach (rawcsv = rawcsvs) %dopar% {
    csv1 <- read.csv(rawcsv, header = TRUE)
    if(nrow(csv1)>0){
      colnames(csv1) <- c("CHROM","POS","ID", "REF", "ALT", "QUAL", 
                          "FILTER", "INFO", "FORMAT", "Code")
      csv1 <- separate(csv1, "Code",
                       into=c("Genotype","Genotype_Quality", "Allele_Deapth", 
                              "Total_Depth", "Variant_Frequency", "NoiseLevel", 
                              "StrandBias", "Uncalled"), sep=":")
      filename <- paste0(rawcsv, ".translated.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
    }
    file.move(rawcsv, "./cache/e_raw_csv", overwrite=TRUE)
  }
  print("All variants info is translated")
}
rm("rawcsvs")
Sys.sleep(2)


# First filter by Qualigy, VAF, and Depth (paralleing)
print("Filtering by Quality, VAF, Depth, and Ref/Alt length...(this might take a while)")
allvariants <- list.files(pattern=".translated.csv")
registerDoParallel(numCores)
foreach (allvariant = allvariants) %dopar% {
  csv1 <- read.csv(allvariant)
  if(nrow(csv1)>0){
    csv2 <- filter(csv1, FILTER == "PASS" & Total_Depth >=50)
    csv3 <- filter(csv2, Variant_Frequency>=0.4 & Variant_Frequency<=0.6)
    csv4 <- filter(csv2, Variant_Frequency>=0.9)
    csv2 <- rbind(csv3, csv4)
    csv5 <- filter(csv1, Total_Depth >=10 & Total_Depth <50 )
    csv6 <- filter(csv5, Variant_Frequency>=0.3 & Variant_Frequency<=0.7)
    csv7 <- filter(csv5, Variant_Frequency>=0.8)
    csv5 <- rbind(csv6, csv7)
    csv1 <- rbind(csv2, csv5)
    csv1$REF <- as.character(csv1$REF)
    csv1$ALT <- as.character(csv1$ALT)
    csv1$Ref_length <- nchar(csv1$REF)
    csv1$Alt_length <- nchar(csv1$ALT)
    large_indel <- filter(csv1, Alt_length>=6 | Ref_length>=6)
    csv1 <- filter(csv1, Alt_length<6 & Ref_length<6)
    filename <- paste0(allvariant, ".filter.csv")
    write.table(csv1, file=filename, sep=",", row.names = FALSE)
    filename <- paste0(allvariant, ".largeindel.csv")
    write.table(large_indel, file=filename, sep=",", row.names = FALSE)
  }
  file.move(allvariant, "./cache/f_all_variants", overwrite=TRUE)
}
print("Filtering completed")
rm("allvariants")
Sys.sleep(2)


# Format variants info
print("Trimming data...")
medcsvs <- list.files(pattern=".csv.filter.csv")
registerDoParallel(numCores)
foreach (medcsv = medcsvs) %dopar% {
  csv1 <- read.csv(medcsv)
  if(nrow(csv1)>0){
    csv1$REF <- gsub("TRUE", "T", csv1$REF)
    csv1$ALT <- gsub("TRUE", "T", csv1$ALT)
    csv1$Sample <- rep(medcsv,nrow(csv1))
    csv1$Sample <- gsub('.genome.vcf.allvariants.csv.translated.csv.filter.csv', '', csv1$Sample)
    csv1$Variant <- paste0(csv1$CHROM, "-", csv1$POS, "-", csv1$REF, "-", csv1$ALT)
    csv1$Variant <- gsub('chr', '', csv1$Variant)
    csv1 <- csv1[c("Sample", "Variant", "ID", "QUAL", "FILTER", "Genotype", "Genotype_Quality", 
                   "Allele_Deapth", "Total_Depth", "Variant_Frequency", "NoiseLevel", "StrandBias", 
                   "Uncalled", "INFO", "CHROM")]
    filename <- paste0(medcsv, ".trim.csv")
    write.table(csv1, file=filename, sep=",", row.names = FALSE)
  }
  file.move(medcsv, "./cache/g_filtered_variants", overwrite=TRUE)
}

medcsvs <- list.files(pattern=".csv.largeindel.csv")
registerDoParallel(numCores)
foreach (medcsv = medcsvs) %dopar% {
  csv1 <- read.csv(medcsv)
  csv1$REF <- gsub("TRUE", "T", csv1$REF)
  csv1$ALT <- gsub("TRUE", "T", csv1$ALT)
  if(nrow(csv1)>0){
    csv1$Sample <- rep(medcsv,nrow(csv1))
    csv1$Sample <- gsub('.genome.vcf.allvariants.csv.translated.csv.largeindel.csv', '', csv1$Sample)
    csv1$Variant <- paste0(csv1$CHROM, "-", csv1$POS, "-", csv1$REF, "-", csv1$ALT)
    csv1$Variant <- gsub('chr', '', csv1$Variant)
    csv1 <- csv1[c("Sample", "Variant", "ID", "QUAL", "FILTER", "Genotype", "Genotype_Quality", 
                   "Allele_Deapth", "Total_Depth", "Variant_Frequency", "NoiseLevel", "StrandBias", 
                   "Uncalled", "INFO")]
    filename <- gsub('.genome.vcf.allvariants.csv.translated.csv.largeindel.csv', '', medcsv)
    filename <- paste0(filename, ".indeltrim.csv")
    write.table(csv1, file=filename, sep=",", row.names = FALSE)
  }
  file.move(medcsv, "./cache/h_all_large_indels", overwrite=TRUE)
}
large_indels <- list.files(pattern=".indeltrim.csv")
large_indel <- lapply(large_indels, read.csv)
large_indel <- do.call(rbind.data.frame, large_indel)
file.move(large_indels, "./cache/i_trimmed_large_indels", overwrite=TRUE)
print("Trimming completed")
rm("medcsvs")
Sys.sleep(2)


# Split variants by chroms
print("Preparing varinats info for GnomAD cross searching...")
medwelcsvs <- list.files(pattern=".trim.csv")
csv1 <- lapply(medwelcsvs, read.csv)
csv1 <- do.call(rbind.data.frame, csv1)
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
          "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
          "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
registerDoParallel(numCores)
foreach (chr = chrs) %dopar% {
  csv2 <- filter(csv1, CHROM == chr)
    if(nrow(csv2)>0){
      csv2 <- separate(csv2,"Sample",into = c("Sample","Suffix"),sep = "_",remove = FALSE,extra = "merge")
      csv2 <- csv2[c("Sample", "Variant", "ID", "QUAL", "FILTER", "Genotype", "Genotype_Quality", 
                     "Allele_Deapth", "Total_Depth", "Variant_Frequency", "NoiseLevel", "StrandBias", 
                     "Uncalled", "INFO")]
      filename <- paste0(chr, ".seq.csv")
      write.table(csv2, file=filename, sep=",", row.names = FALSE)
    }
}
registerDoParallel(numCores)
foreach (medwelcsv = medwelcsvs) %dopar% {
  file.move(medwelcsv, "./cache/j_trimmed_variants", overwrite=TRUE)
}
print("Variants info ready")
rm("medwelcsvs", "csv1", "chrs")
Sys.sleep(1)


# Align and filter with AF and consequences
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
          "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
          "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
registerDoParallel(numCores)
foreach (chr = chrs) %dopar% {
  filename1 <- paste0(chr, ".seq.csv")
  filename2 <- paste0(chr, ".ref.csv")
  filename3 <- paste0(chr, ".align.csv")
  if (file.exists(filename1)){
    csv1 <- read.csv(filename1)
    csv2 <- read.csv(filename2)
    csv2 <- csv2[!duplicated(csv2[,"Variant"]),]
    csv3 <- add_columns(csv1, csv2, by="Variant")
    
    csv3$Global_AF <- replace(csv3$Global_AF, is.na(csv3$Global_AF), 0)
    csv3$African_African_American_AF <- replace(csv3$African_African_American_AF, is.na(csv3$African_African_American_AF), 0)
    csv3$Latino_Admixed_American_AF <- replace(csv3$Latino_Admixed_American_AF, is.na(csv3$Latino_Admixed_American_AF), 0)
    csv3$Ashkenazi_Jewish_AF <- replace(csv3$Ashkenazi_Jewish_AF, is.na(csv3$Ashkenazi_Jewish_AF), 0)
    csv3$East_Asian_AF <- replace(csv3$East_Asian_AF, is.na(csv3$East_Asian_AF), 0)
    csv3$South_Asian_AF <- replace(csv3$South_Asian_AF, is.na(csv3$South_Asian_AF), 0)
    csv3$European_Finnish_AF <- replace(csv3$European_Finnish_AF, is.na(csv3$European_Finnish_AF), 0)
    csv3$European_Non_Finnish_AF <- replace(csv3$European_Non_Finnish_AF, is.na(csv3$European_Non_Finnish_AF), 0)
    csv3$Other_AF <- replace(csv3$Other_AF, is.na(csv3$Other_AF), 0)
    csv3$VEP_Annotation <- replace(csv3$VEP_Annotation, is.na(csv3$VEP_Annotation), "NA")
    
    csv3 <- filter(csv3, Global_AF<0.01 & African_African_American_AF<0.01 & Latino_Admixed_American_AF<0.01 & 
                     Ashkenazi_Jewish_AF<0.01 & East_Asian_AF<0.01 & South_Asian_AF<0.01 & 
                     European_Finnish_AF<0.01 & European_Non_Finnish_AF<0.01 & Other_AF<0.01)
    csv3 <- filter(csv3, !grepl('synonymous_variant', INFO, fixed=TRUE) &
                     !grepl('5_prime_UTR_variant', INFO, fixed=TRUE) &
                     !grepl('3_prime_UTR_variant', INFO, fixed=TRUE))
    csv3 <- filter(csv3, VEP_Annotation != "5_prime_UTR_variant" & VEP_Annotation != "synonymous_variant" & 
                     VEP_Annotation !="3_prime_UTR_variant" & VEP_Annotation != "stop_retained_variant")
    
    # if INFO and VEP both annotated, remove only when both contain "intron_variant"
    csv3 <- filter(csv3, !(grepl('intron_variant', INFO, fixed=TRUE) & VEP_Annotation == "intron_variant"))
    
    # if only INFO annotated and VEP is NA, remove by INFO
    csv3 <- filter(csv3, !(grepl('splice_region_variant', INFO, fixed=TRUE) & grepl('intron_variant', INFO, fixed=TRUE)))
    csv4 <- filter(csv3, grepl('regulatory_region_variant', INFO, fixed=TRUE) | grepl('upstream_gene_variant', INFO, fixed=TRUE) |
                     grepl('downstream_gene_variant', INFO, fixed=TRUE))
    csv4 <- filter(csv4, grepl('missense_variant', INFO, fixed=TRUE))
    csv3 <- filter(csv3, !(grepl('regulatory_region_variant', INFO, fixed=TRUE) | grepl('upstream_gene_variant', INFO, fixed=TRUE) |
                     grepl('downstream_gene_variant', INFO, fixed=TRUE)))
    csv3 <- rbind(csv3, csv4)
    csv3 <- filter(csv3, !(grepl('intron_variant', INFO, fixed=TRUE) & VEP_Annotation == "NA"))
    
        #Retained:
        #NA
        #missense_variant
        #stop_gained
        #inframe_deletion
        #splice_region_variant
        #splice_donor_variant
        #frameshift_variant
        #splice_acceptor_variant
        #stop_lost
        #start_lost
        #inframe_insertion
        #protein_altering_variant
    if(nrow(csv3)>0){
      write.table(csv3, file=filename3, sep=",", row.names = FALSE)
    }
    file.move(filename1, "./cache/k_chr_variants", overwrite=TRUE)
    file.move(filename2, "./cache/l_chr_refs", overwrite=TRUE)
  }
}
rm("chrs")
Sys.sleep(10)


# SNP check
print("Retriving most updated SNP info from NCBI...(This might take a while")
aligns <- list.files(pattern=".align.csv")
for (align in aligns) {
  rsIDs <- read.csv(align)[, "rsID"]
  rsIDs <- rsIDs[!is.na(rsIDs)]
  current_ver <- sapply(binman::list_versions("chromedriver"),"[[",1)
  cDrv <- chrome(version = current_ver, verbose = FALSE)
  eCaps <- list(
    chromeOptions = list(
      args = c('--headless', '--disable-gpu', '--window-size=1280,800'), 
      prefs = list(
        "profile.default_content_settings.popups" = 0L,
        "download.prompt_for_download" = FALSE,
        "download.default_directory" = "."
      )
    )
  )
  for (rsID in rsIDs) {
    remDr <- remoteDriver(browserName = "chrome", port = 4567L, extraCapabilities = eCaps)
    remDr$open()
    Sys.sleep(1)
    remDr$setTimeout(type = "implicit", milliseconds = 10000)
    remDr$setTimeout(type = "page load", milliseconds = 10000)
    new <- rsID
    snpurl <- paste0("https://www.ncbi.nlm.nih.gov/snp/", rsID)
    # snpurl<- "https://www.ncbi.nlm.nih.gov/snp/rs754435604"
    remDr$navigate (snpurl)
    Sys.sleep(2)
    txt<-remDr$findElement(using='css selector',"body")$getElementText()
    if (str_detect(txt, "was merged into") == TRUE){
      elem <- remDr$findElement(using="xpath", "//div[@id='main_content']/main/div[4]/div[2]/dd/a")
      # Sys.sleep(2)
      new <- elem$getElementText()[[1]]
      snpurl <- paste0("https://www.ncbi.nlm.nih.gov/snp/", new)
      remDr$navigate(snpurl)
      Sys.sleep(2)
    }
    elem <-remDr$findElement(using="xpath", '//*[(@id = "DataTables_Table_0")]')
    elemtxt <- elem$getElementAttribute("outerHTML")[[1]] # gets us the HTML
    elemxml <- htmlTreeParse(elemtxt, useInternalNodes=T) # parse string into HTML tree to allow for querying with XPath
    table0 <- readHTMLTable(elemxml)
    table <- as.data.frame(table0)
    colnames(table) <- c("Molecule_type", "Change", "Amino_acid_Codon", "SO_Term")
    error="Error : \t Summary: NoSuchElement\n \t Detail: An element could not be located on the page using the given search parameters.\n\t Further Details: run errorDetails method\n"
    table1 <- try(unlist(remDr$findElement(using='xpath', '//*[(@id = "DataTables_Table_1")]')$getElementAttribute('id')), silent = TRUE) # silent = TRUE
    if(table1[1]!=error){
      elem <-remDr$findElement(using="xpath", '//*[(@id = "DataTables_Table_1")]')
      elemtxt <- elem$getElementAttribute("outerHTML")[[1]] # gets us the HTML
      elemxml <- htmlTreeParse(elemtxt, useInternalNodes=T) # parse string into HTML tree to allow for querying with XPath
      table1 <- readHTMLTable(elemxml)
      table1 <- as.data.frame(table1)
      colnames(table1) <- c("Molecule_type", "Change", "Amino_acid_Codon", "SO_Term")
      table <- rbind(table, table1)
    }
    filename <- paste0(rsID, ">", new, ".snp.csv")
    write.table(table, file=filename, sep=",", row.name=F)
    remDr$close()
    Sys.sleep(10)
  }
  cDrv$stop()
}
csv1 <- lapply(aligns, read.csv)
csv1 <- do.call(rbind.data.frame, csv1)
write.table(csv1, file="all.filtered.csv", sep=",", row.names = FALSE)
file.move(aligns, "./cache/m_variants_aligns", overwrite=TRUE)
rm("cDrv", "eCaps", "table", "table0", "table1", "tables", "elem", "elemtxt", "elemxml", "filename", 
   "snpurl", "remDr", "rsID", "rsIDs", "csv1", "txt", "new", "align", "aligns", "error", "current_ver")
print("Most updated SNP info retrived")
Sys.sleep(10)


# Trim NCBI SNP info
print("Trimming SNP info from NCBI...")
ncbisnps <- list.files(pattern=".snp.csv")
registerDoParallel(numCores)
foreach (ncbisnp = ncbisnps) %dopar% {
  csv1 <- read.csv(ncbisnp)
  csv1$rsID <- rep(ncbisnp,nrow(csv1))
  csv1$rsID <- gsub('.snp.csv', '', csv1$rsID)
  csv1 <- separate(csv1, "rsID",
                   into=c("rsID","NCBI_rsID"), sep=">")
  csv1 <- separate(csv1, "Change",
                   into=c("RefSeq","Change"), sep=":")
  filename <- paste0(ncbisnp, ".ncbisnp.csv")
  write.table(csv1, file=filename, sep=",", row.name=F)
  file.move(ncbisnp, "./original_files/e_ncbi_snp", overwrite=TRUE)
}
Sys.sleep(2)

ncbisnps <- list.files(pattern=".ncbisnp.csv")
csv1 <- lapply(ncbisnps, read.csv)
csv1 <- do.call(rbind.data.frame, csv1)
registerDoParallel(numCores)
foreach (ncbisnp = ncbisnps) %dopar% {
  file.move(ncbisnp, "./cache/n_ncbi_snp", overwrite=TRUE)
}
csv1 <- csv1[c("rsID", "RefSeq", "Change", "SO_Term", "NCBI_rsID")]
names(csv1) <- gsub("SO_Term", "NCBI_consequence", names(csv1))
csv2 <- csv1[c("rsID", "NCBI_rsID")]
write.table(csv1, file="all.refsnp.csv", sep=",", row.name=F)
write.table(csv2, file="rsIDs.csv", sep=",", row.name=F)
rm("csv1", "csv2", "ncbisnps")
Sys.sleep(2)


# rsID confirmation, sequenced ins/del/fs isolation, refsnp ins/del/fs isolation
csv1 <- read.csv("all.filtered.csv", header = TRUE, na.strings=c("","NA"))
unannotated <- filter(csv1, is.na(HGVS_Consequence))
csv2 <- read.csv("rsIDs.csv", header = TRUE, na.strings=c("","NA"))
csv2 <- unique(csv2)
csv1 <- add_columns(csv1, csv2, by=c("rsID"))
csv1$dbSNP <- ifelse(csv1$rsID == csv1$NCBI_rsID, paste0(csv1$rsID), paste0(csv1$NCBI_rsID))
csv1$dbSNP_confirmed_with_NCBI <- ifelse(csv1$rsID == csv1$NCBI_rsID, "Yes", "No")
csv1 <- csv1[c("Sample", "Variant", "QUAL", "FILTER", "Genotype", "Genotype_Quality", "Allele_Deapth", 
               "Total_Depth", "Variant_Frequency", "NoiseLevel", "StrandBias", "Uncalled", "INFO", "dbSNP", 
               "dbSNP_confirmed_with_NCBI", "Gene", "HGVS_Consequence", "VEP_Annotation", "Protein_ID", "Protein_Consequence", 
               "Transcript_ID", "Transcript_Consequence", "Filters_exomes", "Filters_genomes", "Global_AF", "African_African_American_AF", 
               "Latino_Admixed_American_AF", "Ashkenazi_Jewish_AF", "East_Asian_AF", "South_Asian_AF", "European_Finnish_AF", 
               "European_Non_Finnish_AF", "Other_AF", "Source", "Protein_ID_for_PP2", "Flags", "ClinVar_Clinical_Significance", "ClinVar_Variation_ID")]
csv2 <- csv1 %>% 
  filter(str_detect(HGVS_Consequence, "ins") | str_detect(HGVS_Consequence, "del") | str_detect(HGVS_Consequence, "fs"))
csv1 <- csv1[ grep("ins", csv1$HGVS_Consequence, invert = TRUE), ]
csv1 <- csv1[ grep("del", csv1$HGVS_Consequence, invert = TRUE), ]
csv1 <- csv1[ grep("fs", csv1$HGVS_Consequence, invert = TRUE), ]
write.table(csv1, file="nidfs.filtered.csv", sep=",", row.names = FALSE)
write.table(csv2, file="idfs.filtered.csv", sep=",", row.names = FALSE)
file.move("all.filtered.csv", "./cache/o_merge", overwrite=TRUE)
file.move("rsIDs.csv", "./cache/o_merge", overwrite=TRUE)
Sys.sleep(2)

csv1 <- read.csv("all.refsnp.csv")
csv2 <- csv1 %>% 
  filter(str_detect(Change, "ins") | str_detect(Change, "del") | str_detect(Change, "fs"))
csv1 <- csv1[ grep("ins", csv1$Change, invert = TRUE), ]
csv1 <- csv1[ grep("del", csv1$Change, invert = TRUE), ]
csv1 <- csv1[ grep("fs", csv1$Change, invert = TRUE), ]
write.table(csv1, file="nidfs.refsnp.csv", sep=",", row.names = FALSE)
write.table(csv2, file="idfs.refsnp.csv", sep=",", row.names = FALSE)
file.move("all.refsnp.csv", "./cache/o_merge", overwrite=TRUE)
rm("csv1", "csv2")
Sys.sleep(2)


# Organize substitution variants
csv1 <- read.csv("nidfs.filtered.csv")
csv2 <- csv1 %>% 
  filter(str_detect(HGVS_Consequence, "^p."))
csv3 <- csv1 %>%
  filter(str_detect(HGVS_Consequence, "^c."))
csv4 <- read.csv("nidfs.refsnp.csv")
csv4$dbSNP <- ifelse(csv4$rsID == csv4$NCBI_rsID, paste0(csv4$rsID), paste0(csv4$NCBI_rsID))
csv4 <- csv4[c("dbSNP", "RefSeq", "Change", "NCBI_consequence")]
csv5 <- csv4 %>% 
  filter(str_detect(RefSeq, "^NP_"))
csv6 <- csv4 %>%
  filter(str_detect(RefSeq, "^NM_"))
Sys.sleep(1)
# Protein
aachange <- as.matrix(csv2$Protein_Consequence)
aachange <- gsub('[0-9]+', '', aachange)
aachange <- as.data.frame(aachange)
csv2 <- cbind(csv2, aachange)
names(csv2) <- gsub("Protein_ID", "RefSeq", names(csv2))
aachange <- as.matrix(csv5$Change)
aachange <- gsub('[0-9]+', '', aachange)
aachange <- as.data.frame(aachange)
csv5 <- cbind(csv5, aachange)
csv5 <- unique(csv5)
csvp <- add_columns(csv2, csv5, by=c("dbSNP", "RefSeq", "V1"))
names(csvp) <- gsub("Change", "NCBI_Consequence", names(csvp))
names(csvp) <- gsub("NCBI_consequence", "NCBI_Annotation", names(csvp))
names(csvp) <- gsub("RefSeq", "Protein_ID", names(csvp))
csvp <- subset(csvp, select = -c(V1))
Sys.sleep(1)
# Transcript
ntchange <- as.matrix(csv3$Transcript_Consequence)
ntchange <- gsub('[0-9]+', '', ntchange)
ntchange <- as.data.frame(ntchange)
csv3 <- cbind(csv3, ntchange)
names(csv3) <- gsub("Transcript_ID", "RefSeq", names(csv3))
ntchange <- as.matrix(csv6$Change)
ntchange <- gsub('[0-9]+', '', ntchange)
ntchange <- as.data.frame(ntchange)
csv6 <- cbind(csv6, ntchange)
csv6 <- unique(csv6)
csvt <- add_columns(csv3, csv6, by=c("dbSNP", "RefSeq", "V1"))
names(csvt) <- gsub("Change", "NCBI_Consequence", names(csvt))
names(csvt) <- gsub("NCBI_consequence", "NCBI_Annotation", names(csvt))
names(csvt) <- gsub("RefSeq", "Transcript_ID", names(csvt))
csvt <- subset(csvt, select = -c(V1))
Sys.sleep(1)
# export
write.table(csvp, file="nidfs.protein.csv", sep=",", row.names = FALSE)
write.table(csvt, file="nidfs.transcript.csv", sep=",", row.names = FALSE)
file.move("nidfs.filtered.csv", "./cache/p_pro_trans", overwrite=TRUE)
file.move("nidfs.refsnp.csv", "./cache/p_pro_trans", overwrite=TRUE)
txt <- readLines("nidfs.protein.csv")
txt <- gsub(',NA', ',"Not_Found"', txt)
writeLines(txt, "nidfs.protein.csv")
txt <- readLines("nidfs.transcript.csv")
txt <- gsub(',NA', ',"Not_Found"', txt)
writeLines(txt, "nidfs.transcript.csv")
rm("csv1", "csv2", "csv3", "csv4", "csv5", "csv6", "csvp", "csvt", "aachange", "ntchange", "txt")
Sys.sleep(1)


# Organize insertion/deletion variants
csv1 <- read.csv("idfs.filtered.csv")
csv2 <- csv1 %>% 
  filter(str_detect(HGVS_Consequence, "^p."))
csv3 <- csv1 %>%
  filter(str_detect(HGVS_Consequence, "^c."))
csv4 <- read.csv("idfs.refsnp.csv")
csv4$dbSNP <- ifelse(csv4$rsID == csv4$NCBI_rsID, paste0(csv4$rsID), paste0(csv4$NCBI_rsID))
csv4 <- csv4[c("dbSNP", "RefSeq", "Change", "NCBI_consequence")]
csv5 <- csv4 %>% 
  filter(str_detect(RefSeq, "^NP_"))
csv6 <- csv4 %>%
  filter(str_detect(RefSeq, "^NM_"))
Sys.sleep(1)
# Protein
names(csv2) <- gsub("Protein_ID", "RefSeq", names(csv2))
csv5 <- unique(csv5)
csvp <- add_columns(csv2, csv5, by=c("dbSNP", "RefSeq"))
names(csvp) <- gsub("Change", "NCBI_Consequence", names(csvp))
names(csvp) <- gsub("NCBI_consequence", "NCBI_Annotation", names(csvp))
names(csvp) <- gsub("RefSeq", "Protein_ID", names(csvp))
Sys.sleep(1)
# Transcript
names(csv3) <- gsub("Transcript_ID", "RefSeq", names(csv3))
csv6 <- unique(csv6)
csvt <- add_columns(csv3, csv6, by=c("dbSNP", "RefSeq"))
names(csvt) <- gsub("Change", "NCBI_Consequence", names(csvt))
names(csvt) <- gsub("NCBI_consequence", "NCBI_Annotation", names(csvt))
names(csvt) <- gsub("RefSeq", "Transcript_ID", names(csvt))
Sys.sleep(1)
# export
write.table(csvp, file="idfs.protein.csv", sep=",", row.names = FALSE)
write.table(csvt, file="idfs.transcript.csv", sep=",", row.names = FALSE)
file.move("idfs.filtered.csv", "./cache/p_pro_trans", overwrite=TRUE)
file.move("idfs.refsnp.csv", "./cache/p_pro_trans", overwrite=TRUE)
txt <- readLines("idfs.protein.csv")
txt <- gsub(',NA', ',"Not_Found"', txt)
writeLines(txt, "idfs.protein.csv")
txt <- readLines("idfs.transcript.csv")
txt <- gsub(',NA', ',"Not_Found"', txt)
writeLines(txt, "idfs.transcript.csv")
rm("csv1", "csv2", "csv3", "csv4", "csv5", "csv6", "csvp", "csvt", "txt")
Sys.sleep(1)
print("SNP info trimmed")


# Prepare for prediction
print("Formatting variants for prediction...")
provar <- list.files(pattern=".protein.csv")
csv1 <- lapply(provar, read.csv)
csv1 <- do.call(rbind.data.frame, csv1)
travar <- list.files(pattern=".transcript.csv")
csv2 <- lapply(travar, read.csv)
csv2 <- do.call(rbind.data.frame, csv2)
csv3 <- csv1 %>% 
  filter(str_detect(HGVS_Consequence, "ins") | str_detect(HGVS_Consequence, "del") | 
           str_detect(HGVS_Consequence, "fs") | str_detect(HGVS_Consequence, "Ter"))
csv2 <- rbind(csv2, csv3)
csv1 <- csv1[ grep("ins", csv1$HGVS_Consequence, invert = TRUE), ]
csv1 <- csv1[ grep("del", csv1$HGVS_Consequence, invert = TRUE), ]
csv1 <- csv1[ grep("fs", csv1$HGVS_Consequence, invert = TRUE), ]
csv1 <- csv1[ grep("Ter", csv1$HGVS_Consequence, invert = TRUE), ]
write.table(csv1, file="protein.pre.csv", sep=",", row.names = FALSE)
write.table(csv2, file="transcript.pre.csv", sep=",", row.names = FALSE)
file.move(provar, "./cache/p_pro_trans", overwrite=TRUE)
file.move(travar, "./cache/p_pro_trans", overwrite=TRUE)
rm("csv1", "csv2", "csv3", "provar", "travar")
Sys.sleep(1)
# Protein prediction format
csv1 <- read.csv("protein.pre.csv")
csv2 <- csv1[c("Protein_ID", "Protein_ID_for_PP2", "NCBI_Consequence")]
csv2$NCBI_Consequence <- gsub('p.', '', csv2$NCBI_Consequence)
csv2$aapos <- stri_extract_all_regex(csv2$NCBI_Consequence, "[0-9]+")
csv2$NCBI_Consequence <- str_replace_all(csv2$NCBI_Consequence, c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", 
                                                                  "Asx"="B", "Cys"="C", "Glu"="E", "Gln"="Q", 
                                                                  "Glx"="Z", "Gly"="G", "His"="H", "Ile"="I", 
                                                                  "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", 
                                                                  "Pro"="P", "Ser"="S", "Thr"="T", "Trp"="W", 
                                                                  "Tyr"="Y", "Val"="V"))
csv2$NCBI_Consequence <- gsub('[0-9]+', '-', csv2$NCBI_Consequence)
csv2 <- separate(csv2, "NCBI_Consequence", into=c("aaref", "aaalt"), sep="-")
csv2$Provean_format <- paste(csv2$Protein_ID, csv2$aapos, csv2$aaref, csv2$aaalt)
csv2$PolyPhen_format <- paste(csv2$Protein_ID_for_PP2, csv2$aapos, csv2$aaref, csv2$aaalt)
csv3 <- csv2[c("Provean_format", "PolyPhen_format")]
csv1 <- cbind(csv1, csv3)
csv1$Provean_format <- ifelse(str_detect(csv1$Provean_format, "Not_Found"), paste0("Not_Found"), paste0(csv1$Provean_format))
csv1$PolyPhen_format <- ifelse(str_detect(csv1$PolyPhen_format, "Not_Found"), paste0("Not_Found"), paste0(csv1$PolyPhen_format))
write.table(csv1, file="protein.format.csv", sep=",", row.names = FALSE)
Sys.sleep(1)
# Transcript prediction format
csv1 <- read.csv("transcript.pre.csv")
csv2 <- read.csv("transcript.pre.csv")[, "Variant"]
csv2 <- as.data.frame(csv2)
csv2 <- separate(csv2, "csv2", into=c("chr", "pos", "ref", "alt"), sep="-")
csv2$Provean_format <- paste0(csv2$chr, ",", csv2$pos, ",", csv2$ref, ",", csv2$alt)
csv2$PolyPhen_format <- paste0("chr", csv2$chr, ":", csv2$pos, " ", csv2$ref, "/", csv2$alt)
csv3 <- csv2[c("Provean_format", "PolyPhen_format")]
csv1 <- cbind(csv1, csv3)
write.table(csv1, file="transcript.format.csv", sep=",", row.names = FALSE)
file.move("protein.pre.csv", "./cache/p_pro_trans", overwrite=TRUE)
file.move("transcript.pre.csv", "./cache/p_pro_trans", overwrite=TRUE)
rm("csv1", "csv2", "csv3")
Sys.sleep(1)


# Prepare for prediction
pro <- read.csv("protein.format.csv")
pro_polyphen <- pro[c("PolyPhen_format")]
pro_polyphen <- unique(pro_polyphen)
pro <- filter(pro, Provean_format != "Not_Found" & PolyPhen_format != "Not_Found")
pro_provean <- pro[c("Provean_format")]
pro_provean <- unique(pro_provean)
pro_provean <- as.matrix(pro_provean)
writeLines(pro_provean, "protein.provean.txt")
tra <- read.csv("transcript.format.csv")
tra_polyphen <- tra[c("PolyPhen_format")]
tra_polyphen <- unique(tra_polyphen)
tra <- filter(tra, Provean_format != "Not_Found" & PolyPhen_format != "Not_Found")
tra_provean <- tra[c("Provean_format")]
tra_provean <- unique(tra_provean)
tra_provean <- as.matrix(tra_provean)
writeLines(tra_provean, "transcript.provean.txt")
polyphen <- rbind(pro_polyphen, tra_polyphen)
polyphen<- unique(polyphen)
polyphen <- as.data.frame(polyphen)
polyphen <- as.matrix(filter(polyphen, !(grepl("Not_Found", PolyPhen_format, fixed=TRUE))))
writeLines(polyphen, "all.polyphen.txt")
rm("pro", "pro_provean", "pro_polyphen", "tra", "tra_provean", "tra_polyphen", "polyphen")


# Prediction
print("Performing variant prediction...")
current_ver <- sapply(binman::list_versions("chromedriver"),"[[",1)
cDrv <- chrome(version = current_ver, verbose = FALSE)
eCaps <- list(
  chromeOptions = list(
    args = c('--headless', '--disable-gpu', '--window-size=1280,800'),
    prefs = list(
      "profile.default_content_settings.popups" = 0L,
      "download.prompt_for_download" = FALSE,
      "download.default_directory" = "."
    )
  )
)
remDr <- remoteDriver(browserName = "chrome", port = 4567L, extraCapabilities = eCaps)
remDr$open()
Sys.sleep(1)
remDr$setTimeout(type = "implicit", milliseconds = 60000)
remDr$setTimeout(type = "page load", milliseconds = 60000)
# Provean protein prediction
proveanpro <- read_file("protein.provean.txt")
proveanpro <- gsub("\n", "\uE007", proveanpro)
remDr$navigate ("http://provean.jcvi.org/protein_batch_submit.php?species=human")
Sys.sleep(10)
textarea <-remDr$findElement(using="xpath", "//textarea[@id='variants']")
textarea$sendKeysToElement(list(proveanpro))
Sys.sleep(2)
submit <-remDr$findElement(using="xpath", "//input[@type='submit']")
submit$clickElement()
Sys.sleep(60)
download <-remDr$findElement(using="xpath", "(//a[contains(text(),'Download')])[2]")
download$clickElement()
Sys.sleep(2)
result <- list.files(pattern=".result.tsv")
file.rename(result, "protein.result.provean.tsv")
file.move("protein.provean.txt", "./cache/q_pre_format", overwrite=TRUE)
Sys.sleep(2)
# Provean transcription prediction
proveantra <- read_file("transcript.provean.txt")
proveantra <- gsub("\n", "\uE007", proveantra)
remDr$navigate ("http://provean.jcvi.org/genome_submit_2.php?species=human")
Sys.sleep(10)
textarea <-remDr$findElement(using="xpath", "//textarea[@id='CHR']")
textarea$sendKeysToElement(list(proveantra))
Sys.sleep(2)
submit <-remDr$findElement(using="xpath", "//input[@type='submit']")
submit$clickElement()
Sys.sleep(60)
download <-remDr$findElement(using="xpath", "(//a[contains(text(),'Download')])[3]")
download$clickElement()
Sys.sleep(2)
result <- list.files(pattern=".result.one.tsv")
file.rename(result, "transcript.result.provean.tsv")
file.move("transcript.provean.txt", "./cache/q_pre_format", overwrite=TRUE)
Sys.sleep(2)
# PolyPhen prediction
polyphen <- read_file("all.polyphen.txt")
polyphen <- gsub("\n", "\uE007", polyphen)
remDr$navigate ("http://genetics.bwh.harvard.edu/pph2/bgi.shtml")
Sys.sleep(10)
textarea <-remDr$findElement(using="xpath", "//textarea[@name='_ggi_batch']")
textarea$sendKeysToElement(list(polyphen))
Sys.sleep(2)
submit <-remDr$findElement(using="xpath", "//input[@name='_ggi_target_pipeline']")
submit$clickElement()
Sys.sleep(60)
refresh <-remDr$findElement(using="xpath", "//input[@name='_ggi_target_manage']")
refresh$clickElement()
Sys.sleep(2)
download <-remDr$findElement(using = "link text", "Short")
download$clickElement()
Sys.sleep(2)
pp2short <-remDr$findElement(using='css selector',"body")$getElementText()
write.table(pp2short, file="pp2short.txt", sep="\t", row.name=F)
Sys.sleep(2)
file.move("all.polyphen.txt", "./cache/q_pre_format", overwrite=TRUE)
remDr$close()
cDrv$stop()
# Export
pp2 <- readLines("pp2short.txt")
pp2 <- as.data.frame(pp2)
row <- nrow(pp2)
row_4 <- row-4
pp2 <- pp2[-c(1:2,row_4:row),]
write.table(pp2, file="pp2.csv", sep=",", row.name=F)
pp2 <- read.csv("pp2.csv", col.name="PP2")
pp2$PP2 <- str_replace_all(pp2$PP2, c("probably damaging"="probably_damaging", 
                                      "possibly damaging"="possibly_damaging"))
pp2$PP2 <- gsub('[ ]+', ',', pp2$PP2)
pp2 <- separate(pp2, "PP2", into=c("o_acc","o_pos","o_aa1","o_aa2","rsid","acc","pos","aa1","aa2",
                                   "prediction","pph2_prob","pph2_FPR","pph2_TPR"), sep=",")
write.table(pp2, file="all.result.polyphen.csv", sep=",", row.names = FALSE)
file.move("pp2short.txt", "./original_files/f_predictions", overwrite=TRUE)
file.move("pp2.csv", "./cache/r_predictions", overwrite=TRUE)
print("All prediction completed")
rm("cDrv", "eCaps", "download", "refresh", "submit", "proveanpro", "proveantra", "polyphen", "remDr", "result", 
   "textarea", "row", "row_4", "pp2short", "pp2", "current_ver")
Sys.sleep(2)


# Merge prediction results (unique, add_columns)
print("Final trimming...")
csv1 <- read.table("protein.result.provean.tsv", sep="\t", na.strings=c("","NA"))
colnames(csv1) <- c("#ROW_NO.", "Provean_format", "PROTEIN_ID", "POSITION", "RESIDUE_REF", "RESIDUE_ALT", "SCORE", 
                    "Provean_prediction", "#SEQ", "#CLUSTER", "SCORE", "SIFT_prediction", 
                    "MEDIAN_INFO", "#SEQ")
csv1 <- csv1[c("Provean_format", "Provean_prediction", "SIFT_prediction")]
csv2 <- read.table("transcript.result.provean.tsv",sep="\t", na.strings=c("","NA"))
colnames(csv2) <- c("#ROW_NO.", "Provean_format", "PROTEIN_ID", "LENGTH", "STRAND", "CODON_CHANGE", "POS", "RESIDUE_REF", 
                    "RESIDUE_ALT", "TYPE", "SCORE", "Provean_prediction", "#SEQ", "#CLUSTER", "SCORE", 
                    "SIFT_prediction", "MEDIAN_INFO", "#SEQ", "dbSNP_ID")
csv2 <- csv2[c("Provean_format", "Provean_prediction", "SIFT_prediction")]
csv3 <- rbind(csv1,csv2)
csv4 <- read.csv("all.result.polyphen.csv", header = TRUE, na.strings=c("","NA"))
csv4$PolyPhen_format <- paste(csv4$o_acc, csv4$o_pos, csv4$o_aa1, csv4$o_aa2)
names(csv4) <- gsub("prediction", "PolyPhen_prediction", names(csv4))
csv4 <- csv4[c("PolyPhen_format", "PolyPhen_prediction")]
csv1 <- read.csv("protein.format.csv", header=TRUE, na.strings=c("","NA"))
csv2 <- read.csv("transcript.format.csv", header=TRUE, na.strings=c("","NA"))
csv1 <- rbind(csv1, csv2)
csv1 <- add_columns(csv1, csv3, by=c("Provean_format"))
csv1 <- add_columns(csv1, csv4, by=c("PolyPhen_format"))
write.table(csv1, file="all.csv", sep=",", row.names = FALSE)
file.move(c("protein.format.csv", "protein.result.provean.tsv", "transcript.format.csv", 
            "transcript.result.provean.tsv", "all.result.polyphen.csv"), 
          "./cache/s_prediction_merge", overwrite=TRUE)
print("Combining all data...")
Sys.sleep(1)

csv1 <- read.csv("all.csv")
csv1$"dbSNP_confirmed_with_NCBI" <- gsub('Not_Found', 
                                         "Check variant location @ https://www.ncbi.nlm.nih.gov/variation/view/, it may be a novel variant", 
                                         csv1$"dbSNP_confirmed_with_NCBI")
csv1$"NCBI_Consequence" <- gsub('Not_Found', 
                                "Check manually, it may be a novel variant or a novel change at a reported location", 
                                csv1$"NCBI_Consequence")
csv1$"Filters_exomes" <- gsub('Not_Found', "NO", csv1$"Filters_exomes")
csv1$"Filters_genomes" <- gsub('Not_Found', "NO", csv1$"Filters_genomes")
csv1$"Flags" <- gsub('Not_Found', NA, csv1$"Flags")
csv1$"NCBI_Annotation" <- gsub('Not_Found', NA, csv1$"NCBI_Annotation")
csv1$"Provean_format" <- gsub('Not_Found', NA, csv1$"Provean_format")
csv1$"PolyPhen_format" <- gsub('Not_Found', NA, csv1$"PolyPhen_format")
csv1$"ClinVar_Clinical_Significance" <- gsub('Not_Found', NA, csv1$"ClinVar_Clinical_Significance")
csv1$"ClinVar_Variation_ID" <- gsub('Not_Found', NA, csv1$"ClinVar_Variation_ID")
csv1$PolyPhen_prediction <- str_replace_all(csv1$PolyPhen_prediction, c("benign"="Benign", 
                                                                        "probably_damaging"="Probably Damaging", "possibly_damaging"="Possibly Damaging"))
csv1$VEP_Annotation <- str_replace_all(csv1$VEP_Annotation, c("missense_variant"="Missense Variant", 
                                                      "stop_gained"="Stop Gained", "inframe_deletion"="Inframe Deletion", "splice_region_variant"="Splice Region Variant", 
                                                      "splice_donor_variant"="Splice Donor Variant", "frameshift_variant"="Frameshift Variant", "splice_acceptor_variant"="Splice Acceptor Variant", 
                                                      "stop_lost"="Stop Lost", "start_lost"="Start Lost", "inframe_insertion"="Inframe Insertion", "protein_altering_variant"="Protein Altering Variant"))
pool <- unique(pool)
csv1 <- add_columns(csv1, pool, by="Sample")
csv1 <- csv1[c("Sample", "Diagnosis", "Gene", "Variant", "dbSNP", "dbSNP_confirmed_with_NCBI", "QUAL", "FILTER", 
               "Genotype", "Genotype_Quality", "Allele_Deapth", "Total_Depth", "Variant_Frequency", "NoiseLevel", 
               "StrandBias", "Uncalled", "INFO", "HGVS_Consequence", "VEP_Annotation", "NCBI_Consequence", "NCBI_Annotation", 
               "Protein_ID", "Protein_Consequence", "Transcript_ID", "Transcript_Consequence", "Filters_exomes", 
               "Filters_genomes", "Global_AF", "African_African_American_AF", "Latino_Admixed_American_AF", "Ashkenazi_Jewish_AF", "East_Asian_AF", 
               "South_Asian_AF", "European_Finnish_AF", "European_Non_Finnish_AF", "Other_AF", "Source", "Flags", "ClinVar_Clinical_Significance", 
               "ClinVar_Variation_ID", "Provean_format", "PolyPhen_format", "Provean_prediction", "SIFT_prediction", "PolyPhen_prediction")]
Sys.sleep(1)


# Count
print("Exporting aggregate report...")
csv3 <- filter(csv1, Total_Depth >=50) # good
names(csv3) <- gsub("Total_Depth", "Total Depth", names(csv3))
csv4 <- filter(csv1, Total_Depth <50) # rerun
names(csv4) <- gsub("Total_Depth", "Total Depth", names(csv4))

csv2 <- csv1[c("Sample", "Diagnosis")]
csv2 <- unique(csv2)
csv2 <- csv2[c("Diagnosis")]
variantcount <- table(csv2)
variantcount <- as.data.frame(variantcount)
colnames(variantcount) <- c("Diagnosis", "Patients_with_variants")

samples <- as.data.frame(samplesinthisbatch)
samples <- separate(samples, "samplesinthisbatch", into=c("Sample", "other"), sep="_")
samples <- add_columns(samples, pool, by="Sample")
samples <- unique(samples)
samples <- samples[c("Diagnosis")]
samplecount <- table(samples)
samplecount <- as.data.frame(samplecount)
colnames(samplecount) <- c("Diagnosis", "Patients_in_this_analysis")

csv1 <- add_columns(samplecount, variantcount, by="Diagnosis")
csv1 <- filter(csv1, Patients_in_this_analysis != 0)
csv1 <- mutate(csv1, Percentage=Patients_with_variants/Patients_in_this_analysis*100)
colnames(csv1) <- c("Diagnosis", "Patients in this analysis", "Patients with variants", "Percentage (MAF<1%, including re-run)")


# Export
csv2 <- filter(csv3 ,Global_AF<0.005 & African_African_American_AF<0.005 & Latino_Admixed_American_AF<0.005 & 
                 Ashkenazi_Jewish_AF<0.005 & East_Asian_AF<0.005 & South_Asian_AF<0.005 & 
                 European_Finnish_AF<0.005 & European_Non_Finnish_AF<0.005 & Other_AF<0.005)
csv3 <- filter(csv3, Global_AF>=0.005 | African_African_American_AF>=0.005 | Latino_Admixed_American_AF>=0.005 | 
                 Ashkenazi_Jewish_AF>=0.005 | East_Asian_AF>=0.005 | South_Asian_AF>=0.005 | 
                 European_Finnish_AF>=0.005 | European_Non_Finnish_AF>=0.005 | Other_AF>=0.005)

colnames <- c("Sample", "Diagnosis", "Gene", "Variant", "dbSNP", "dbSNP confirmed with NCBI", "Quality", "Filter", 
              "Genotype", "Genotype Quality", "Allele Deapth", "Total_Depth", "Variant Frequency", "NoiseLevel", 
              "StrandBias", "Uncalled", "Other Info", "HGVS Consequence", "VEP Annotation", "NCBI Consequence", "NCBI Annotation", 
              "Protein ID", "GnomAD Protein Consequence", "Transcript ID", "GnomAD Transcript Consequence", "Filters exomes", 
              "Filters genomes", "Global AF", "African/African American AF", "Latino/Admixed American AF", "Ashkenazi Jewish AF", "East Asian AF", 
              "South Asian AF", "European (Finnish) AF", "European (non-Finnish) AF", "Other AF", "AF data source", "Flags", "ClinVar Clinical Significance", 
              "ClinVar Variation ID", "Provean format", "PolyPhen format", "Provean prediction", "SIFT prediction", "PolyPhen prediction")
colnames(csv2) <- colnames
colnames(csv3) <- colnames
colnames(csv4) <- colnames
print("Writing excel file...")
Sys.sleep(1)
now <- Sys.time()
filename <- paste0(format(now, "%Y-%m-%d_%H%M_"), "BRStudio_export.xlsx")
sheets <- list("Summary" = csv1, "AF<0.5%"= csv2, "0.5%<AF<1%" = csv3, "Re-run (lowDP)" = csv4, 
               "Unannotated by GnomAD (check manually)" = unannotated, "Large indels (check BAM)" = large_indel)
write_xlsx(sheets, filename)
file.move("all.csv", "./cache/t_all_merge", overwrite=TRUE)
rm("csv1", "csv2", "csv3", "csv4", "samples", "samplecount", "variantcount", "colnames", 
   "sheets", "now", "panel", "pool", "numCores", "filename", "samplesinthisbatch",
   "unannotated", "large_indel", "large_indels")
rm(list = ls())
print("All done! Good luck analysing!")
Sys.sleep(1)
