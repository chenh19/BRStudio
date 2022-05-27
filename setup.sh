#!/bin/bash
# BRStudio setup v0.1.1
# This shell script is intented for Kubuntu 22.04 LTS and other linux systems based on Ubuntu 22.04 LTS
# To set up the environment for BRStudio, simply ```bash ./setup.sh``` in terminal

# update and install packages required by R and R packages installing
sudo apt-get update
sudo apt-get dist-upgrade -y
sudo apt-get upgrade -y
sudo apt install default-jre default-jdk libxml2-dev libssl-dev libcurl4-openssl-dev libnlopt-dev -y

# install R from cran
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list
sudo apt-get update
sudo apt-get install r-base -y

# configure java for R
sudo R CMD javareconf

# install R packages and download web scraping driver
echo "install.packages(c('devtools', 'BiocManager', 'tidyverse', 'readxl', 'writexl', 'expss', 'vcfR', 'filesstrings', 'R.utils', 'car', 'foreach', 'doParallel', 'rJava', 'RSelenium'))" > packages.R
echo "wdman::chrome(version = 'latest')" > webdriver.R
sudo Rscript ./packages.R
Rscript ./webdriver.R
rm ./packages.R ./webdriver.R

# install RStudio and packages that no longer available in ubuntu 22.04
wget http://archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.0g-2ubuntu4_amd64.deb
wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2022.02.2-485-amd64.deb
sudo apt-get install -f -y ./libssl1.1_1.1.0g-2ubuntu4_amd64.deb ./rstudio-2022.02.2-485-amd64.deb
rm ./libssl1.1_1.1.0g-2ubuntu4_amd64.deb ./rstudio-2022.02.2-485-amd64.deb

# configure RStudio exec command in desktop file so that it can run
sudo sed -i 's+Exec=/usr/lib/rstudio/bin/rstudio %F+Exec=/usr/lib/rstudio/bin/rstudio --no-sandbox %F+g' /usr/share/applications/rstudio.desktop
sudo chmod +x /usr/share/applications/rstudio.desktop

# cleanup
sudo apt-get autoremove -y
sudo apt-get clean
if [ -e "examples" ]; then rm -rf "examples"; fi
if [ -f ./LICENSE ]; then rm ./LICENSE; fi
if [ -f ./README.md ]; then rm ./README.md; fi
if [ -f ./setup.sh ]; then rm ./setup.sh; fi

# output
echo 'All done!'
if [ -f /var/run/reboot-required ]; then echo 'Please reboot the system.'; fi
