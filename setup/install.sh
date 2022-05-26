#!/bin/bash
sudo apt-get update
sudo apt-get dist-upgrade -y
sudo apt-get upgrade -y
sudo apt install default-jre default-jdk libxml2-dev libssl-dev libcurl4-openssl-dev libnlopt-dev -y

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list

sudo apt-get update
sudo apt-get install r-base -y
sudo R CMD javareconf
sudo Rscript ./packages.R
Rscript ./webdriver.R

wget http://archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.0g-2ubuntu4_amd64.deb
wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2022.02.2-485-amd64.deb
sudo apt-get install -f -y ./libssl1.1_1.1.0g-2ubuntu4_amd64.deb ./rstudio-2022.02.2-485-amd64.deb
rm ./libssl1.1_1.1.0g-2ubuntu4_amd64.deb ./rstudio-2022.02.2-485-amd64.deb

sudo sed -i 's+Exec=/usr/lib/rstudio/bin/rstudio %F+Exec=/usr/lib/rstudio/bin/rstudio --no-sandbox %F+g' /usr/share/applications/rstudio.desktop
sudo chmod +x /usr/share/applications/rstudio.desktop

sudo apt-get autoremove -y
sudo apt-get clean
reboot