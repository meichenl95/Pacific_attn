#!/bin/bash

eventFile="../events_2000_2020_6.5"
stationFile1="USArray.sr2"
stationFile2="Alaska.sr"

gmt begin map png
gmt basemap -Rg -JH8i -Bxa60 -Bya30
gmt coast -Wthinnest -Df -A10000 
gawk -F',' 'NR>1{print $3,$2}' $eventFile | gmt plot -Sa0.3c -Wblack -Gred
gawk '{print $8,$7}' $stationFile1 | uniq | gmt plot -Sc0.05c -Wblack -Glightblue
gawk '{print $8,$7}' $stationFile2 | uniq | gmt plot -Sc0.05c -Wblack -Glightblue
gmt end
