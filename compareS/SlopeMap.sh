#!/bin/bash

## Alaska stations
region="Alaska"
Trange="-0.3/0.3"
minNum=1000

gmt begin $region png
gmt set COLOR_FOREGROUND RED4
gmt set COLOR_BACKGROUND BLUE4
gmt basemap -Rg -JA-160/60/90/5i -Bg
gmt coast -Wthinnest -Df -A10000 -Ggrey -Sgrey
gmt makecpt -Cpolar -T$Trange -Ic -Z -H > icpt.cpt
gawk -F , '{if($5>"'"$minNum"'") print $4,$3,$6}' ${region}_evtInfo.csv | gmt plot -Sa0.4c -Cicpt.cpt
gmt colorbar -Cicpt.cpt -Ba0.1 -Dn0.5/-0.05+jCM+w10c/0.3c+e+h
gmt end

## USArray_WUS stations
region="USArray_WUS"
Trange="-0.3/0.3"
minNum=1000

gmt begin $region png
gmt set COLOR_FOREGROUND RED4
gmt set COLOR_BACKGROUND BLUE4
gmt basemap -Rg -JA-100/55/90/5i -Bg
gmt coast -Wthinnest -Df -A10000 -Ggrey -Sgrey
gmt makecpt -Cpolar -T$Trange -Ic -Z -H > icpt.cpt
gawk -F , '{if($5>"'"$minNum"'") print $4,$3,$6}' ${region}_evtInfo.csv | gmt plot -Sa0.4c -Cicpt.cpt
gmt colorbar -Cicpt.cpt -Ba0.1 -Dn0.5/-0.05+jCM+w10c/0.3c+e+h
gmt end

## USArray_EUS stations
region="USArray_EUS"
Trange="-0.3/0.3"
minNum=1000

gmt begin $region png
gmt set COLOR_FOREGROUND RED4
gmt set COLOR_BACKGROUND BLUE4
gmt basemap -Rg -JA-100/55/90/5i -Bg
gmt coast -Wthinnest -Df -A10000 -Ggrey -Sgrey
gmt makecpt -Cpolar -T$Trange -Ic -Z -H > icpt.cpt
gawk -F , '{if($5>"'"$minNum"'") print $4,$3,$6}' ${region}_evtInfo.csv | gmt plot -Sa0.4c -Cicpt.cpt
gmt colorbar -Cicpt.cpt -Ba0.1 -Dn0.5/-0.05+jCM+w10c/0.3c+e+h
gmt end

## Europe stations
region="Europe"
Trange="-0.1/0.1"
minNum=1000

gmt begin $region png
gmt set COLOR_FOREGROUND RED4
gmt set COLOR_BACKGROUND BLUE4
gmt basemap -Rg -JA20/70/90/5i -Bg
gmt coast -Wthinnest -Df -A10000 -Ggrey -Sgrey
gmt makecpt -Cpolar -T$Trange -Ic -Z -H > icpt.cpt
gawk -F , '{if($5>"'"$minNum"'") print $4,$3,$6}' ${region}_evtInfo.csv | gmt plot -Sa0.4c -Cicpt.cpt
gmt colorbar -Cicpt.cpt -Ba0.1 -Dn0.5/-0.05+jCM+w10c/0.3c+e+h
gmt end
