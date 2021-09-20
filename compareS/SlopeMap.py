#!/home/meichen/anaconda3/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plotSlopeMap(region):
    f = pd.read_csv("{}_evtInfo.csv".format(region),skipinitialspace=True)
    index = np.where((np.array(f['numPairs'])>1000))[0]
    lat = f['eventLat'][index]
    lon = f['eventLon'][index]
    slope = f['slope'][index]

    fig,ax = plt.subplots(1,1,figsize=[6,4])
    cbar = ax.scatter(lon,lat,c=slope,cmap="bwr",vmin=-0.05,vmax=0.05,lw=0.5,edgecolor="k",alpha=0.5)
    ax.set_xlabel("Longitude",size=12)
    ax.set_ylabel("Latitude",size=12)
    ax.set_title("{}".format(region),size=14)
    plt.colorbar(cbar,ax=ax)
    ax.grid()
    plt.show()
    plt.close()

def main():
    regionList = ["Alaska", "Europe", "USArray_WUS", "USArray_EUS", "Indonesia"]
    
    for region in regionList:
        plotSlopeMap(region)
    
main()
