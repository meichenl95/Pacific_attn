#!/home/meichen/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import HuberRegressor, LinearRegression

region    = "Indonesia"
USArraySR = np.genfromtxt('{}.sr'.format(region))
events    = pd.read_csv("/home/meichen/Research/Pacific_attn/events_2000_2020_6.5",sep=",",skipinitialspace=True)
eventId   = events['id']
eventLat  = events['latitude']
eventLon  = events['longitude']

## ---- Exclude no results stations
USArraySR = USArraySR[abs(USArraySR[:,9])>0.0001]
evla   = USArraySR[:,0]
evlo   = USArraySR[:,1]
evdp   = USArraySR[:,2]
stla1  = USArraySR[:,3]
stlo1  = USArraySR[:,4]
gcarc1 = USArraySR[:,5]
stla2  = USArraySR[:,6]
stlo2  = USArraySR[:,7]
gcarc2 = USArraySR[:,8]
dtstar = USArraySR[:,9]
rmean1 = USArraySR[:,10]
rmax1  = USArraySR[:,11]
rmean2 = USArraySR[:,12]
rmax2  = USArraySR[:,13]
dgcarc = gcarc1 - gcarc2
## ---- Convert gcarc difference to positive
dtstar = dtstar * np.sign(dgcarc)
dgcarc = dgcarc * np.sign(dgcarc)

def main():
    ## ---- All events recorded in the region
    fig,ax = plt.subplots(1,1,figsize=[6,4])
    uniqEvt = np.unique([evlo,evla],axis=1).T
    ax.scatter(uniqEvt[:,0],uniqEvt[:,1],s=1,c='k')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    fig.savefig('event.png')
    ## ---- For each event, fit a slope with uncertainty, record the number of pairs
    numPairs = np.zeros(uniqEvt.shape[0])
    slope    = np.zeros(uniqEvt.shape[0])
    std      = np.zeros(uniqEvt.shape[0])
    eventName = []
    for i in np.arange(uniqEvt.shape[0]):
        traceIdx = np.where((evlo==uniqEvt[i,0]) & (evla==uniqEvt[i,1]))[0]
        nameIdx = np.where((abs(eventLon-uniqEvt[i,0])<0.001) & (abs(eventLat-uniqEvt[i,1])<0.001))[0][0]
        eventName.append(eventId[nameIdx])
        try:
#            p,cov = np.polyfit(dgcarc[traceIdx],dtstar[traceIdx],1,cov=True)
             huber = HuberRegressor(epsilon=1.35).fit(dgcarc[traceIdx].reshape(-1,1), dtstar[traceIdx])
             p = [huber.coef_,huber.intercept_]
        except:
            numPairs[i] = 0
            slope[i] = 0
            std[i] = 0
            print("Failed to perform linear fit on {}".format(eventID[nameIdx]))
            continue
        numPairs[i] = len(traceIdx)
        slope[i] = p[0]
#        std[i] = np.sqrt(np.diag(cov))[0]
        plotCountMap(traceIdx,eventId[nameIdx],p,numPairs[i],slope[i],region)

    ## ---- Save to csv
    df = pd.DataFrame()
    df['eventId']   = pd.Series(eventName)
    df['eventLat']  = pd.Series(uniqEvt[:,1])
    df['eventLon']  = pd.Series(uniqEvt[:,0])
    df['numPairs']  = pd.Series(numPairs)
    df['slope']     = pd.Series(slope)
    df['std']       = pd.Series(std)
    df.to_csv("{}_evtInfo.csv".format(region),float_format="%.3f")

def func(x,slope,intercept):
    return x*slope+intercept

def plotCountMap(traceIdx,eventName,p,numPairs,slope,region):
    fig,ax = plt.subplots(1,1,figsize=[6,4])
    x = np.arange(0.5,25,1)
    y = np.arange(-4.0,4.1,0.2)
    X,Y = np.meshgrid(x,y)                                                      
    Z = np.zeros((X.shape))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            idx = np.where((dgcarc[traceIdx]>x[i]-0.5) & (dgcarc[traceIdx]<x[i]+0.5) & (dtstar[traceIdx]>y[j]-0.1) & (dtstar[traceIdx]<y[j]+0.1))[0]
            Z[j,i] = len(idx)
    cm = ax.imshow(Z,aspect='auto',cmap='Reds',origin='bottom')
    cbar = plt.colorbar(cm,ax=ax)  
    cbar.ax.set_title('log N',size=10)
    cbar.ax.tick_params(labelsize=6)  
    ax.set_yticks([0,5,10,15,20,25,30,35,40])
    ax.set_yticklabels([-4,-3,-2,-1,0,1,2,3,4])
    ax.set_xlabel(r'$\delta\Delta$',size=12)   
    ax.set_ylabel(r'$\delta t^*$',size=12)     
    ax.set_title('Event_{}, num={:.0f}'.format(eventName,numPairs),size=11)      
    ax.tick_params(axis='both',which='both',labelsize=8)
    ax.plot([0,5],[5*(p[1])+20,5*(5*p[0]+p[1])+20],lw=1,c='k',ls='--',label="slope={:.3f}".format(slope))
    ax.legend()
    plt.savefig("countMap/{}_{}.png".format(region,eventName))
    plt.close()

main()
