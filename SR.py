#!/home/meichen/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import obspy
from obspy.taup import TauPyModel
from scipy.signal import hilbert
from scipy.signal import tukey
from mtspec import mtspec
from sklearn.metrics import r2_score
from scipy.interpolate import interp1d

def rsm(x,y,freqmin,freqmax):                                               
    for i in np.arange(len(x)):
        if x[i] == 0.0:
            x[i] = x[i+1]
    f = interp1d(x,y)
    xn = np.linspace(freqmin,freqmax,num=800)
    yn = f(xn)
    return xn,yn

def distaz2latlon(gcarc,az,inlat,inlon):
    ct       = np.cos(inlat*np.pi/180.0)
    st       = np.sin(inlat*np.pi/180.0)
    cp       = np.cos(inlon*np.pi/180.0)
    sp       = np.sin(inlon*np.pi/180.0)
    cb       = np.cos(az*np.pi/180.0)
    sb       = np.sin(az*np.pi/180.0)
    cd       = np.cos(gcarc*np.pi/180.0)
    sd       = np.sin(gcarc*np.pi/180.0)
    ez       = cd * st + sd * cb * ct
    ey       = cd * ct * cp + sd * (-cb * st * cp - sb * sp)
    ex       = cd * ct * sp + sd * (-cb * st * sp + sb * cp)

    outlat   = np.arctan2(ez,np.sqrt(ex*ex + ey*ey)) * 180/np.pi
    outlon   = np.arctan2(ex,ey) * 180/np.pi

    return outlat,outlon

def smooth(y,W):                     
    yy = np.zeros(y.shape)           
    for i in np.arange(W):           
        yy[i]  = np.mean(y[0:i+W])   
        yy[-i-1] = np.mean(y[-i-W::])
    for i in np.arange(W,len(y)-W):  
        yy[i]  = np.mean(y[i-W:i+W]) 
                                     
    return yy

def specRatio(tr,minTr,freqmin,freqmax):

    ## ---- cut waveforms, normalize and taper using tukey windows
    wftr     = tr.data[tr.stats.sac.NS-int(25/tr.stats.delta):tr.stats.sac.NS+int(25/tr.stats.delta)]/tr.data[tr.stats.sac.NS]
    wftr     = wftr * tukey(len(wftr),alpha=0.1)
    wfminTr  = minTr.data[minTr.stats.sac.NS-int(25/tr.stats.delta):minTr.stats.sac.NS+int(25/minTr.stats.delta)]/minTr.data[minTr.stats.sac.NS]
    wfminTr  = wfminTr * tukey(len(wfminTr),alpha=0.1)

    ## ---- multitatper fourier transform for the trace
    ntr          = (int) (2**np.ceil(np.log2(wftr.shape[0]*10)))
    amptr,freqtr = mtspec(data=wftr,delta=tr.stats.delta,time_bandwidth=1.5,number_of_tapers=2,nfft=ntr,adaptive=False)
    amptr        = np.sqrt(amptr)
    freqtr,amptr = rsm(freqtr,amptr,freqmin,freqmax)
    amptr        = smooth(amptr,100)
    ## ---- multitatper fourier transform for the minTrace
    nminTr             = (int) (2**np.ceil(np.log2(wfminTr.shape[0]*10)))
    ampminTr,freqminTr = mtspec(data=wfminTr,delta=minTr.stats.delta,time_bandwidth=1.5,number_of_tapers=2,nfft=nminTr,adaptive=False)
    ampminTr           = np.sqrt(ampminTr)
    freqminTr,ampminTr = rsm(freqminTr,ampminTr,freqmin,freqmax)
    ampminTr           = smooth(ampminTr,100)

    ## ---- calculate spectral ratio
    ratio      = ampminTr/amptr

    ## ---- fit the slope
    # the weight here is relative, it does not matter if we multiply by 2 or divide by 2.
    slope,intercept = np.polyfit(freqminTr,np.log(ratio),1,w=np.sqrt(amptr+ampminTr))

    return slope/np.pi

def s2n(tr,NS,model):

    data = tr.data

    ## ---- compute signal noise ratio
    if (NS-int(120/tr.stats.delta)>0):
        rMax  = data[NS]/np.max(np.abs(data[NS-int(120/tr.stats.delta):NS-int(60/tr.stats.delta)]))
        rMean  = data[NS]/np.mean(np.abs(data[NS-int(120/tr.stats.delta):NS-int(60/tr.stats.delta)]))
    else:
        rMax = 0
        rMean = 0

    return rMean,rMax

def Smax(tr,model):
    
    data = tr.data
    flag = 1

    ## ---- expected S arrival times from PREM
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['S'])
    if (arrivals == []):
        flag = 0
        return 0,flag
    tS       = arrivals[0].time
    ## ---- to avoid overlapping of S and ScS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['ScS'])
    if (arrivals != [] and np.abs(arrivals[0].time-tS)<25):
        flag = 0
        return 0,flag
    ## ---- to avoid overlapping of S and SKS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['SKS'])
    if (arrivals != [] and np.abs(arrivals[0].time-tS)<25):
        flag = 0;
        return 0,flag
    ## ---- to avoid overlapping of S and sS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['sS'])
    if (arrivals != [] and np.abs(arrivals[0].time-tS)<25):
        flag = 0
        return 0,flag 

    ## ---- find the real S arrival times
    NS       = int((tS-tr.stats.sac.b)/tr.stats.delta)
    Nwin     = int(25 / tr.stats.delta)
    if (NS-Nwin<0 or NS+Nwin>=len(data)):
        return 0,0
    N        = np.argmax(np.abs(data[NS-Nwin:NS+Nwin]))
    NS       = NS - Nwin + N

    return NS,flag

def inRegion(tr,Regionname):
    if(Regionname == "Alaska"):
        latmin = 50
        latmax = 72
        lonmin = 174
        lonmax = -130+360
    elif(Regionname == "USArray_WUS"):
        latmin = 25
        latmax = 50
        lonmin = -125+360
        lonmax = -110+360
    elif(Regionname == "USArray_EUS"):
        latmin = 25
        latmax = 50
        lonmin = -95+360
        lonmax = -65+360
    elif(Regionname == "Europe"):   
        latmin = 30                 
        latmax = 70                 
        lonmin = -10+360
        lonmax = 30                 
    elif(Regionname == "Indonesia"):
        latmin = -9                 
        latmax = 6                  
        lonmin = 95                 
        lonmax = 141

    stla = tr.stats.sac.stla
    stlo = tr.stats.sac.stlo+360 if tr.stats.sac.stlo<0 else tr.stats.sac.stlo
    if (stla>=latmin and stla<=latmax and stlo>=lonmin and stlo<=lonmax):
        return True
    elif (stla>=latmin and stla<=latmax and stlo>=lonmin and lonmin>lonmax):
        return True
    elif (stla>=latmin and stla<=latmax and stlo<=lonmax and lonmin>lonmax):
        return True
    else:
        return False
        

def main():

    events   = pd.read_csv('/home/meichen/Research/Pacific_attn/events_temp',sep=',')
 ##   eventids = ['usp000fsnz']
 ##   eventmag = [6.7]
    eventids = events['id']
    eventmag = events['mag']
    model    = obspy.taup.TauPyModel(model='prem')
    ## ---- the region names include "Alaska" "USArray_WUS" "USArray_EUS" "Europe" "China" "Indonesia"
    RegionName = "USArray_WUS"
    freqmin  = 0.02
    freqmax  = 0.8

    ## ---- initialize dtstar, evla,evlo,stla,stlo,evdp,rfla,rflo
    dtstar = []
    evla   = []
    evlo   = []
    stla1  = []
    stlo1  = []
    gcarc1 = []
    rmean1 = []
    rmax1  = []
    evdp   = []
    stla2  = []
    stlo2  = []
    gcarc2 = []
    rmean2 = []
    rmax2  = []

    ## ---- loop events
    for l,ids in enumerate(eventids):
        path_sacfile = '/home/meichen/work1/Pacific_attn/event_{}/waveforms/SAC_files'.format(ids)
        st = obspy.read('{}/*.BHT.SAC.dis'.format(path_sacfile))

        ## ---- find the trace with minimum gcarc
        st.filter('bandpass',freqmin=0.02,freqmax=0.8,zerophase=True,corners=2)
        for tr in st:
            if (not inRegion(tr,RegionName)):
                st.remove(tr)
                continue
            if (not np.all(tr.data)):
                st.remove(tr)
                continue
            NS,flag = Smax(tr,model)
            if(flag == 0):
                st.remove(tr)
                continue
            rMean,rMax = s2n(tr,NS,model)
            if (rMax<2 or rMean<4):
                st.remove(tr)
                continue
            tr.stats.sac.rMean = rMean
            tr.stats.sac.rMax = rMax
            tr.stats.sac.NS = NS

        ## ---- at least two traces in a stream
        if (st.count()<2):
            continue

        ## ---- compute spectral ratios
        for i in np.arange(st.count()-1):
            for j in np.arange(i,st.count()):
                tr1 = st[i]
                tr2 = st[j]
                ## ---- calculate spectral ratios
                slope = specRatio(tr1,tr2,freqmin,freqmax)
    
                # save results
                print("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}".format(tr1.stats.sac.evla,tr1.stats.sac.evlo,tr1.stats.sac.evdp,tr1.stats.sac.stla,tr1.stats.sac.stlo,tr1.stats.sac.gcarc,tr2.stats.sac.stla,tr2.stats.sac.stlo,tr2.stats.sac.gcarc,slope,tr1.stats.sac.rMean,tr1.stats.sac.rMax,tr2.stats.sac.rMean,tr2.stats.sac.rMax))
                dtstar.append(slope)
                evla.append(tr1.stats.sac.evla)
                evlo.append(tr1.stats.sac.evlo)
                evdp.append(tr1.stats.sac.evdp)
                stla1.append(tr1.stats.sac.stla)
                stlo1.append(tr1.stats.sac.stlo)
                gcarc1.append(tr1.stats.sac.gcarc)
                rmean1.append(tr1.stats.sac.rMean)
                rmax1.append(tr1.stats.sac.rMax)
                stla2.append(tr2.stats.sac.stla)
                stlo2.append(tr2.stats.sac.stlo)
                gcarc2.append(tr2.stats.sac.gcarc)
                rmean2.append(tr2.stats.sac.rMean)
                rmax2.append(tr2.stats.sac.rMax)

    np.savez("output",dtstar=dtstar,evla=evla,evlo=evlo,evdp=evdp,stla1=stla1,stlo1=stlo1,gcarc1=gcarc1,stla2=stla2,stlo2=stlo2,gcarc2=gcarc2,rmean1=rmean1,rmax1=rmax1,rmean2=rmean2,rmax2=rmax2)

main()
