#!/home/meichen/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import obspy
from obspy.taup import TauPyModel
from scipy.signal import hilbert
from scipy.signal import tukey

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

def instan_freq(tr,tstar):

    data = tr.data

    ## ---- cut waveforms, normalize and taper using tukey windows
    wf   = data[tr.stats.sac.NS-int(25/tr.stats.delta):tr.stats.sac.NS+int(25/tr.stats.delta)].copy()
    wf   = wf/np.max(wf)
    wf   = wf * tukey(len(wf),alpha=0.1)

    ## ---- attenuate waveforms
#    if (tstar!=0):
    if (1):
        nwf       = (int) (2**np.ceil(np.log2(wf.shape[0]*10)))
        amp       = np.fft.rfft(wf,n=nwf)
        freq      = np.fft.rfftfreq(nwf,d=tr.stats.delta)
        omega     = 2*np.pi*freq
        omegar    = 2*np.pi
        D         = np.ones(freq.shape) + 0j
        index     = np.where(abs(freq)>0)[0]
#        D[index]  = np.exp(-omega[index]*tstar/2-1j*omega[index]*tstar/np.pi*np.log(omega[index]/omegar))
        D[index]  = np.exp(-omega[index]*tstar/2)
        amp       = amp*D
        #print(np.sum(freq*abs(amp)*(freq[1]-freq[0]))/np.sum(abs(amp)*(freq[1]-freq[0])))
        wf        = np.fft.irfft(amp)[0:len(wf)]
    
    ## ---- instantaneous phase
    wfhilb   = hilbert(wf)
    wfquad   = wfhilb.imag

    ## ---- envelop
    wfenv    = np.abs(wfhilb)
    amphilb = np.abs(np.fft.fft(wfhilb,n = len(wfhilb)*100))
    freq = np.fft.fftfreq(len(wfhilb)*100,d=tr.stats.delta)
    #print(np.sum(freq*amphilb*(freq[1]-freq[0]))/np.sum(amphilb*(freq[1]-freq[0])))

    ## ---- instantaneous frequency
    wfdydt      = np.zeros(wf.shape)
    wfdyqdt     = np.zeros(wf.shape)
    # first point
    wfdydt[0]   = (wf[1]-wf[0])/tr.stats.delta
    wfdyqdt[0]  = (wfquad[1]-wfquad[0])/tr.stats.delta
    # last point
    wfdydt[-1]  = (wf[-1]-wf[-2])/tr.stats.delta
    wfdyqdt[-1] = (wfquad[-1]-wfquad[-2])/tr.stats.delta
    # other points
    for i in np.arange(1,wf.shape[0]-1):
        wfdydt[i]    = (wf[i+1]-wf[i-1])/tr.stats.delta/2
        wfdyqdt[i]   = (wfquad[i+1]-wfquad[i-1])/tr.stats.delta/2
    # calculate instantaneous frequency list
    wfInstF     = (wf*wfdyqdt-wfquad*wfdydt) / (2*np.pi*(wf**2+wfquad**2 + np.max(wf**2+wfquad**2)*0.001))
    # smooth
#    W          = int(5/tr.stats.delta)
    W          = 1
    wfInstFsm  = np.zeros(wfInstF.shape)
    for i in np.arange(1,W):
        wfInstFsm[i]  = np.sum(wfInstF[1:i+W]*wfenv[1:i+W]**2)/np.sum(wfenv[1:i+W]**2)
        wfInstFsm[-i] = np.sum(wfInstF[-1-i-W::]*wfenv[-1-i-W::]**2)/np.sum(wfenv[-1-i-W::]**2)
    for i in np.arange(W,wfInstF.shape[0]-W+1):
        wfInstFsm[i]    = np.sum(wfInstF[i-W:i+W]*wfenv[i-W:i+W]**2)/np.sum(wfenv[i-W:i+W]**2)
    wfInstFsm[0] = wfInstFsm[1]
    wfInstFsm[-1] = wfInstFsm[-2]

    ## ---- instantaneous frequency at envelope peaks
    InstF        = wfInstFsm[np.argmax(wfenv)]
    #print(InstF)

    return InstF

def s2n(tr,NS):

    data = tr.data

    ## ---- compute signal noise ratio
    if (NS-int(120/tr.stats.delta)>0):
        rMax  = data[NS]/np.max(np.abs(data[NS-int(120/tr.stats.delta):NS-int(60/tr.stats.delta)]))
        rMean = data[NS]/np.mean(np.abs(data[NS-int(120/tr.stats.delta):NS-int(60/tr.stats.delta)]))
    else:
        rMax = 0
        rMean = 0

    return rMean,rMax

def Smax(tr,model):
    
    data = tr.data
    flag = 1

    ## ---- expected S and ScS arrival times from PREM
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['S'])
    if (arrivals == []):
        flag = 0
        return 0,flag
    tS       = arrivals[0].time
    ## ---- to avoid overlapping of S and ScS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['ScS'])
    if (arrivals!=[] and np.abs(arrivals[0].time-tS)<25):
        flag = 0
        return 0,flag
    ## ---- to avoid overlapping of S and SKS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['SKS'])
    if (arrivals!=[] and np.abs(arrivals[0].time-tS)<25):
        flag = 0
        return 0,flag
    ## ---- to avoid overlapping of S and sS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['sS'])
    if (arrivals!=[] and np.abs(arrivals[0].time-tS)<25):
        flag = 0
        return 0,flag

    ## ---- find the real S arrival times
    NS       = int((tS-tr.stats.sac.b)/tr.stats.delta)
    Nwin     = int(25 / tr.stats.delta)
    N        = np.argmax(np.abs(data[NS-Nwin:NS+Nwin]))
    NS       = NS - Nwin + N

    return NS,flag

def inRegion(tr,Regionname):                                                
    if(Regionname == "Alaska"):
        latmin = 50
        latmax = 72
        lonmin = 174
        lonmax = -130+360
    elif(Regionname == "USArray"):
        latmin = 25
        latmax = 50
        lonmin = -125+360
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
    else:
        return False

def main():

    events   = pd.read_csv('/home/meichen/Research/Pacific_attn/events_temp',sep=',')
 ##   eventids = ['usp000fsnz']
 ##   eventmag = [6.7]
    eventids = events['id']
    eventmag = events['mag']
    model    = obspy.taup.TauPyModel(model='prem')
    RegionName = "Alaska"

    ## ---- initialize dtstar, evla,evlo,stla,stlo,gcarc,evdp,stlaminTr,stlominTr,gcarcminTr
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
        # According to trials, the filtering frequency range doesn't affect the dtstar significantly.
        st.filter('bandpass',freqmin=0.02,freqmax=0.8,zerophase=True,corners=2)
        for tr in st:
            if(not inRegion(tr,RegionName)):
                st.remove(tr)
                continue
            if (not np.all(tr.data)):
                st.remove(tr)
                continue
            NS,flag = Smax(tr,model)
            if (flag == 0):
                st.remove(tr)
                continue
            rMean,rMax = s2n(tr,NS)
            if (rMax<2 or rMean<4):
                st.remove(tr)
                continue
            tr.stats.sac.rMean = rMean
            tr.stats.sac.rMax = rMax
            tr.stats.sac.NS = NS

        ## ---- compute instantaneous frequency of minTr
        if (st.count()<2):
            continue

        for i in np.arange(st.count()-1):
            tr1 = st[i]
            InstFtr1 = instan_freq(tr1,0)
            for j in np.arange(i,st.count()):
                tr2 = st[j]
                InstFtr2 = instan_freq(tr2,0)
                ## ---- matching tstar
                # trace 1 has a higher instantaneous frequency
                if (InstFtr2<InstFtr1):
                    n   = 1
                    dt  = 0.0
                    ddt = 1.0
                    while(n<=100):
                        dt  = dt + ddt
                        InstFtr1_temp = instan_freq(tr1,dt)
                        n = n + 1
                        if (abs(InstFtr1_temp-InstFtr2)<0.00001):
                            break
                        elif (InstFtr1_temp<InstFtr2):
                            dt  = dt - ddt
                            ddt = ddt/10.0
                    dt = -dt
                # minTr has a higher instantaneous frequency
                if (InstFtr1<=InstFtr2):
                    n   = 1
                    dt  = 0.0
                    ddt = 1.0
                    while(n<=100):
                        dt  = dt + ddt
                        InstFtr2_temp = instan_freq(tr2,dt)
                        n = n + 1
                        if (abs(InstFtr2_temp-InstFtr1)<0.00001):
                            break
                        elif (InstFtr2_temp<InstFtr1):
                            dt  = dt - ddt
                            ddt = ddt/10.0
                if (n>100):
                    dt    = 0
                # save results
                print("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}".format(tr1.stats.sac.evla,tr1.stats.sac.evlo,tr1.stats.sac.evdp,tr1.stats.sac.stla,tr1.stats.sac.stlo,tr1.stats.sac.gcarc,tr2.stats.sac.stla,tr2.stats.sac.stlo,tr2.stats.sac.gcarc,dt,tr1.stats.sac.rMean,tr1.stats.sac.rMax,tr2.stats.sac.rMean,tr2.stats.sac.rMax))
                dtstar.append(dt)
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
