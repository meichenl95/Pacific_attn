#!/home/meichen/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import obspy
from obspy.taup import TauPyModel
from scipy.signal import hilbert
from scipy.signal import tukey
from mtspec import mtspec
from scipy.interpolate import interp1d

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

def rsm(x,y,freq):
    for i in np.arange(len(x)):
        if x[i] == 0.0:
            x[i] = x[i+1]
    f = interp1d(x,y)
    yn = f(freq)
    return yn

def sr(tr,envNS,envNScS,leftlen,rightlen,freq):

    data = tr.data

    ## ---- cut waveforms, normalize and taper using tukey windows
    cutS    = data[envNS-leftlen:envNS+rightlen]/data[envNS]
    cutScS  = data[envNScS-leftlen:envNScS+rightlen]/data[envNS]
    cutN    = data[envNS-int(150/tr.stats.delta):envNS-int(50/tr.stats.delta)]/data[envNS]

    ## ---- fourier transform
    # S wave
    nfft       = (int) (2**np.ceil(np.log2(cutS.shape[0]*10)))
    ampS,frequency   = mtspec(data=cutS,delta=tr.stats.delta,time_bandwidth=1.5,number_of_tapers=2,nfft=nfft)
    ampS       = np.sqrt(ampS)
    ampS       = rsm(frequency,ampS,freq)
    # ScS wave
    ampScS,frequency = mtspec(data=cutScS,delta=tr.stats.delta,time_bandwidth=1.5,number_of_tapers=2,nfft=nfft)
    ampScS     = np.sqrt(ampScS)
    ampScS     = rsm(frequency,ampScS,freq)
    # noise
    nfftN      = (int) (2**np.ceil(np.log2(cutN.shape[0]*10)))
    ampN,frequency = mtspec(data=cutN,delta=tr.stats.delta,time_bandwidth=1.5,number_of_tapers=2,nfft=nfftN)
    ampN       = np.sqrt(ampN)
    ampN       = rsm(frequency,ampN,freq)

    return ampS,ampScS,ampN
#    plt.plot(np.arange(len(cutS))*tr.stats.delta,cutS,'b-',lw=1,label='S')
#    plt.plot(np.arange(len(cutScS))*tr.stats.delta,cutScS,'r-',lw=1,label='ScS')
#    plt.xlabel('Time (s)',size=14)
#    plt.legend()
#    plt.show()
#    plt.close()
#    plt.plot(frequency,np.log(ampS),'b-',lw=1,label='S')
#    plt.plot(frequency,np.log(ampScS),'r-',lw=1,label='ScS')
#    plt.plot(frequency[indx],np.log(ratio)[indx],'k-',lw=1,label='ratio')
#    plt.xlim([0,1.0])
#    plt.xlabel('Frequency (Hz)',size=14)
#    plt.ylabel('log(amplitude)',size=14)
#    plt.legend()
#    plt.show()
#    plt.close()

    return slope/np.pi

def s2n(tr,envNS,envNScS,model):

    data = tr.data
    envelope = np.abs(hilbert(data))
    flag = 1

    ## ---- find S wave duration
    # left length
    i = 1
    while(envelope[envNS-i]>0.5*envelope[envNS] or envelope[envNScS-i]>0.5*envelope[envNScS]):
        i = i + 1
        if envNS-i<0:
            flag = 0
            break
    leftlen = 3*i
    # right length
    i = 1
    while(envelope[envNS+i]>0.5*envelope[envNS] or envelope[envNScS+i]>0.5*envelope[envNScS]):
        i = i + 1
        if envNScS+i>=tr.stats.npts:
            flag = 0
            break
    rightlen = 3*i
                  
    ## ---- exclude large perturbations
    if (leftlen<envNS and rightlen+envNScS<tr.stats.npts and leftlen>3 and rightlen>3):
        if (np.max(np.abs(envelope[envNS-leftlen:envNS-int(leftlen/3)-2]))>0.8*envelope[envNS] or np.max(np.abs(envelope[envNS+int(rightlen/3)+2:envNS+rightlen]))>0.8*envelope[envNS]):
            flag = 0
        if (np.max(np.abs(envelope[envNScS-leftlen:envNScS-int(leftlen/3)-2]))>0.8*envelope[envNScS] or np.max(np.abs(envelope[envNScS+int(rightlen/3)+2:envNScS+rightlen]))>0.8*envelope[envNScS]):
            flag = 0
    else:
        flag = 0

    ## ---- ensure separation of S and ScS
    if (envNScS-envNS<=rightlen+leftlen):
        flag = 0

    ## ---- compute signal noise ratio based on envelope
    dur = leftlen+rightlen
    if (envNS-5*dur>0):
        maxratio1  = envelope[envNS]/np.max(envelope[envNS-5*dur:envNS-2*dur])
        meanratio1 = envelope[envNS]/np.mean(envelope[envNS-5*dur:envNS-2*dur])
        maxratio2  = envelope[envNScS]/np.max(envelope[envNS-5*dur:envNS-2*dur])
        meanratio2 = envelope[envNScS]/np.mean(envelope[envNS-5*dur:envNS-2*dur])
    else:
        maxratio1 = 0
        meanratio1 = 0
        maxratio2 = 0
        meanratio2 = 0
        midratio = 0
        flag = 0

    ## ---- exlude phase interference of sS
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['S','ScS','sS'])
    if (abs(arrivals[2].time-arrivals[0].time)<tr.stats.delta*(leftlen+rightlen) or abs(arrivals[1].time-arrivals[2].time)<tr.stats.delta*(leftlen+rightlen)):
        flag = 0

    return maxratio1,meanratio1,maxratio2,meanratio2,flag,leftlen,rightlen

def envelopemax(tr,model):
    
    data = tr.data
    envelope = np.abs(hilbert(data))
    flag = 1

    ## ---- expected S and ScS arrival times from PREM
    arrivals = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp,distance_in_degree=tr.stats.sac.gcarc,phase_list=['S','ScS'])
    tS       = arrivals[0].time
    tScS     = arrivals[1].time

    ## ---- find the real S and ScS arrival times
    NS       = int((tS-tr.stats.sac.b)/tr.stats.delta)
    NScS     = int((tScS-tr.stats.sac.b)/tr.stats.delta)
    # find S
    Nwin     = int(20 / tr.stats.delta)
    depmax   = envelope[NS]
    envNS    = NS
    for i in np.arange(Nwin):
        if(abs(envelope[NS-i])>depmax):
            depmax    = envelope[NS-i]
            envNS     = NS-i
        if(abs(envelope[NS+i])>depmax):
            depmax    = envelope[NS+i]
            envNS     = NS+i
    # find ScS
    Nwin     = int(20 / tr.stats.delta)
    NScS     = NScS + envNS - NS
    depmax   = envelope[NScS]
    envNScS  = NScS
    for i in np.arange(Nwin):
        if(envelope[NScS-i]>depmax and envelope[NScS-i]>envelope[NScS-i-1] and envelope[NScS-i]>envelope[NScS-i+1]):
            depmax    = envelope[NScS-i]
            envNScS   = NScS-i
        if(envelope[NScS+i]>depmax and envelope[NScS+i]>envelope[NScS+i-1] and envelope[NScS+i]>envelope[NScS+i+1]):
            depmax    = envelope[NScS+i]
            envNScS   = NScS+i
    if (envNScS == NScS):
        flag = 0

    return envNS,envNScS,flag

def main():

    events   = pd.read_csv('events_temp',sep=',')
 ##   eventids = ['usp000fgsv']
 ##   eventmag = [6.7]
    eventids = events['id']
    eventmag = events['mag']
    model    = obspy.taup.TauPyModel(model='prem')

    ## ---- initialize dtstar, evla,evlo,stla,stlo,evdp,rfla,rflo
    dtstar = []
    evla   = []
    evlo   = []
    stla   = []
    stlo   = []
    evdp   = []
    rfla   = []
    rflo   = []
    freq   = np.linspace(0.1,2,num=8000)
    spectrum_S = np.zeros(8000)
    spectrum_ScS = np.zeros(8000)
    spectrum_N = np.zeros(8000)

    ## ---- loop events
    for l,ids in enumerate(eventids):
        path_sacfile = '/home/meichen/work1/Pacific_attn/event_{}/waveforms/SAC_files'.format(ids)
        st = obspy.read('{}/*.BHT.SAC.raw'.format(path_sacfile))

        ## ---- loop traces
        for tr in st:
            ## ---- lowpass
#            tr.filter('bandpass',freqmin=0.02,freqmax=0.12,zerophase=True)

            ## ---- check if nan value exist
            if (not np.all(tr.data)):
                continue

            ## ---- find S and ScS using envelope
            envNS,envNScS,flag = envelopemax(tr,model)
            if flag == 0:
                continue

            ## ---- signal noise ratio
            r1,r2,r3,r4,flag,leftlen,rightlen = s2n(tr,envNS,envNScS,model)
            if r1<3 or r2<6 or r3<2 or r4<4 or flag == 0:
                continue
            if(envNS-int(200/tr.stats.delta)<0 or envNScS+leftlen+rightlen>=tr.stats.npts) or tr.stats.delta*(leftlen+rightlen)<=10:
                continue

            ## ---- calculate spectral ratios
            specS,specScS,specN = sr(tr,envNS,envNScS,leftlen,rightlen,freq)

            ## ---- find reflection point of ScS at CMB
            arrivals = model.get_pierce_points(tr.stats.sac.evdp,tr.stats.sac.gcarc,phase_list=["ScS"])
            for i in np.arange(arrivals[0].pierce.shape[0]):
                if(arrivals[0].pierce[i][3] == 2891.0):
                    reflgcarc  = arrivals[0].pierce[i][2]*tr.stats.sac.gcarc
            if (i == arrivals[0].pierce.shape[0]):
                print("Reflection point wrong")
                continue
            reflat,reflon = distaz2latlon(reflgcarc,tr.stats.sac.az,tr.stats.sac.evla,tr.stats.sac.evlo)

            ## ---- stack spectrum
            spectrum_S = spectrum_S + np.log(specS)
            spectrum_ScS = spectrum_ScS + np.log(specScS)
            spectrum_N = spectrum_N + np.log(specN)

    fig,ax = plt.subplots(1,1,figsize=[5,4])
    ax.plot(freq,spectrum_S,c='k',lw=1,label='S')
    ax.plot(freq,spectrum_ScS,c='r',lw=1,label='ScS')
    ax.plot(freq,spectrum_N,c='b',lw=1,label='Noise')
    ax.legend()
    ax.set_xlabel('Frequency (Hz)',size=12)
    ax.set_ylabel('ln amp',size=12)
    fig.tight_layout()
    plt.savefig('s2n.pdf')

main()
