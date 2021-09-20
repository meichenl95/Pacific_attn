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

def instan_freq(tr,tstar,envNS,envNScS,leftlen,rightlen):

    data = tr.data                 
                                   
    ## ---- cut waveforms, normalize and taper using tukey windows
    cut    = data[envNS-int(200/tr.stats.delta):envNScS+leftlen+rightlen]
                                
    ## ---- attenuate waveforms 
    if (1):              
        nfft       = (int) (2**np.ceil(np.log2(cut.shape[0]*10)))
        amp        = np.fft.rfft(cut,n=nfft)                     
        frequency  = np.fft.rfftfreq(nfft,d=tr.stats.delta)      
        omega      = 2*np.pi*frequency
        omegar     = 2*np.pi          
        D          = np.ones(frequency.shape) + 0j
        index      = np.where(abs(frequency)>0)[0]
        D[index]   = np.exp(-omega[index]*tstar/2-1j*omega[index]*tstar/np.pi*np.log(omega[index]/omegar))
        amp        = amp*D 
        cut        = np.fft.irfft(amp)[0:len(cut)]
                                   
    ## ---- instantaneous phase    
    cuthilb   = hilbert(cut)
    cutquad   = cuthilb.imag
                   
    ## ---- envelop
    cutenv    = np.abs(cuthilb)
                               
    ## ---- instantaneous frequency
    cutdydt      = np.zeros(cut.shape)
    cutdyqdt     = np.zeros(cut.shape)
    # first point
    cutdydt[0]   = (cut[1]-cut[0])/tr.stats.delta
    cutdyqdt[0]  = (cutquad[1]-cutquad[0])/tr.stats.delta
    # last point
    cutdydt[-1]  = (cut[-1]-cut[-2])/tr.stats.delta
    cutdyqdt[-1] = (cutquad[-1]-cutquad[-2])/tr.stats.delta
    # other points
    for i in np.arange(1,cut.shape[0]-1):
        cutdydt[i]    = (cut[i+1]-cut[i-1])/tr.stats.delta/2
        cutdyqdt[i]   = (cutquad[i+1]-cutquad[i-1])/tr.stats.delta/2
    # calculate instantaneous freq
    cutf         = (cut*cutdyqdt-cutquad*cutdydt) / (2*np.pi*(cut**2+cutquad**2 + np.max(cut**2+cutquad**2)*0.001))
    # smooth
    W         = int(6/tr.stats.delta)
    cutfsm   = np.zeros(cutf.shape)
    for i in np.arange(1,W):
        cutfsm[i]  = np.sum(cutf[1:i+W]*cutenv[1:i+W]**2)/np.sum(cutenv[1:i+W]**2)
        cutfsm[-i] = np.sum(cutf[-1-i-W::]*cutenv[-1-i-W::]**2)/np.sum(cutenv[-1-i-W::]**2)
    for i in np.arange(W,cutf.shape[0]-W+1):
        cutfsm[i]    = np.sum(cutf[i-W:i+W]*cutenv[i-W:i+W]**2)/np.sum(cutenv[i-W:i+W]**2)
    cutfsm[0] = cutfsm[1]
    cutfsm[-1] = cutfsm[-2]

    ## ---- instantaneous frequency at envelope peaks
    fS        = cutfsm[np.argmax(cutenv[int(200/tr.stats.delta)-leftlen:int(200/tr.stats.delta)+rightlen])+int(200/tr.stats.delta)-leftlen]
    fScS      = cutfsm[np.argmax(cutenv[int(200/tr.stats.delta)+envNScS-envNS-leftlen:int(200/tr.stats.delta)+envNScS-envNS+rightlen])+int(200/tr.stats.delta)+envNScS-envNS-leftlen]

    return fS,fScS

def cc(tr,tstar,envNS,envNScS,leftlen,rightlen,attnphase):

    data = tr.data

    ## ---- cut waveforms, normalize and taper using tukey windows
    cut    = data
    envelope = np.abs(hilbert(data))

    ## ---- attenuate waveforms
    if (tstar!=0):
        nfft       = (int) (2**np.ceil(np.log2(cut.shape[0]*2)))
        amp        = np.fft.rfft(cut,n=nfft)
        frequency  = np.fft.rfftfreq(nfft,d=tr.stats.delta)
        omega      = 2*np.pi*frequency
        omegar     = 2*np.pi
        D          = np.ones(frequency.shape) + 0j
        index      = np.where(abs(frequency)>0)[0]
        D[index]   = np.exp(-omega[index]*tstar/2-1j*omega[index]*tstar/np.pi*np.log(omega[index]/omegar))
        amp        = amp*D
        cut        = np.fft.irfft(amp)[0:len(cut)]

    ## ---- calculate cross-correlation coefficient
    cut    = np.abs(hilbert(cut))
    if attnphase == "S":
        ampmaxS   = 0
        for ii in np.arange(envNS-leftlen,envNS+rightlen):
            if cut[ii]>cut[ii-1] and cut[ii]>cut[ii+1] and cut[ii]>ampmaxS:
                ampmaxS    = cut[ii]
                maxidxS    = ii
        if(ampmaxS == 0.0):
            return 0.0
        ampmaxScS   = 0
        for ii in np.arange(envNScS-leftlen,envNScS+rightlen):
            if envelope[ii]>envelope[ii-1] and envelope[ii]>envelope[ii+1] and envelope[ii]>ampmaxScS:
                ampmaxScS    = envelope[ii]
                maxidxScS = ii
        if(ampmaxScS ==0.0 ):
            return 0.0
        iil       = 1
        while( cut[maxidxS-iil]>0.5*ampmaxS and envelope[maxidxScS-iil]>0.5*ampmaxScS ):
            iil = iil + 1
        iir       = 1
        while( cut[maxidxS+iir]>0.5*ampmaxS and envelope[maxidxScS+iir]>0.5*ampmaxScS ):
            iir = iir + 1
        cutS        = cut[maxidxS-iil:maxidxS+iir]/cut[maxidxS]
        cutScS      = envelope[maxidxScS-iil:maxidxScS+iir]/envelope[maxidxScS]
    if attnphase == "ScS":
        ampmaxS   = 0
        for ii in np.arange(envNS-leftlen,envNS+rightlen):
            if envelope[ii]>envelope[ii-1] and envelope[ii]>envelope[ii+1] and envelope[ii]>ampmaxS:
                ampmaxS    = envelope[ii]
                maxidxS   = ii
        if (ampmaxS == 0.0):
            return 0.0
        ampmaxScS   = 0
        for ii in np.arange(envNScS-leftlen,envNScS+rightlen):
            if cut[ii]>cut[ii-1] and cut[ii]>cut[ii+1] and cut[ii]>ampmaxScS:
                ampmaxScS    = cut[ii]
                maxidxScS = ii
        if (ampmaxScS == 0):
            return 0.0
        iil      = 1
        while( envelope[maxidxS-iil]>0.5*ampmaxS and cut[maxidxScS-iil]>0.5*ampmaxScS ):
            iil = iil + 1
        iir      = 1
        while( envelope[maxidxS+iir]>0.5*ampmaxS and cut[maxidxScS+iir]>0.5*ampmaxScS ):
            iir = iir + 1
        cutS     = envelope[maxidxS-iil:maxidxS+iir]/envelope[maxidxS]
        cutScS   = cut[maxidxScS-iil:maxidxScS+iir]/cut[maxidxScS]
#    coef        = np.corrcoef(np.vstack((cutS,cutScS)))[0,1]
    coef        = 1/(np.sum((cutS-cutScS)**2*cutS)/np.sum(cutS))
#    print(tstar,coef)
#    plt.plot(np.arange(len(cutS))*tr.stats.delta,cutS,label='S')
#    plt.plot(np.arange(len(cutScS))*tr.stats.delta,cutScS,label='ScS')
#    plt.xlabel('Time (s)',size=14)
#    plt.legend()
#    plt.show()
#    plt.close()

    return coef

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
                  
    ## ---- exclude large perturbation
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
    eventids = ['usp000fsnz']
    eventmag = [6.7]
 ##   eventids = events['id']
 ##   eventmag = events['mag']
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

    ## ---- loop events
    for l,ids in enumerate(eventids):
        path_sacfile = '/home/meichen/work1/Pacific_attn/event_{}/waveforms/SAC_files'.format(ids)
        st = obspy.read('{}/*.BHT.SAC.raw'.format(path_sacfile))

        ## ---- loop traces
        for tr in st:
            ## ---- lowpass
            tr.filter('bandpass',freqmin=0.001,freqmax=2,zerophase=True)

            ## ---- check if nan value exist
            if (not np.all(tr.data)):
                continue

            ## ---- cut waveforms
            envNS,envNScS,flag = envelopemax(tr,model)
            if flag == 0:
                continue

            ## ---- signal noise ratio
            r1,r2,r3,r4,flag,leftlen,rightlen = s2n(tr,envNS,envNScS,model)
            if r1<3 or r2<6 or r3<2 or r4<4 or flag == 0:
                continue

            ## ---- calculate original instantaneous frequency
            fs,fscs = instan_freq(tr,0,envNS,envNScS,leftlen,rightlen)

            ## ---- find reflection point of ScS at CMB
            arrivals = model.get_pierce_points(tr.stats.sac.evdp,tr.stats.sac.gcarc,phase_list=["ScS"])
            for i in np.arange(arrivals[0].pierce.shape[0]):
                if(arrivals[0].pierce[i][3] == 2891.0):
                    reflgcarc  = arrivals[0].pierce[i][2]*tr.stats.sac.gcarc
            if (i == arrivals[0].pierce.shape[0]):
                print("Reflection point wrong")
                continue
            reflat,reflon = distaz2latlon(reflgcarc,tr.stats.sac.az,tr.stats.sac.evla,tr.stats.sac.evlo)

            ## ---- matching tstar
            coefmax = 0
            for ttstar in np.arange(0,10,0.1):
                coef = cc(tr,ttstar,envNS,envNScS,leftlen,rightlen,"S")
                if coef > coefmax:
                    tstar = ttstar
                    coefmax = coef
            for ttstar in np.arange(0,10,0.1):
                coef = cc(tr,ttstar,envNS,envNScS,leftlen,rightlen,"ScS")
                if coef > coefmax:
                    tstar = -ttstar
                    coefmax = coef
            if coefmax<0.98:
                tstar = 0.0
            # save results
            print("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}".format(tr.stats.sac.evla,tr.stats.sac.evlo,tr.stats.sac.evdp,tr.stats.sac.stla,tr.stats.sac.stlo,reflat,reflon,tstar,tr.stats.sac.gcarc))
            dtstar.append(tstar)
            evla.append(tr.stats.sac.evla)
            evlo.append(tr.stats.sac.evlo)
            evdp.append(tr.stats.sac.evdp)
            stla.append(tr.stats.sac.stla)
            stlo.append(tr.stats.sac.stlo)
            rfla.append(reflat)
            rflo.append(reflon)

#    np.savez("output",dtstar=dtstar,evla=evla,evlo=evlo,evdp=evdp,stla=stla,stlo=stlo,rfla=rfla,rflo=rflo)

main()
