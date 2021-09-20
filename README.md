This project is to investigate the attenuation factor "tstar" in the mantle.

Computing global attenuation profile is not the aim of this project. Since our understanding of the attenuation in the mantle is poorly known, distinct attenuation factors along different corridors is prominent to demonstrate the lateral variation of the mantle attenuation. Here we use deep earthquakes (>100 km) to avoid crustal attenuation from the source side. We choose clusters of stations in Alaska, USArray, and Europe, hoping that the receiver-side crustal attenuation are faily averaged out by dense stations.

data_download.py includes following processes: 
    Use obspy.mass_downloader for downloading a large amount of data. Note: a java converter is needed.
    The downloaded miniseed files are then converted to SAC files.
    Event info are not included in miniseed files to minimize data size, so we add them back
    Remove instrument response and rotate to gcarc
    Remove unnecessary files

SR.py:
    The spectral ratio method to calculate tstar between S waves of two seismograms. It utilizes obspy to manipulate seismic traces in a event-wise way.

IFM.py:
    The instantaneous frequency matching method to calculate tstar.

events_2000_2020_6.5:
    A csv file includes informations of earthquakes occurred in 2000-2020 with magnitude larger than Mw6.5. Downloaded from USGS.
