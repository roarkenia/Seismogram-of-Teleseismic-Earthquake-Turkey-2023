# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 19:43:49 2024

@author: lenovo
"""
#PREVIOSLY INSTALATION 
# pip install haversine
# conda install conda-forge::obspy
# conda install pygmt

dir().clear()

#Libreries
import pygmt
from obspy import read
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.inventory import read_inventory
from matplotlib import pyplot as plt
from IPython import get_ipython
import numpy as np
import haversine as hs   
from haversine import haversine, Unit
get_ipython().run_line_magic('matplotlib', 'qt')

###########################################################################
###########################################################################
#################    TURKEY EARTHQUAKE M7.8    ############################
###########################################################################
###########################################################################

#Azimuthal projection
fig = pygmt.Figure()
#fig.basemap(region=[0,360,-90,0], projection="S 0/-90/5i", frame="g")
fig.coast(projection="E 37.22/37.01/180/7i", region="g", land="olivedrab",water="royalblue3")
fig.coast(shorelines=True)
fig.basemap(region="g", projection="S 0/-90/7i", frame=["ag10", "+tPROYECCIÓN AZIMUTAL EQUIDISTANTE"])
# fig.basemap(frame="a")

client= Client("IRIS")
network="IU"
location="10"
channel="LHZ"

Earthquakelocation=[37.014,37.226]
starttime=UTCDateTime("2023-02-06T01:17:34")
endtime=UTCDateTime(starttime+5000)
inv= client.get_stations(network=network, channel=channel, starttime= starttime, endtime=endtime)

latitude=[]
longitude=[]
stations=[]
y=[]

for k in range(len(inv[0].stations)):
    stations.append(inv[0].stations[k].code)
    latitude.append(inv[0].stations[k].latitude)
    longitude.append(inv[0].stations[k].longitude)
    y.append(latitude[k]+0.7)
    

# fig.plot(x=longitude, y=latitude, style="i0.3c",fill="indianred",region="g",projection="E 37.22/37.01/180/7i")
# font = "7p,Helvetica-Bold"
# fig.text(x=longitude, y=y, text=stations, font=font)
# fig.show()

stationnames=["ANTO", "KIEV", "PAB", "KMBO", "MA2", "ADK", "NWAO", "POHA", "PMSA", "PTCN"]

stationselect=[]
latitudeselect=[]
longitudeselect=[]

for k in range(len(stationnames)):
    for j in range(len(inv[0].stations)):
        if stations[j]==stationnames[k]:
            stationselect.append(stations[j])
            latitudeselect.append(latitude[j])
            longitudeselect.append(longitude[j])
            
distdegree=[]
distkilometers=[]
# aconstant=[]
# cconstant=[]

for k in range(len(stationnames)):
    Location2=[latitudeselect[k],longitudeselect[k]]
    distdegree.append(round(hs.haversine(Earthquakelocation,Location2,unit=Unit.DEGREES),1))
    distkilometers.append(round(hs.haversine(Earthquakelocation,Location2,unit=Unit.KILOMETERS),2))
#         dlon = -(Earthquakelocation[0]*np.pi/180)+(longitudeselect[k]*np.pi/180)
#         dlat = -(Earthquakelocation[1]*np.pi/180)+(latitudeselect[k]*np.pi/180)
#         aconstant.append((np.sin(dlat/2)**2 + np.cos(Earthquakelocation[0]) * np.cos(latitudeselect[k]) * np.sin(dlon / 2)**2))
#         cconstant.append(2*np.arctan2(np.sqrt(aconstant[k]),np.sqrt(1-aconstant[k]))*180/np.pi)


pygmt.makecpt(cmap="plasma", series=[min(longitude), max(longitude)])           
fig.plot(x=longitudeselect, y=latitudeselect, style="i0.8c",fill=latitudeselect,region="g",projection="E 37.22/37.01/180/7i", cmap=True)
font = "12p,Helvetica-Bold"
fig.text(x=np.array(longitudeselect)+0.1*(1-np.array(distdegree)/max(distdegree))+1, y=np.array(latitudeselect)+0.1*(1-np.array(distdegree)/max(distdegree))+1, text=stationselect, font=font)
fig.show()

stt2=[]
inv2=[]  
tr=[]
trP=[]
pre_filt = (0.01, 0.02, 40.0, 45.0)
Parrives=[88,205,407,453,663,784,852]
Sarrives=[161,422,720,846,1232,1463,1727]
Pkparrives=[1214]

fig, ax=plt.subplots(10,1)
fig2, ax2=plt.subplots(4,2,constrained_layout = True)
for k in range(len(stationnames)):
    stt2.append(client.get_waveforms(network=network, station=stationselect[k], location=location, channel=channel, starttime=starttime, endtime=endtime))
    inv2.append(client.get_stations(network=network, station= stationselect[k], location=location, channel=channel, level= "response", starttime= starttime, endtime=endtime))
    stt2[k].attach_response(inv2[k])
    stt2[k].remove_response(output="DISP", inventory=inv2[k], pre_filt = pre_filt)
    stt2[k].detrend('linear')
    stt2[k].detrend('demean')
    stt2[k].taper(type='cosine',max_percentage=0.05)
    tr.append(stt2[k].traces[0])
    tr[k].filter("lowpass", freq=0.45)
    temporalcopy=stt2[k].copy()
    ax[k].set_xlim(0, 7000)
    ax[k].set_ylim(-1, 1)
    ax[k].spines['top'].set_visible(False)
    ax[k].spines['right'].set_visible(False)
    ax[k].spines['bottom'].set_visible(False)
    ax[k].spines['left'].set_visible(False)
    ax[k].get_xaxis().set_ticks([])
    ax[k].get_yaxis().set_ticks([])
    ax[k].set_ylabel(stationselect[k])
    ax[k].plot(tr[k].times(), tr[k].data/max(tr[k].data),color="darkorange")
    ax[k].annotate(str(distdegree[k])+"° "+str(distkilometers[k])+"km" ,xy=(5500,0.5) , xytext=(0, 2), textcoords='offset points', ha='center', va='bottom')
    
    if distdegree[k]<110:
        temporalcopy.trim(starttime+Parrives[k]-20,starttime+Parrives[k]+80)
        trP.append(temporalcopy.traces[0])
        trP[k].filter("lowpass", freq=0.1)
        if k<=3:
            i=0
            j=k 
        if k>3:
            i=1
            j=k-4

        ax2[j,i].plot(trP[k].times(),trP[k].data)
        ax2[j,i].annotate("Estación "+(stationnames[k]) ,xy=(50,0) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        ax[k].axvline(Parrives[k], ymin=0.5, ymax=0.8, label="Fase P", color="black",linestyle='-')
        ax[k].axvline(Sarrives[k], ymin=0.5, ymax=0.8, label="Fase S", color="black", linestyle='-')
        ax[k].annotate("P" ,xy=(Parrives[k],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        ax[k].annotate("S" ,xy=(Sarrives[k],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")

    if distdegree[k]>150:
        ax[k].axvline(Pkparrives[0], ymin=0.5, ymax=0.8, label="Fase PkP", color="black",linestyle='-')
        ax[k].annotate("PkP" ,xy=(Pkparrives[0],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        
        temporalcopy.trim(starttime+Pkparrives[0]-20,starttime+Pkparrives[0]+80)
        trP.append(temporalcopy.traces[0])
        trP[k].filter("lowpass", freq=0.1)
        ax2[3,1].plot(trP[k].times(),trP[k].data)
        ax2[3,1].annotate("Estación "+(stationnames[k]) ,xy=(50,0) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        
    if distdegree[k]<150 and distdegree[k]>110:
        ax[k].annotate("Zona de sombra" ,xy=(Pkparrives[0],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        trP.append(temporalcopy.traces[0])
        trP[k].filter("lowpass", freq=0.1)

Trace1=tr.copy()        
fig, ax=plt.subplots(5,2,constrained_layout = True)
for k in range(len(stationnames)):
    data = np.fft.rfft(tr[k].data)
    freqs = np.fft.rfftfreq(tr[k].stats.npts, d=tr[k].stats.delta)
    if k<=4:
        i=0
        j=k 
    if k>4:
        i=1
        j=k-5
        
    ax[j,i].set_xlim(0.005, 0.5)
    ax[j,i].set_yscale('symlog')
    ax[j,i].set_xlabel("Frequencia (Hz)")
    ax[j,i].plot(freqs, 2*np.abs(data)**2)
    ax[j,i].annotate("Estación "+(stationnames[k]), xy=(0.35,0.00000001) , xytext=(0.35, 0.00000001), textcoords='offset points', ha='center', va='bottom')
    fig.tight_layout(pad=0.1)
    # ax.text(x=4500, y=0, text=distdegree[k], font=font)
    # ax.text(x=5500, y=0, text=distkilometers[k], font=font)

#stt= client.get_waveforms(network="IU", station="ANTO", location="10", channel=channel, starttime=starttime, endtime=endtime)
#stt[0].plot()

###########################################################################
###########################################################################
#################    TURKEY EARTHQUAKE M7.5    ############################
###########################################################################
###########################################################################

#Azimuthal projection
fig = pygmt.Figure()
#fig.basemap(region=[0,360,-90,0], projection="S 0/-90/5i", frame="g")
fig.coast(projection="E 37.22/37.01/180/7i", region="g", land="olivedrab",water="royalblue3")
fig.coast(shorelines=True)
fig.basemap(region="g", projection="S 0/-90/7i", frame=["ag10", "+tPROYECCIÓN AZIMUTAL EQUIDISTANTE"])
# fig.basemap(frame="a")

client= Client("IRIS")
network="IU"
location="10"
channel="LHZ"

Earthquakelocation=[38.011,37.196]
starttime=UTCDateTime("2023-02-06T10:24:48")
endtime=UTCDateTime(starttime+5000)
inv= client.get_stations(network=network, channel=channel, starttime= starttime, endtime=endtime)

distdegree=[]
distkilometers=[]
# aconstant=[]
# cconstant=[]

for k in range(len(stationnames)):
    Location2=[latitudeselect[k],longitudeselect[k]]
    distdegree.append(round(hs.haversine(Earthquakelocation,Location2,unit=Unit.DEGREES),1))
    distkilometers.append(round(hs.haversine(Earthquakelocation,Location2,unit=Unit.KILOMETERS),2))
#         dlon = -(Earthquakelocation[0]*np.pi/180)+(longitudeselect[k]*np.pi/180)
#         dlat = -(Earthquakelocation[1]*np.pi/180)+(latitudeselect[k]*np.pi/180)
#         aconstant.append((np.sin(dlat/2)**2 + np.cos(Earthquakelocation[0]) * np.cos(latitudeselect[k]) * np.sin(dlon / 2)**2))
#         cconstant.append(2*np.arctan2(np.sqrt(aconstant[k]),np.sqrt(1-aconstant[k]))*180/np.pi)


stt2=[]
inv2=[]  
tr=[]
trP=[]
pre_filt = (0.01, 0.02, 40.0, 45.0)
Parrives=[68,168,386,436,642,763,832]
Sarrives=[122,365,706,819,1212,1436,1615]
Pkparrives=[1206]

fig, ax=plt.subplots(10,1)
fig2, ax2=plt.subplots(4,2,constrained_layout = True)
for k in range(len(stationnames)):
    stt2.append(client.get_waveforms(network=network, station=stationselect[k], location=location, channel=channel, starttime=starttime, endtime=endtime))
    inv2.append(client.get_stations(network=network, station= stationselect[k], location=location, channel=channel, level= "response", starttime= starttime, endtime=endtime))
    stt2[k].attach_response(inv2[k])
    stt2[k].remove_response(output="DISP", inventory=inv2[k], pre_filt = pre_filt)
    stt2[k].detrend('linear')
    stt2[k].detrend('demean')
    stt2[k].taper(type='cosine',max_percentage=0.05)
    tr.append(stt2[k].traces[0])
    tr[k].filter("lowpass", freq=0.45)
    temporalcopy=stt2[k].copy()
    ax[k].set_xlim(0, 7000)
    ax[k].set_ylim(-1, 1)
    ax[k].spines['top'].set_visible(False)
    ax[k].spines['right'].set_visible(False)
    ax[k].spines['bottom'].set_visible(False)
    ax[k].spines['left'].set_visible(False)
    ax[k].get_xaxis().set_ticks([])
    ax[k].get_yaxis().set_ticks([])
    ax[k].set_ylabel(stationselect[k])
    ax[k].plot(tr[k].times(), tr[k].data/max(tr[k].data),color="darkorange")
    ax[k].annotate(str(distdegree[k])+"° "+str(distkilometers[k])+"km" ,xy=(5500,0.5) , xytext=(0, 2), textcoords='offset points', ha='center', va='bottom')
    
    if distdegree[k]<110:
        temporalcopy.trim(starttime+Parrives[k]-20,starttime+Parrives[k]+80)
        trP.append(temporalcopy.traces[0])
        trP[k].filter("lowpass", freq=0.1)
        if k<=3:
            i=0
            j=k 
        if k>3:
            i=1
            j=k-4

        ax2[j,i].plot(trP[k].times(),trP[k].data)
        ax2[j,i].annotate("Estación "+(stationnames[k]) ,xy=(50,0) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        ax[k].axvline(Parrives[k], ymin=0.5, ymax=0.8, label="Fase P", color="black",linestyle='-')
        ax[k].axvline(Sarrives[k], ymin=0.5, ymax=0.8, label="Fase S", color="black", linestyle='-')
        ax[k].annotate("P" ,xy=(Parrives[k],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        ax[k].annotate("S" ,xy=(Sarrives[k],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        
    if distdegree[k]>150:
        ax[k].axvline(Pkparrives[0], ymin=0.5, ymax=0.8, label="Fase PkP", color="black",linestyle='-')
        ax[k].annotate("PkP" ,xy=(Pkparrives[0],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        
        temporalcopy.trim(starttime+Pkparrives[0]-20,starttime+Pkparrives[0]+80)
        trP.append(temporalcopy.traces[0])
        trP[k].filter("lowpass", freq=0.1)
        ax2[3,1].plot(trP[k].times(),trP[k].data)
        ax2[3,1].annotate("Estación "+(stationnames[k]) ,xy=(50,0) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        
    if distdegree[k]<150 and distdegree[k]>110:
        ax[k].annotate("Zona de sombra" ,xy=(Pkparrives[0],-1) , xytext=(0, 0), textcoords='offset points', ha='center', va='bottom', color="black")
        trP.append(temporalcopy.traces[0])
        trP[k].filter("lowpass", freq=0.1)
        
Trace2=tr.copy()           
fig, ax=plt.subplots(5,2,constrained_layout = True)
for k in range(len(stationnames)):
    data = np.fft.rfft(tr[k].data)
    freqs = np.fft.rfftfreq(tr[k].stats.npts, d=tr[k].stats.delta)
    if k<=4:
        i=0
        j=k 
    if k>4:
        i=1
        j=k-5
        
    ax[j,i].set_xlim(0.005, 0.5)
    ax[j,i].set_yscale('symlog')
    ax[j,i].set_xlabel("Frequencia (Hz)")
    ax[j,i].plot(freqs, 2*np.abs(data)**2)
    ax[j,i].annotate("Estación "+(stationnames[k]), xy=(0.35,0.00000001) , xytext=(0.35, 0.00000001), textcoords='offset points', ha='center', va='bottom')
    fig.tight_layout(pad=0.1)
    # ax.text(x=4500, y=0, text=distdegree[k], font=font)
    # ax.text(x=5500, y=0, text=distkilometers[k], font=font)

#stt= client.get_waveforms(network="IU", station="ANTO", location="10", channel=channel, starttime=starttime, endtime=endtime)
#stt[0].plot()