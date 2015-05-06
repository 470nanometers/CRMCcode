#!/usr/bin/python

#Need this to read the crmc HepMC outputs and report the number and kind of particles produced for each event

#Particles are in the format:
# P int      int  double  double  double  double  double         int   double  double        int            int          int      int
# P barcode PDGid   px      py      pz    energy  generatedMass status PolTheta PolPhi barcodeForVertex EntriesInFlow CodeIndex CodeinIndex

# A particle with the same id, px, py, pz and energy is the same particle, but has multiple outputs from the vertex

#Vertex is labeled:
# V int    int double double double double     int           int         int          double
# V barcode id   x      y       z     ctau  numberOrphanIn numberOut EntriesinWeight ListWeight

#New events are headed by lines starting with: E, U, C, H, F
# E - general event info
# U - Momentum and position units
# C - Cross-section info
# H - Heavy Ion information
# F - Pdf Info (partons)
# E eventNumber NumberOfInteractions EventScale alphaQCD alphaQED SignalProcessID barcodeSPvertex NumberVerticies barcodeParticle1 barcodeParticle2 numberEntriesRandomState ListRandomState numberWeightEntries WeightEntries
# U  MomentumUnits(MeV/GeV) LengthUnits(mm/cm)
# C crossSection(pb) error(pb)
# H see pg 21 of HepMC user manual
# F see pg 21 of HepMC user manual

# Will want for each event: entries 10 and 13 from H; all of C; number of unique PDG91 from P and count of all unique P for each event

#make object:
# Object name is a count of the event (i.e object 1, 2, 3, ...,N), initialize when starting code?
# has CrossSection, CSerror, Impact parameter, NN inelastic CrossSection, Number of unique P produced, Number of unique(?) 91(QGP-like)

# argument against unique 91 events, want to know how large of an effect it has on particles, more particles == more NON-unique events of 91
# easiest: count number of P and number of 91P
import matplotlib as mlib
import sys
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages


class Event:
  def __init__(self, CS, CSe,Im,iCS,numP,numQGP,numpi,numK,nummu):
    self.CS = CS #cross section
    self.CSe = CSe # cross section erroe
    self.Im = Im # impact parameter
    self.iCS = iCS #inelastic Cross section
    self.numP = numP#number of non unique particles
    self.numQGP = numQGP#number of non unique QGP-like particles
    self.numpi = numpi
    self.numK = numK
    self.nummu = nummu

class Sim: #take energy and averages,  1000 events per sim
  def __init__(self,energy,summu,sumpi,sumK,sumP,sumQGP,col,events):
    self.E = energy #inGeV
    self.aveP = sumP/float(events)#number of non unique particles
    self.aveQGP = sumQGP/float(events)#number of non unique QGP-like particles
    self.avepi = sumpi/float(events)
    self.aveK = sumK/float(events)
    self.avemu = summu/float(events)
    self.color = col
    #still need errorbars



#####main

###### get file, gather data ######
if len(sys.argv)>1:
  i=1
  simarray=[]
  while i < len(sys.argv): #get multiple files
    filename = str(sys.argv[i])
    #eventarray=[]
    fileout = filename.strip('.hepmc')
    HEPfile=open(filename, 'r')# cant use loadtxt from numpy due to lines being formatted differently
    #output=open(fileout,'w+')
    getenergy = fileout.split('_')
    energy = float(getenergy[5])
    if getenergy[3]=='p':
      color = 'red'
    elif getenergy[3]=='C':
      color = 'green'
    else: #'pdg'
      color = 'blue'
    '''
    Im= 0
    iCS= 0
    CS= 0
    CSe= 0
    numP= 0
    numQGP = 0
    numpi = 0
    numK = 0
    nummu = 0
    '''
    events=1
    sumP= 0
    sumQGP = 0
    sumpi = 0
    sumK = 0
    summu = 0
    px=0
    py=0
    pz=0
    iterator =0
    for line in HEPfile:
      iterator +=1
#since we start in the first event
      if iterator >4: #dumb workaround to get past the empty line and the first E line so 1st event is not zeros
        parts = line.split() #P,V,E, U, C, H,F
        if parts[0]=='E': #start of new event reset, but where to record the info?
          #eventarray.append(Event(CS,CSe,Im,iCS,numP,numQGP,numpi,numK,nummu)) 
          events+=1 #count all events to scale averaging
          px=0
          py=0
          pz=0
        if parts[0] == 'H':
          Im = parts[10]
          iCS = parts[13]
        if parts[0] == 'C':
          CS = parts[1]
          CSe = parts[2]
        if parts[0] == 'P':
          if px!= parts[3] or py!= parts[4] or pz != parts[5]: #equal px,py,pz==same particle, only do if not same particle
            px=parts[3]
            py=parts[4]
            pz=parts[5]
            sumP+=1
            if parts[2] == '91':
              sumQGP+=1
            if parts[2] == '211':# or  parts[2] =='111' or  parts[2] =='-211': #pion codes for pi0pi+pi-
              sumpi+=1
            if parts[2] == '321':# or  parts[2] =='311' or  parts[2] =='-321':
              sumK+=1
            if parts[2] == '13' or  parts[2] =='-13':
              summu+=1
        #if parts[0] == 'HepMC::IO_GenEvent-END_EVENT_LISTING': #record last event
          #output.write(str(numP)+' '+str(numQGP)+' '+str(Im)+' '+str(CSe)+' '+str(iCS)+' '+str(CS)+'\n')
          #eventarray.append(Event(CS,CSe,Im,iCS,numP,numQGP,numpi,numK,nummu))
    simarray.append(Sim(energy,summu,sumpi,sumK,sumP,sumQGP,color,events))
    #print energy
    HEPfile.close()
    #output.close()
    i+=1

### now to graphing #####
# will need errorbars eventually
  energy=[]
  numhist=[]
  pct=[]
  kpi=[]
  QGP=[]
  kk=[]
  pii=[]
  col = []
  mm=[]
  for item in simarray:
    energy.append(item.E)
    pi = item.avepi
    k = item.aveK
    #print pi, k, item.aveQGP, item.E
    kpi.append(k/pi)
    kk.append(k)
    pii.append(pi)
    QGP.append(item.aveQGP)
    col.append(item.color)
    mm.append(item.avemu)
    

  pdf =PdfPages('X_C.pdf')


  xlab="Energy (GeV)"
  ylab="Average Counts unique QGP"
  tle = "QGP vs Energy"
  fig = figure()
  ax = fig.add_subplot(111)
  plt.scatter(energy,QGP,c=col)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  ax.set_xscale('log')
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Energy (GeV)"
  ylab="Average K/pi ratio"
  tle = "Energy vs K/pi ratio"
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(energy, kpi, c=col)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  ax2.set_xscale('log')
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Energy (GeV)"
  ylab="Average pi "
  tle = "Energy vs pi "
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(energy, pii, c=col)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  ax2.set_xscale('log')
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Energy (GeV)"
  ylab="Average K"
  tle = "Energy vs K"
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(energy, kk, c=col)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  ax2.set_xscale('log')
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Energy (GeV)"
  ylab="Number Mu"
  tle = "Energy vs Mu"
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(energy, mm, c=col)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  ax2.set_xscale('log')
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()
  

  pdf.close()







#### if the file is not input ####
else:
  print "Usage: python CRMCreader.py [PATH/FILE]"
  print "File must be in HepMC format"

  




