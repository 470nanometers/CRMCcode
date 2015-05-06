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




#####main

###### get file, gather data ######
if len(sys.argv)>1:
  #filename = raw_input("Please enter full path name of the HepMC file: ")
  filename = str(sys.argv[1])
  eventarray=[]
  fileout = filename.strip('.hepmc')+'.txt'
  HEPfile=open(filename, 'r')# cant use loadtxt from numpy due to lines being formatted differently
  #output=open(fileout,'w+')
  Im= 0
  iCS= 0
  CS= 0
  CSe= 0
  numP= 0
  numQGP = 0
  numpi = 0
  numK = 0
  nummu = 0
  px=0
  py=0
  pz=0
  iterator =0
  for line in HEPfile:
    iterator +=1
    #print iterator
    if iterator >4: #dumb workaround to get past the empty line and the first E line so 1st event is not zeros
      #print 'test'
      parts = line.split() #P,V,E, U, C, H,F
      #print parts
      if parts[0]=='E': #start of new event reset, but where to record the info?
        #output.write(str(numP)+' '+str(numQGP)+' '+str(Im)+' '+str(CSe)+' '+str(iCS)+' '+str(CS)+'\n')
        eventarray.append(Event(CS,CSe,Im,iCS,numP,numQGP,numpi,numK,nummu)) 
        Im= 0 #reset
        iCS= 0
        CS= 0
        CSe= 0
        numP= 0
        numQGP = 0
        numpi = 0
        numK = 0
        nummu = 0
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
          numP+=1
          if parts[2] == '91':
            numQGP+=1
          if parts[2] == '211' or  parts[2] =='111' or  parts[2] =='-211': #pion codes for pi0pi+pi-
            numpi+=1
          if parts[2] == '321' or  parts[2] =='311' or  parts[2] =='-321':
            numK+=1
          if parts[2] == '13' or  parts[2] =='-13':
            print 'mu'
            nummu+=1
      if parts[0] == 'HepMC::IO_GenEvent-END_EVENT_LISTING': #record last event
        #output.write(str(numP)+' '+str(numQGP)+' '+str(Im)+' '+str(CSe)+' '+str(iCS)+' '+str(CS)+'\n')
        eventarray.append(Event(CS,CSe,Im,iCS,numP,numQGP,numpi,numK,nummu))
  HEPfile.close()
  #output.close()

### now to graphing #####
# what to graph?  want numQGP but against what?, just histogram?
  QGPhist=[]
  numhist=[]
  pct=[]
  kpi=[]
  mu=[]
  k=[]
  pii=[]
  for item in eventarray:
    if item.numP!=0 and item.numpi !=0:
      QGPhist.append(item.numQGP)
      numhist.append(item.numP)
      qp=float(item.numQGP)
      pp = float(item.numP)
      pct.append(qp/pp)
      kk=float(item.numK)
      pi = float(item.numpi)
      k.append(kk)
      pii.append(pi)
      kpi.append(kk/pi)
      m=float(item.nummu)
      mu.append(m)

  pdf =PdfPages(filename.strip('.hepmc')+'_U.pdf')
  xlab="Number of unique QGP-like events"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(QGPhist, bins=50, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Number of unique QGP-like events"
  ylab="Number of all unique particles"
  tle = filename.strip('.hepmc')+'_QGPvsP'
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(QGPhist, numhist)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()
  
  xlab="Number of unique QGP-like events"
  ylab="Fraction of QGP-like"
  tle = filename.strip('.hepmc')+'_QGPFrac'
  fig3 = figure()
  ax3 = fig3.add_subplot(111)
  plt.scatter(QGPhist, pct)
  #ylim(0,0.0005)
  #xlim(0,1000)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()
    
  xlab="Fraction of unique QGP-like events"
  ylab="Counts"
  tle = filename.strip('.hepmc')+'_QGPfrac_hist'
  fig4 = figure()
  ax4 = fig4.add_subplot(111)
  plt.hist(pct, bins=50, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Number of unique QGP-like events"
  ylab="K/pi ratio"
  tle = filename.strip('.hepmc')+'_Kpi'
  fig3 = figure()
  ax3 = fig3.add_subplot(111)
  plt.scatter(QGPhist, kpi)
  #ylim(0,0.0005)
  #xlim(0,100)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Number of unique QGP-like events"
  ylab="number of K"
  tle = filename.strip('.hepmc')+'_K'
  fig3 = figure()
  ax3 = fig3.add_subplot(111)
  plt.scatter(QGPhist, k)
  #ylim(0,0.0005)
  #xlim(0,100)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Number of unique QGP-like events"
  ylab="number of pi"
  tle = filename.strip('.hepmc')+'_pi'
  fig3 = figure()
  ax3 = fig3.add_subplot(111)
  plt.scatter(QGPhist, pii)
  #ylim(0,0.0005)
  #xlim(0,100)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Number of unique QGP-like events"
  ylab="Number of mu"
  tle = filename.strip('.hepmc')+'_mu'
  fig3 = figure()
  ax3 = fig3.add_subplot(111)
  plt.scatter(QGPhist, mu)
  #ylim(0,0.0005)
  #xlim(0,100)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()


  pdf.close()







#### if the file is not input ####
else:
  print "Usage: python CRMCreader.py [PATH/FILE]"
  print "File must be in HepMC format"

  




