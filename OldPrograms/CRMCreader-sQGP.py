#!/usr/bin/python

#############################
#Need this to read the crmc HepMC outputs and report the number and kind of particles produced for each event

#Events are in the format:
# E    int    int     double    double    double    int     int     int     int     int     int     long    int     double
# E  eventnum numinteractions eventscale alphaQCD alphaQED signalprocessID barcodeSignalProcessV numVerticies barBeam1 barBeam2 numRandState listRandState numWeight listWeight 
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

#TODO:
## for starters, count number of strange baryons vs mesons and total strange particles (as well as total strangeness possibly, should be zero), separate strange vs antistrange

## secondary, find if the strange particles are originating from the mini QGPs (ID 91)

##################
import matplotlib as mlib
import sys
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import itertools

class Event:
  def __init__(self, numP,numQGP, numQGPnu, sBaryon, sMeson): #numpi,numK,CS, CSe,Im,iCS,
    #self.CS = CS #cross section
    #self.CSe = CSe #cross section error
    #self.Im = Im #impact parameter
    #self.iCS = iCS #inelastic Cross section
    self.numP = numP #number of unique particles per event
    self.numQGP = numQGP #number of unique QGP-like particles per event
    self.numQGPD = numQGPnu #number of non unique QGP-like particles aka number of QGP daughters per event
    #self.numpi = numpi
    #self.numK = numK
    self.numsBaryon = sBaryon
    self.numsMeson = sMeson
    #self.ecc = ecc #eccentricity of participating nucleons, doesn't show any change in any situation, not correct ecc

class Particle:
  def __init__(self, px, py, pz, ID, energy, phi, pEnergy):
    self.px = px
    self.py = py
    self.pz = pz
    self.pT = np.sqrt(float(px)**2 + float(py)**2)
    self.eta=(0.5*np.log((np.sqrt(float(px)**2+float(py)**2+float(pz)**2)+float(pz))/(np.sqrt(float(px)**2+float(py)**2+float(pz)**2)-float(pz))))
    self.ID = int(ID)
    self.xF = float(pz)/float(energy) #xF=pL/P(max) , so this will; energy here is total energy
    self.phi = phi
    self.pEnergy = pEnergy #energy of particle
    ## as far as I can tell, pz is along the 'beam' i.e. direction of impact.
    ## pT = sqrt(px^2 + py^2)
    ## eta = 1/2 * ln ((P + pL)/(P - pL)) ;  pseudorapidity; pL==pz



#####main

###### get file, gather data ######
if len(sys.argv)>1:
  #filename = raw_input("Please enter full path name of the HepMC file: ")
  filename = str(sys.argv[1])
  eventarray=[]
  particlearray=[]
  fileout = filename.strip('.hepmc')+'.txt'
  getenergy = fileout.split('_')
  #print getenergy
  energy = float(getenergy[5].strip('.txt'))
  HEPfile=open(filename, 'r')# cant use loadtxt from numpy due to lines being formatted differently
  ##initialize variables
  #Im= 0
  #iCS= 0
  #CS= 0
  #CSe= 0
  numP= 0
  numQGP = 0
  numQGPnu = 0
  sBaryon = 0
  sMeson = 0
  #numpi = 0
  #numK = 0
  #ecc = 0
  px = 0
  py = 0
  pz = 0
  phi = 0
  pEnergy = 0
  iterator = 0
  for line in HEPfile:
    iterator +=1
    if iterator >4: #dumb workaround to get past the empty line and the first E line so 1st event is not zeros
      parts = line.split() #types of lines prefixed with: P,V,E, U, C, H,F
      if parts[0]=='E': #start of new event -> reset
        eventarray.append(Event(numP,numQGP,numQGPnu, sBaryon, sMeson)) #numpi,numK,CS,CSe,Im,iCS,
        #Im= 0 #reset
        #iCS= 0
        #CS= 0
        #CSe= 0
        numP= 0
        numQGP = 0
        numQGPnu = 0
        sBaryon = 0
        sMeson = 0
        #numpi = 0
        #numK = 0
        #ecc = 0
        px = 0
        py = 0
        pz = 0
        phi = 0
        pEnergy = 0
      if parts[0] == 'H': #get impact parameter and inelastic cross section
        Im = parts[10]
        #ecc = float(parts[12])
        iCS = parts[13]
      if parts[0] == 'C': #get cross section and error
        CS = parts[1]
        CSe = parts[2]
      if parts[0] == 'P':
        if parts[2] == '91': #get non-unique count of QGP-like
          numQGPnu+=1
        if px!= parts[3] or py!= parts[4] or pz != parts[5]: #equal px,py,pz==same particle, only do the following if not same particle
          px = parts[3]
          py = parts[4]
          pz = parts[5]
          pEnergy = parts[6]
          phi = parts[10]
          particlearray.append(Particle(px, py, pz, parts[2], energy, phi, pEnergy))
          numP+=1 #count number of unique particles
          if parts[2] == '91':
            numQGP+=1#get unique QGP count aka don't count multiple decays
          if abs(int(parts[2]))==3122 or abs(int(parts[2]))==3222 or abs(int(parts[2]))==3212 or abs(int(parts[2]))==3112 or abs(int(parts[2]))==3322 or abs(int(parts[2]))==3312 or abs(int(parts[2]))==3334: #strange and multistrange baryons, abs(x) because ID is negative for anti-particles. Lambda, Sigma(-/0/+), Xi(0/-), Omega(-)
            sBaryon+=1#strange baryons
          if abs(int(parts[2]))==310 or abs(int(parts[2]))==311 or abs(int(parts[2]))==321: #abs(item.ID)==130(0L) strange mesons aka Kaons(0L/0S/0/+) abs(x) because ID is negative for anti-particles
            sMeson+=1#strange mesons

      if parts[0] == 'HepMC::IO_GenEvent-END_EVENT_LISTING': #record last event
        eventarray.append(Event(numP,numQGP, numQGPnu, sBaryon, sMeson))#numpi,numK,CS,CSe,Im,iCS,
  HEPfile.close()

### now to graphing #####


  pdf =PdfPages(filename.strip('.hepmc')+'_strange.pdf') #put ALL graphs into one PDF

  ##Particle number graphs
  QGPhist=[]
  daughtershist=[]
  dratio=[]
  numhist=[]
  pct=[]
  mesons=[]
  baryons=[]
  #ecc=[]
  #kpi=[]
  #k=[]
  #pii=[]
  ## if there is a better way to do this I don't know it.  Wish I did.
  for item in eventarray:
    if item.numP!=0 :#and item.numpi !=0:
      QGPhist.append(item.numQGP)
      numhist.append(item.numP)
      qp=float(item.numQGP)
      pp = float(item.numP)
      pct.append(qp/pp)
      mesons.append(item.numsMeson)
      baryons.append(item.numsBaryon)
      #ecc.append(item.ecc)
      #pi = float(item.numpi)
      #k.append(kk)
      #pii.append(pi)
      #kpi.append(kk/pi)
      daughtershist.append(item.numQGPD)
      if item.numQGP !=0:
        dratio.append(float(item.numQGPD)/float(item.numQGP))
      else:
        dratio.append(-1)




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
  tle = filename.strip('.hepmc')
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
  ylab="Number of strange mesons"
  tle = filename.strip('.hepmc')
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(QGPhist, mesons)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Number of unique QGP-like events"
  ylab="Number of strange baryons"
  tle = filename.strip('.hepmc')
  fig2 = figure()
  ax2 = fig2.add_subplot(111)
  plt.scatter(QGPhist, baryons)
  #ylim(0,3000)
  #xlim(0,500)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()
  
  xlab="Number of unique QGP-like events"
  ylab="numQGP/totalParticles"
  tle = filename.strip('.hepmc')
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
    
  xlab="numQGP/totalParticles"
  ylab="Counts"
  tle = filename.strip('.hepmc')
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
  ylab="Number of QGP daughters"
  tle = filename.strip('.hepmc')
  fig3 = figure()
  ax3 = fig3.add_subplot(111)
  plt.scatter(QGPhist, daughtershist)
  #ylim(0,0.0005)
  #xlim(0,100)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Ratio of number of daughters to unique QGP-like events"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(dratio, bins=50, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  #pdf.close()

  ##Lets empty things to free up memory
  QGPhist=[]
  daughtershist=[]
  dratio=[]
  numhist=[]
  pct=[]
  mesons=[]
  baryons=[]
  eventarray=[] #we are done with it


  ##Momentum graphs
  x=[]
  y=[]
  z=[]
  pT=[]
  spT=[] #pT of strange Baryons/mesons
  sxF=[]
  eta=[]
  seta=[]
  #deta=[]
  #ID=[]
  xF=[]
  sBaryon=[]
  asBaryon=[]
  sMeson=[]
  asMeson=[]
  phi = []
  for item in particlearray:
    #x.append(float(item.px))
    #y.append(float(item.py))
    #z.append(float(item.pz))
    #print item.eta
    #phi.append(item.phi)#phi == 0 for some reason in teh polarization data from HepMC
    pT.append(float(item.pT))
    #ID.append(item.ID)
    if str(item.eta) == 'inf' :
      eta.append(50) #if inf, is a divide by zero, meaing pz >>1, meaning very far forward)
    elif str(item.eta) =='-inf':
      eta.append(-15)
    else:
      eta.append(item.eta)
    xF.append(item.xF)
    if abs(item.ID)==3122 or abs(item.ID)==3222 or abs(item.ID)==3212 or abs(item.ID)==3112 or abs(item.ID)==3322 or abs(item.ID)==3312 or abs(item.ID)==3334: #strange and multistrange baryons, abs(x) because ID is negative for anti-particles. Lambda, Sigma(-/0/+), Xi(0/-), Omega(-)
      #print item.ID
      spT.append(item.pT)
      sxF.append(item.xF)
      if str(item.eta) == 'inf' :
        seta.append(50) #if inf, is a divide by zero, meaing pz >>1, meaning very far forward)
      elif str(item.eta) =='-inf':
        seta.append(-15)
      else:
        seta.append(item.eta)
      if item.ID>0:
        sBaryon.append(item.ID)#strange baryons
      else:
        asBaryon.append(abs(item.ID))#antistrange
    if abs(item.ID)==310 or abs(item.ID)==311 or abs(item.ID)==321: #abs(item.ID)==130(0L) strange mesons aka Kaons(0L/0S/0/+) abs(x) because ID is negative for anti-particles
      #spT.append(item.pT)
      #sxF.append(item.xF)
      if item.ID>0:
        sMeson.append(item.ID)#strange mesons
      else:
        #print item.ID
        asMeson.append(abs(item.ID))#antistrange
  #print eta

#TODO:  Need eta1 vs eta2 and phi1 vs phi2  OR a deltaEta = eta1-eta2 histo

  #this is dumb and will take forever to do:
  '''
  f=0
  while f < len(eta):
    d=0
    while d < len(eta):
      if d !=f: #so not subtracting the same 'particle' from itself, need pairs
        deta.append(eta[f] - eta[d])
      d+=1
    f+=1
  '''
  # ^^^BAD PLAN!!!! ABORT! ABORT!  EATS ALL MEMORY!


  #empty particle array as we are done with it:
  particlearray=[]

  #pdf =PdfPages(filename.strip('.hepmc')+'QGPmomentums.pdf')
  '''
  xlab="Number of unique QGP-like events by px"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(x, bins=50, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  plt.show() 
  close()


  xlab="Number of unique QGP-like events by py"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(y, bins=50, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  plt.show() 
  close()

  xlab="Number of unique QGP-like events by pz"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(z, bins=50, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  plt.show() 
  close()
  '''

  xlab="pT of strange baryons"# and mesons"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(spT, bins=50, align='left')
  ax.set_yscale('log')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="xF of strange baryons"# and mesons"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(sxF, bins=50, align='left')
  ax.set_yscale('log')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Pseudorapidity of strange baryons"# and mesons"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(seta, bins=50, align='left')
  #ax.set_yscale('log')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Particle ID of strange baryons"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(sBaryon, bins=75, align='left')
  plt.hist(asBaryon, bins=75, align='right', color='red')#shifted one to the right for better viewing
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="Particle ID of Kaon species"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(sMeson, bins=50, align='left')
  plt.hist(asMeson, bins=50, align='right', color='red')#width is not showing up
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  
  xlab="Pseudorapidity (eta)"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(eta, bins=100, align='left')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="pT"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(pT, bins=100, align='left')
  #ylim(1,1000000)
  #ax.set_yscale('log')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="xF"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(xF, bins=100, align='left')
  #ylim(1,1000000)
  #ax.set_yscale('log')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  #Empty everything else just in case:
  x=[]
  y=[]
  z=[]
  pT=[]
  spT=[] #pT of strange Baryons/mesons
  sxF=[]
  seta=[]
  #ID=[]
  xF=[]
  sBaryon=[]
  asBaryon=[]
  sMeson=[]
  asMeson=[]
  phi = []

  #etacorr = np.correlate(eta, eta, mode='same') #auto correlate eta with eta?, unnormalized, what are x and y axes?
  #etadiff = [x-y for x in eta for y in eta if y <20 and y >-5 and x <20 and x >-5] #NOPE!
  #eta =[]
  #print len(eta), len(etacorr)
  #etadiff => want pairwise difference
  #s=[eta,eta] #NOPE!!!
  #s=list(itertools.product(eta,eta)) #this does the right thing BUT EATS ALL MEMORY YET AGAIN!
  #print s

  xlab="Eta1"
  ylab="Eta2"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.acorr(eta, usevlines=False)#, bins=100) #autocorrelation
  #plt.xcorr(eta, eta)#, bins=100) #cross correlation
  #plt.scatter(eta, eta,c=etacorr, cmap='jet')#this is meaningless really
  #plt.hist(s)
  #xlim(-5,20)
  #ylim(-5,20)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  '''
  xlab="eta1"
  ylab="eta2"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.scatter(eta,eta) #this is derp and doesn't work
  ylim(2, 10)
  xlim(2, 10)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()

  xlab="phi"
  ylab="phi"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.scatter(phi,phi,c=z)
  #ylim(0, 5000)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()


  xlab="eta"
  ylab="Particle ID #"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.scatter(eta,ID)
  ylim(0, 5000)
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()


  xlab="Particle ID #"
  ylab="Counts"
  tle = filename.strip('.hepmc')
  fig = figure()
  ax = fig.add_subplot(111)
  plt.hist(ID,range =(-500, 500), bins=1000, align='left')
  #ylim(0.00001,1000000)
  #ax.set_yscale('log')
  xlabel(xlab)
  ylabel(ylab)
  title(tle)
  plt.savefig(pdf,format='pdf')
  #plt.show() 
  close()
  '''

  pdf.close()





#### if the file is not input ####
else:
  print "Usage: python CRMCreader-strange.py [PATH/FILE]"
  print "File must be in HepMC format"

  



