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
from scipy import stats
import itertools
from mpl_toolkits.mplot3d import Axes3D
import ROOT as root #pyROOT, works on bradbury;  to use elsewhere (force)recompile root from APE with python enabled in the root config file

class Event:
  def __init__(self, numP,numQGP, numQGPnu, particlelist): #numpi,numK,CS, CSe,Im,iCS,
    #self.CS = CS #cross section
    #self.CSe = CSe #cross section error
    #self.Im = Im #impact parameter
    #self.iCS = iCS #inelastic Cross section
    self.numP = numP #number of unique particles per event
    self.numQGP = numQGP #number of unique QGP-like particles per event
    self.numQGPD = numQGPnu #number of non unique QGP-like particles aka number of QGP daughters per event
    self.numsBaryon = sBaryon
    self.numsMeson = sMeson
    self.particlelist = particlelist

class Particle:
  def __init__(self, px, py, pz, ID, energy, phi, pEnergy):
    self.px = px
    self.py = py
    self.pz = pz
    self.pT = np.sqrt(float(px)**2 + float(py)**2)
    self.eta=(0.5*np.log((np.sqrt(float(px)**2+float(py)**2+float(pz)**2)+float(pz))/(np.sqrt(float(px)**2+float(py)**2+float(pz)**2)-float(pz))))
    self.ID = int(ID)
    self.xF = float(pz)/float(energy) #xF=2*pL/P(max) = 2*pL/sqrt(sNN) , so this will; energy here is total energy
    #self.phix = np.arccos(float(px)/np.sqrt(float(px)**2 + float(py)**2)) #in radians
    #self.phiy = np.arcsin(float(py)/np.sqrt(float(px)**2 + float(py)**2)) #phix and phiy should match, but they don't.  Why?
    self.phi = np.arctan2(float(py),float(px))
    self.pEnergy = pEnergy #energy of particle
    ## as far as I can tell, pz is along the 'beam' i.e. direction of impact.
    ## pT = sqrt(px^2 + py^2)
    ## eta = 1/2 * ln ((P + pL)/(P - pL)) ;  pseudorapidity; pL==pz



#####main

###### get file, gather data ######
if len(sys.argv)>1:
  filename = str(sys.argv[1])
  eventarray=[]
  fileout = filename.strip('.hepmc')+'.txt'
  getenergy = fileout.split('_')
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
  px = 0
  py = 0
  pz = 0
  phi = 0
  pEnergy = 0
  particlelist=[]
  iterator = 0
  for line in HEPfile:
    iterator +=1
    if iterator >4: #dumb workaround to get past the empty line and the first E line so 1st event is not zeros
      parts = line.split() #types of lines prefixed with: P,V,E, U, C, H,F
      if parts[0]=='E': #start of new event -> reset
        eventarray.append(Event(numP,numQGP,numQGPnu, particlelist)) #numpi,numK,CS,CSe,Im,iCS,
        #Im= 0 #reset
        #iCS= 0
        #CS= 0
        #CSe= 0
        numP= 0
        numQGP = 0
        numQGPnu = 0
        sBaryon = 0
        px = 0
        py = 0
        pz = 0
        phi = 0
        pEnergy = 0
        particlelist=[]
      if parts[0] == 'H': #get impact parameter and inelastic cross section
        Im = parts[10]
        #ecc = float(parts[12])
        iCS = parts[13]
      if parts[0] == 'C': #get cross section and error
        CS = parts[1]
        CSe = parts[2]
      if parts[0] == 'P':
        if parts[2] == '91': #get non-unique count of QGP-like aka number of QGP daughters roughly
          numQGPnu+=1
        if px!= parts[3] or py!= parts[4] or pz != parts[5]: #equal px,py,pz==same particle, only do the following if not same particle
          px = parts[3]
          py = parts[4]
          pz = parts[5]
          pEnergy = parts[6]
          phi = parts[10]
          #particlearray.append(Particle(px, py, pz, parts[2], energy, phi, pEnergy))
          particlelist.append(Particle(px, py, pz, parts[2], energy, phi, pEnergy))
          numP+=1 #count number of unique particles
          if parts[2] == '91':
            numQGP+=1#get unique QGP count aka don't count multiple decays

      if parts[0] == 'HepMC::IO_GenEvent-END_EVENT_LISTING': #record last event
        eventarray.append(Event(numP,numQGP, numQGPnu, particlelist))#numpi,numK,CS,CSe,Im,iCS,
  HEPfile.close()

### now to graphing #####


  #pdf =PdfPages(filename.strip('.hepmc')+'_events.pdf') #put ALL graphs into one PDF

  ##Particle number graphs
  #QGPhist=[]
  #daughtershist=[]
  #dratio=[]
  #numhist=[]
  #pct=[]
  #mesons=[]
  #baryons=[]
  ## if there is a better way to do this I don't know it.  Wish I did.
  histlist=[]
  hists=[]
  length=[]
  lentot=0
  i =0
  for item in eventarray:
    i+=1
    if item.numP!=0 :#and item.numpi !=0:
      QGP=item.numQGP
      part=item.numP
      #qp=float(item.numQGP)
      #pp = float(item.numP)
      #pct.append(qp/pp)
      #daughtershist.append(item.numQGPD)
      #if item.numQGP !=0:
        #dratio.append(float(item.numQGPD)/float(item.numQGP))
      #else:
        #dratio.append(-1)

      #pT=[]
      eta=[]
      #phi=[]
      #ID=[]
      #xF=[]


      if part >0: #multiplicity cuts
        for item2 in item.particlelist: #particle lists for each event
          #pT.append(float(item2.pT))
          if str(item2.eta) != 'inf' and str(item2.eta) != '-inf':
            eta.append(item2.eta) #skip infinities, i.e. very far forward/backwards or errors
          #xF.append(item2.xF)
          #ID.append(item2.ID)
          #if item2.ID >95: #not photons, EM, or QGP/strings.  Note: muons are ID +- 13, but we'd see late muons from decays not in early shower
            #print item2.phi # I think this works

        '''
        xlab="pT"
        ylab="ID"
        tle = 'Event_'+str(i)+'_QGP='+str(QGP)+'_NumP='+str(part)
        fig = figure()
        ax = fig.add_subplot(111)
        plt.scatter(pT, ID)
        ylim(-3500,3500)
        #ax.set_yscale('log')
        xlabel(xlab)
        ylabel(ylab)
        title(tle)
        plt.savefig(pdf,format='pdf')
        #plt.show() 
        close()
        ''' 

        #pT=[]
        ID=[]
        #xF=[]
  
        #s = np.array(list(itertools.permutations(eta,2))) #repeats AB, AC, AD, BA, BC, BD, CA, CB, CD, etc.
        s = np.array(list(itertools.combinations(eta,2))) #no repeats : AB, AC, AD, BC, BD, CD, etc.
        etadiff = root.TH1F('etadiff_'+str(i)+'_QGP='+str(QGP)+'_N='+str(part), 'etadiff', 200, -6, 6)
        for item in s:
          etadiff.Fill(item[0] - item[1]) #slow filling, is there a faster way?, can't seem to find one
        hists.append(etadiff)
        histlist.append(eta)
        length.append(len(s))
        #lentot+=len(s)
  
        eta=[] #dump eta, we don't need it anymore
        '''
        xlab="Eta1-Eta2"
        ylab="Counts"
        tle = 'Event_'+str(i)+'_Kurtosis='+str(kurt)
        fig = figure()
        ax = fig.add_subplot(111)
        
        #plt.imshow(s[:,0] - s[:,1], cmap = 'jet')#, bins=100)
        #plt.scatter(s[:,0], s[:,1])
        plt.hist(s[:,0] - s[:,1], bins=200, range=(-5,5), log=True, histtype='step')#, normed=True)#)
        #plt.hexbin(s[:,0], s[:,1],cmap=cm.jet,bins='log', gridsize=300) 
        #plt.hexbin(X, Y,cmap=cm.jet,bins='log', gridsize=300) 
        #ax.plot_surface(X, Y, X-Y,cmap=cm.jet)
        #plt.colorbar()
        xlim(-5,5)
        #ylim(-5,20)
        xlabel(xlab)
        ylabel(ylab)
        title(tle)
        plt.savefig(pdf,format='pdf')
        #plt.show() 
        close()
        '''


  #print histlist
  #pdf.close()

  ### Now to get the mixed event background i.e. eta(event1)1-eta(event2)2'...etc.; mixed events shouldn't be correlated
  ## this will take forever, but will hopefully not eat all the memory
  totalmixed = root.TH1F('totalmixed', 'etadiff', 200, -6, 6)
  f=0
  while f < len(histlist): #get first list of etas
    j=f+1 #start on second eta list
    while j < len(histlist): #get second list of etas; repeat for all etas, no repeats i.e. 1-2, 1-3, 1-4, but no 1-1 or 4-1
      s = np.array(list(itertools.product(histlist[f],histlist[j]))) #get permutations of one event etas with all the other event's etas
      #print s
      for item in s:
        totalmixed.Fill(item[0]-item[1]) # get deltas and fill histogram, only one hist needed
      lentot+=len(s) #number of mixed event pairs
      j+=1
    f+=1
  #print totalmixed




  i=0
  c = root.TCanvas()
  c.Print(filename.strip('.hepmc')+'_eventsROOT.pdf[') #start pdf
  for item in hists: #graph all events
    coeff = lentot/length[i] #total mixed-event pairs/ one event pairs
    #print coeff, lentot, length[i]
    #graphme = root.TH1F('graph', 'etadiff', 200, -6, 6)
    graphme = item 
    graphme.Divide(totalmixed) # signal pairs hist/ mixed events hist
    graphme.Scale(coeff) #multiply by num mixed pairs/num signal pairs
    root.gROOT.Reset()
    root.gROOT.SetStyle("Plain")
    c.Modified()
    #root.gStyle.SetOptStat()
    #root.gStyle.SetOptFit()
    root.gStyle.SetPalette(1)
    #graphme.Fit("gaus")  
    graphme.Draw()
    graphme.GetXaxis().SetTitle("etadiff")
    graphme.GetYaxis().SetTitle("C(deltaEta)")
    graphme.SetLabelSize(0.025);
    graphme.SetLabelSize(0.025,"Y");
    #c.WaitPrimitive()
    c.Print(filename.strip('.hepmc')+'_eventsROOT.pdf') #fill multipage pdf
    i+=1

  c.Print(filename.strip('.hepmc')+'_eventsROOT.pdf]') #close pdf
  



#### if the file is not input ####
else:
  print "Usage: python CRMCreader-strange.py [PATH/FILE]"
  print "File must be in HepMC format"

  




