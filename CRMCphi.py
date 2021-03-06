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

#H is heavy ion collision information
# H   int   int   int   int   int   int   int   int   int   float   float   float   float
# H hardScatterings projectileParticipants targetPatricipants NNcollisions spectatorNeutrons spectatorProtons N-Nwounded Nwounded-N Nwound-Nwound ImpactParameter AzimuthalAngle Eccentricity NNineasticCrossSection
#one or more of these shoudl tell us the centrality

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

###############################################################################################################
##
##
## This program will count the number of QGP, Particles from a .hepmc file 
## and calculate the Delta phi distribution between them
##      Written March 2014
##        By Danielle LaHurd
##
##
###############################################################################################################
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
import random
import CRMC2CORSIKA as C2C



def set_palette(name="palette", ncontours = 100):#ncontours=999):

    from array import array
    from ROOT import TColor, gStyle

    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
        # elif name == "whatever":
        # (define more palettes)
    elif name =='wtf':
        stops = [0.00, 0.25, 0.50, 0.75, 1.00]
        red   = [0.00, 0.00, 0.67, 1.00, 0.51]
        green = [0.00, 0.61, 0.80, 0.20, 0.00]
        blue  = [0.51, 0.50, 0.08, 0.00, 0.00]

    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)

    # For older ROOT versions
    #gStyle.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

class Event:
  def __init__(self, numP,numQGP, numQGPnu, particlelist, energy, parti): #numpi,numK,CS, CSe,Im,iCS,
    #self.CS = CS #cross section
    #self.CSe = CSe #cross section error
    #self.Im = Im #impact parameter
    #self.iCS = iCS #inelastic Cross section
    self.energy = energy
    self.numP = numP #number of unique particles per event
    self.numQGP = numQGP #number of unique QGP-like particles per event
    self.numQGPD = numQGPnu #number of non unique QGP-like particles aka number of QGP daughters per event
    self.numsBaryon = sBaryon
    self.numsMeson = sMeson
    self.particlelist = particlelist
    self.parti = parti

class Particle:
  def __init__(self, px, py, pz, ID, energy, phi, pEnergy, status):
    self.px = px
    self.py = py
    self.pz = pz
    self.pT = np.sqrt(float(px)**2 + float(py)**2)
    self.eta=(0.5*np.log((np.sqrt(float(px)**2+float(py)**2+float(pz)**2)+float(pz))/(np.sqrt(float(px)**2+float(py)**2+float(pz)**2)-float(pz))))
    self.ID = int(ID)
    self.xF = float(pz)/float(energy) #xF=2*pL/P(max) = 2*pL/sqrt(sNN) , so this will; energy here is total energy
    self.phi = np.arctan2(float(py),float(px))
    self.pEnergy = pEnergy #energy of particle
    self.status = status
    ## pz is along the 'beam' i.e. direction of impact.
    ## pT = sqrt(px^2 + py^2)
    ## eta = 1/2 * ln ((P + pL)/(P - pL)) ;  pseudorapidity; pL==pz
    #self.v2 = (float(px)**2-float(py)**2)/(float(px)**2+float(py)**2) #its actually supposed to be the expectation value, but this should give us some idea



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
  parti = 0
  sBaryon = 0
  sMeson = 0
  px = 0
  py = 0
  pz = 0
  phi = 0
  pEnergy = 0
  status =0
  particlelist=[]
  iterator = 0
  for line in HEPfile:
    iterator +=1
    if iterator >4: #dumb workaround to get past the empty line and the first E line so 1st event is not zeros
      parts = line.split() #types of lines prefixed with: P,V,E, U, C, H,F
      if parts[0]=='E': #start of new event -> reset
        eventarray.append(Event(numP,numQGP,numQGPnu, particlelist, energy, parti)) #numpi,numK,CS,CSe,Im,iCS,
        #Im= 0 #reset
        #iCS= 0
        #CS= 0
        #CSe= 0
        numP= 0
        numQGP = 0
        numQGPnu = 0
        parti = 0
        sBaryon = 0
        px = 0
        py = 0
        pz = 0
        phi = 0
        pEnergy = 0
        status = 0
        particlelist=[]
      if parts[0] == 'H': #get impact parameter and inelastic cross section
        Im = parts[10]
        iCS = parts[13]
        parti = int(parts[2]) + int(parts[3])
      if parts[0] == 'C': #get cross section and error
        CS = parts[1]
        CSe = parts[2]
      if parts[0] == 'P':
        if parts[2] == '91': #get non-unique count of QGP-like aka number of direct QGP daughters, roughly
          numQGPnu+=1
        if px!= parts[3] or py!= parts[4] or pz != parts[5]: #equal px,py,pz==same particle, only do the following if not same particle
          px = parts[3]
          py = parts[4]
          pz = parts[5]
          status = parts[8] #should be 1
          pEnergy = parts[6]
          phi = parts[10]
          if parts[2] == '91':
            numQGP+=1#get unique QGP count aka don't count multiple decays
          if str(status) == '1': # count only if final state particle
            particlelist.append(Particle(px, py, pz, parts[2], energy, phi, pEnergy, status))
            numP+=1 #count number of unique, final state, particles


      if parts[0] == 'HepMC::IO_GenEvent-END_EVENT_LISTING': #record last event
        eventarray.append(Event(numP,numQGP, numQGPnu, particlelist, energy, parti))#numpi,numK,CS,CSe,Im,iCS,
  HEPfile.close()

### now to graphing #####



  ##Particle graphs
  ## if there is a better way to do this I don't know it.  Wish I did.
  histlist=[]
  hists=[]
  length=[]
  lentot=0
  i =0
  for item in eventarray: #write output files here?
    i+=1
    if item.numP!=0 :
      QGP=item.numQGP
      part=item.numP
      parti = item.parti
      phi=[]

      if part > 800: #this is here for multiplicity cuts, most ridges appear above N=800, setting teh cut to ~ 750 should get all of them


        #stackout = open('/corsikaqgp/stacks/testing2/'+filename.strip('.hepmc')+'_'+str(i)+'.part', "w+")  
        #stackout.write("    "+str(part)+" "+ '%6e' %float(item.energy)+'\n') 
        #1st line, num particles, primary energy
        p=0
        for item2 in item.particlelist: #particle lists for each event
          p+=1
          #corid = C2C.convert(item2.ID)

          #if abs(item2.ID) ==13 or abs(item2.ID) ==111 or abs(item2.ID) ==211 or abs(item2.ID) ==130 or abs(item2.ID) ==321 or abs(item2.ID) ==310 or abs(item2.ID) ==3122 or abs(item2.ID) ==3212 or abs(item2.ID) ==3112 or abs(item2.ID) ==3222 or abs(item2.ID) ==3322 or abs(item2.ID) ==3312 or abs(item2.ID) ==3334:
          if item2.ID >95: #not photons, EM, or QGP/strings.  Note: muons are ID +- 13, but we'd see late muons from decays not in early shower
            #if str(item2.eta) != 'inf' and str(item2.eta) != '-inf':
            phi.append(item2.phi) # I think this works
  
        f = np.array(list(itertools.permutations(phi,2))) #repeats AB, AC, AD, BA, BC, BD, CA, CB, CD, etc.
        #f = np.array(list(itertools.combinations(etaphi,2))) #combination of permutations of eta/phi pairs
        etaphidiff = root.TH1F('phidiff_'+str(i)+'_QGP='+str(QGP)+'_N='+str(part), '#Delta#phi', 300, -4, 4)#used to be 100x100 bins, now 50x50
        for item in f:
          etaphidiff.Fill(item[0]-item[1]) #should be phi1-phi2 ([p1, p2]); hopefully this is the correct phi
        hists.append(etaphidiff)
        histlist.append(phi)
        length.append(len(f))
        #lentot+=len(s)

  ### Now to get the mixed event background i.e. eta(event1)1-eta(event2)2'...etc.; mixed events shouldn't be correlated
  ## this will take forever, but will hopefully not eat all the memory
  #totalmixed = root.TH1F('totalmixed', 'etadiff', 200, -6, 6)
  totalmixed = root.TH1F('totalmixed', '#Delta#phi', 300, -4, 4)#used to be 100x100 bins, now 50x50
  f=0
  while f < len(histlist): #get first list of etas
    j=f+1 #start on second eta list
    while j < len(histlist): #get second list of etas; repeat for all etas, no repeats i.e. 1-2, 1-3, 1-4, but no 1-1 or 4-1
      s = np.array(list(itertools.product(histlist[f],histlist[j]))) #get permutations of one event etas with all the other event's etas
      for item in s:
        totalmixed.Fill(item[0]-item[1]) # get deltas and fill histogram, only one hist needed
      lentot+=len(s) #number of mixed event pairs
      j+=1
    f+=1



  i=0
  fit=root.TF1("fit","1+2*[0]*cos(2*x)")
  #fit.SetLineColor(kRed)
  c = root.TCanvas()
  #c.Print(filename.strip('.hepmc')+'_ROOTetaphimix.pdf[') #start pdf
  c.Print(filename.strip('.hepmc')+'_ROOTphi2.pdf[') #start pdf
  for item in hists: #graph all events
    if length[i] ==0:
      coeff = 1 
    else:
      coeff = lentot/length[i] #total mixed-event pairs/ one event pairs
    graphme = item 
    graphme.Divide(totalmixed) # signal pairs hist/ mixed events hist, needs the coeff
    graphme.Scale(coeff) #multiply by num mixed pairs/num signal pairs
    root.gROOT.Reset()
    root.gROOT.SetStyle("Plain")
    c.Modified()
    root.gStyle.SetOptStat()
    root.gStyle.SetOptFit()
    #set_palette('wtf') #still not very visible around 2-3, should find some way to 'set' teh z axis scale
    #root.gStyle.SetPalette(30)
    #root.gPad.SetLogz() 
    #graphme.SetMinimum(0)#testing z axis min
    #graphme.SetMaximum(5)#testing z axis max
    graphme.Draw()
    graphme.Fit("fit")
    #graphme.Draw("contz")
    #graphme.SetMinimum(1e-15)
    graphme.GetXaxis().SetTitle("#Delta#phi (radians)")
    #graphme.GetYaxis().SetTitle("#Delta#phi (radians)")
    graphme.SetLabelSize(0.025);
    graphme.SetLabelSize(0.025,"Y");
    #c.WaitPrimitive() #uncomment if you want to look at the graph before its put in the pdf
    c.Print(filename.strip('.hepmc')+'_ROOTphi2.pdf') #fill multipage pdf
    i+=1

  c.Print(filename.strip('.hepmc')+'_ROOTphi2.pdf]') #close pdf
  

  

#### if the file is not input ####
else:
  print "Usage: python CRMCphi.py [PATH/FILE]"
  print "File must be in HepMC format"

  




