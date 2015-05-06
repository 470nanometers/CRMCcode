#!/usr/bin/python

import sys
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import stats
import itertools
from mpl_toolkits.mplot3d import Axes3D
#import ROOT as root #pyROOT, works on bradbury;  to use elsewhere (force)recompile root from APE with python enabled in the root config file
from ROOT import *



def set_palette(name="palette", ncontours = 99):#ncontours=999):

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
    elif name == "hot":
        stops = [ 0.00, 0.25, 0.50, 0.75, 1.00]
        red   = [ 0.00, 0.50, 1.00, 1.00, 1.00]
        green = [ 0.00, 0.00, 0.55, 1.00, 1.00]
        blue  = [ 0.00, 0.00, 0.00, 0.00, 1.00]
    elif name== "invhot":
        stops = [ 0.00, 0.25, 0.50, 0.75, 1.00]
        red   = [ 1.00, 1.00, 1.00, 0.50, 0.00]
        green = [ 1.00, 1.00, 0.55, 0.00, 0.00]
        blue  = [ 1.00, 0.00, 0.00, 0.00, 0.00]



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




if len(sys.argv)>1:
  filename = str(sys.argv[1])
  eventarray=[]
  #HEPfile=open(filename, 'r')# cant use loadtxt from numpy due to lines being formatted differently
  f = TFile(filename, 'read')


#need to import root file, then graph various combinations of variables
  c = TCanvas()
  #c.Print(filename.strip('.hepmc')+'_ROOTetaphimix.pdf[') #start pdf
  c.Print(filename.strip('.root')+'_ROOT.pdf[') #start pdf
  filename2=filename.strip('.root')+'_ROOT.pdf' #fill multipage pdf

  #use && to join multiple cuts, use || for or

  ### Graphing 1D histos: x,y,r,px,py,pt,id,gen,theta,thetaP histos
  gPad.SetLogz()

  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  gPad.SetLogy()
  pNtuple.SetLineWidth(2)
  pNtuple.SetLineColor(4)
  pNtuple.Draw("abs(x)","","")
  pNtuple.SetLineColor(2)
  pNtuple.Draw("abs(y)","","same")
  pNtuple.SetLineColor(1)
  pNtuple.Draw("sqrt(x*x+y*y)","","same")
  #GetXaxis().SetTitle("Distance(m)")
  ##pNtuple.GetYaxis().SetTitle("N")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  gPad.SetLogy()
  pNtuple.SetLineWidth(2)
  pNtuple.SetLineColor(4)
  pNtuple.Draw("abs(px)","","")
  pNtuple.SetLineColor(2)
  pNtuple.Draw("abs(py)","","same")
  pNtuple.SetLineColor(1)
  pNtuple.Draw("sqrt(px*px+py*py)","","same")
  ##pNtuple.GetXaxis().SetTitle("Momentum(GeV)")
  ##pNtuple.GetYaxis().SetTitle("N")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)


  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  gPad.SetLogy()
  pNtuple.SetLineColor(4)
  pNtuple.SetLineWidth(2)
  pNtuple.Draw("gen>>h(100)","","")
  ##pNtuple.GetXaxis().SetTitle("Generation")
  ##pNtuple.GetYaxis().SetTitle("N")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  gPad.SetLogy()
  pNtuple.SetLineColor(4)
  pNtuple.SetLineWidth(2)
  pNtuple.Draw("id>>h(100)","","")
  ##pNtuple.GetXaxis().SetTitle("ID")
  ##pNtuple.GetYaxis().SetTitle("N")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  gPad.SetLogy(0)
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  gPad.SetLogy()
  pNtuple.SetLineColor(4)
  pNtuple.SetLineWidth(2)
  pNtuple.Draw("TMath::ATan2(y,x)>>h(100,-3.14,3.14)","","")

  #pNtuple.GetXaxis().SetTitle("#theta(rad)")
  #pNtuple.GetYaxis().SetTitle("N")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  gPad.SetLogy()
  pNtuple.SetLineColor(4)
  pNtuple.SetLineWidth(2)
  pNtuple.Draw("TMath::ATan2(py,px)>>h(100,-3.14,3.14)","","")
  #pNtuple.GetXaxis().SetTitle("#theta(rad)")
  #pNtuple.GetYaxis().SetTitle("N")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  gPad.SetLogy(0)
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetOptStat()
  #gPad.SetLogy()
  pNtuple.Draw("id:sqrt(x*x+y*y)>>h2(100, 0, 7500, 80, 0, 80)","","colz")
  #pNtuple.GetXaxis().SetTitle("Corsika ID")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)
  ## 1D histos done

  ##Graphing 2D histos:
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetNumberContours(99)

  pNtuple.Draw("id:sqrt(x*x+y*y)>>h2(100,0,3000,80,0,80)","","colz") ## id:radius, want Ncontour higher
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Corsika ID")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  pNtuple.Draw("gen:sqrt(x*x+y*y)","","colz") ## gen:radius
  #pNtuple.GetYaxis().SetTitle("Generation")
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  pNtuple.Draw("gen:id","","colz") ## gen:id
  #pNtuple.GetXaxis().SetTitle("Corsika ID")
  #pNtuple.GetYaxis().SetTitle("Generation")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  #pNtuple.Draw("px:x>>h2(100,-5000,5000,100,-3,3)","","colz") ## px:x
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Momentum(GeV)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  #c.Print(filename2)

  #pNtuple.Draw("py:y>>h2(100,-5000,5000,100,-3,3)","","colz") ## py:y
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Momentum(GeV)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  #c.Print(filename2)

  #pNtuple.Draw("sqrt(px*px+py*py):sqrt(x*x+y*y)>>h2(100,0,5000,100,0,3)","","colz") ## pT:r
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Momentum(GeV)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  #c.Print(filename2)

  #pNtuple.Draw("sqrt(px*px+py*py+pz*pz):sqrt(x*x+y*y)>>h2(100,0,1000,100,0,500)","","colz")## P:r, may be useless
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Momentum(GeV)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  #c.Print(filename2)

  pNtuple.Draw("py:px>>h2(100,-3,3,100,-3,3)","","colz") ## py:px, may be useless, may want to change binning and range
  #pNtuple.GetXaxis().SetTitle("Momentum(GeV)")
  #pNtuple.GetYaxis().SetTitle("Momentum(GeV)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  pNtuple.Draw("py:px>>h2(200,-0.01,0.01,200,-0.01,0.01)","","colz") ## py:px, zoom in on central feature
  #pNtuple.GetXaxis().SetTitle("Momentum(GeV)")
  #pNtuple.GetYaxis().SetTitle("Momentum(GeV)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  pNtuple.Draw("y:x>>h2(100,-5000,5000,100,-5000,5000)","","colz") ## y:x, may want to change binning and range ~1000-1500 meters
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  pNtuple.Draw("y:x>>h2(200,-2000,2000,200,-2000,2000)","","colz") ## y:x, smaller binning and smaller range
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  #######Polar graphs, need an axis display, still cant get "pol colz" to work
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gStyle.SetNumberContours(99)
  pNtuple.Draw("sqrt(x*x+y*y):TMath::ATan2(y,x)>>h2(36,-3.14,3.14,30,0,3000)","","pollego20z") 
  #pNtuple.Draw("sqrt(x*x+y*y):TMath::ATan2(y,x)>>h2(36,-3.14,3.14,50,0,2000)","","pollego20") ##theta:R polar, adjust binning and radial distance
  gPad.SetTheta(90)
  gPad.SetLogz()
  #pNtuple.Draw("sqrt(x*x+y*y):TMath::ATan2(y,x)>>h2(36,-3.14,3.14,50,0,100)","","pollego20")
##theta:R polar, adjust binning and radial distance
  ## interesting note: if log y is set, radius is log scale
  c.Print(filename2)

  #########
  ##2D histos done


  ##Graphing 3D histos
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  #pNtuple.Draw("id:y:x>>h3(100,-2000,2000,100,-2000,2000,80,0,80)","","iso") #interesting but no idea how to intrepret
  gPad.SetLogz(0)
  pNtuple.Draw("y:x:id","","colz") ## works but overwhelmed by muons everywhere (the mu with extra info)
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)
  ##3D histos done

  set_palette('wtf')
  ## x:y by particle type??
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gPad.SetLogz()
  gPad.SetTheta(90)
  pNtuple.Draw("y:x>>h(50,-500,500,50,-500,500)","id==7||id==8||id==9","colz") # pi, binsize == 10m
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)
  pNtuple.Draw("y:x>>h(150,-1500,1500,150,-1500,1500)","id==5||id==6||id==75||id==76","colz")# mu
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)
  pNtuple.Draw("y:x>>h(150,-1500,1500,150,-1500,1500)","id==3||id==2","colz")# electron/positron
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)
  pNtuple.Draw("y:x>>h(150,-1500,1500,150,-1500,1500)","id==10||id==11||id==12||id==16","colz") # K
  #pNtuple.GetXaxis().SetTitle("Distance(m)")
  #pNtuple.GetYaxis().SetTitle("Distance(m)")
  #pNtuple.SetLabelSize(0.025)
  #pNtuple.SetLabelSize(0.025,"Y")
  c.Print(filename2)

  ## radial, need to show axis for radius
  gROOT.Reset()
  gROOT.SetStyle("Plain")
  c.Modified()
  gPad.SetLogz()
  gPad.SetTheta(90)
  gStyle.SetOptStat()
  pNtuple.Draw("sqrt(x*x+y*y):TMath::ATan2(y,x)>>h(36,-3.14,3.14,20,0,200)","id==9||id==7||id==8&&sqrt(x*x+y*y)<200","pollego20z") #pi, zoomed in
  c.Print(filename2)
  pNtuple.Draw("sqrt(x*x+y*y):TMath::ATan2(y,x)>>h(36,-3.14,3.14,15,0,1500)","id==5||id==6||id==75||id==76&&sqrt(x*x+y*y)<1500","pollego20z") #mu, needs logz
  c.Print(filename2)
  ## K too sparse
  #set_palette('invhot')
  pNtuple.Draw("sqrt(x*x+y*y):TMath::ATan2(y,x)>>h(36,-3.14,3.14,20,0,2000)","id==3||id==2&&sqrt(x*x+y*y)<2000","pollego20z")#EM
  c.Print(filename2)




  #c.Print(filename.strip('.hepmc')+'_ROOTetaphimix.pdf]') #close pdf
  c.Print(filename.strip('.root')+'_ROOT.pdf]') #close pdf
  
#maybe have select options for what graphs to do?, do this later
  

#### if the file is not input ####
else:
  print "Usage: python ROOTgraphing.py [PATH/FILE]"
  print "File must be in .root format"

  




