//Need this to read the crmc HepMC outputs and report the number and kind of particles produced for each event

//Particles are in the format:
// P int      int  double  double  double  double  double         int   double  double        int            int          int      int
// P barcode PDGid   px      py      pz    energy  generatedMass status PolTheta PolPhi barcodeForVertex EntriesInFlow CodeIndex CodeinIndex

// A particle with the same id, px, py, pz and energy is the same particle, but has multiple outputs from the vertex

//Vertex is labeled:
// V int    int double double double double     int           int         int          double
// V barcode id   x      y       z     ctau  numberOrphanIn numberOut EntriesinWeight ListWeight

//New events are headed by lines starting with: E, U, C, H, F
// E - general event info
// U - Momentum and position units
// C - Cross-section info
// H - Heavy Ion information
// F - Pdf Info (partons)
// E eventNumber NumberOfInteractions EventScale alphaQCD alphaQED SignalProcessID barcodeSPvertex NumberVerticies barcodeParticle1 barcodeParticle2 numberEntriesRandomState ListRandomState numberWeightEntries WeightEntries
// U  MomentumUnits(MeV/GeV) LengthUnits(mm/cm)
// C crossSection(pb) error(pb)
// H see pg 21 of HepMC user manual
// F see pg 21 of HepMC user manual

// Will want for each event: entries 10 and 13 from H; all of C; number of unique PDG91 from P and count of all unique P for each event

//make object:
// Object name is a count of the event (i.e object 1, 2, 3, ...,N), initialize when starting code?
// has CrossSection, CSerror, Impact parameter, NN inelastic CrossSection, Number of unique P produced, Number of unique(?) 91(QGP-like)

// argument against unique 91 events, want to know how large of an effect it has on particles, more particles == more NON-unique events of 91
// easiest: count number of P and number of 91P


#include <cstdio>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<sstream>
#include<math.h>
//#include<TROOT.h>
//#include<TStyle.h>
//#include<TCanvas.h>
//#include<TTree.h>
//#include<TFile.h>
//#include<TH2F.h>
//#include<TMath.h>
//#include<TApplication.h>
#include <sstream>
#include <numeric>
#include <regex>

using namespace std;

class Event
{
  private:
    double CrossSection;
    double CSError;
    double Impact;
    double inelasticCS;
    double numParticles; //non-unique
    double numQGPlike; //non-unique



  public:
    Event(double CS, double CSe, double Im, double iCS, double numP, double num91)
      {
      CrossSection = CS;
      CSError = CSe;
      Impact = Im;
      inelasticCS = iCS;
      numParticles = numP;
      numQGPlike = num91;
      }


    double GetCS()
      {
      return CrossSection;
      }
    double GetCSe()
      {
      return CSError;
      }
    double GetIm()
      {
      return Impact;
      }
    double GetiCS()
      {
      return inelasticCS;
      }
    double GetnumPart()
      {
      return numParticles;
      }
    double GetQGP()
      {
      return numQGPlike;
      }

};

int main(int argc, char *argv[])
  {  
  char* filename;
  if(argc == 2 )
    {
      filename = argv[1];
    }
  else
    {
      //cout << "Takes 1 arguments, " << argc-1 << " given." << endl;
      cout << "USAGE: ./CRMCreader [/PATH/FILENAME]" << endl;
      return -1;
    }
  vector<Event> events;
  ifstream file;
  file.open(filename, ios::in);
  file.getline();
  







  return 0;
  };

