#include "LHEF.h"
#include "TH1.h"
#include "TFile.h"
#include "cmath"
#include "TLorentzVector.h"

// c++ -o testPowheg01 `root-config --glibs --cflags` testPowheg01.cpp 

using namespace std ;


TLorentzVector getPart (LHEF::Reader & reader ,int iPart)
{
  TLorentzVector part (
    reader.hepeup.PUP.at (iPart).at (0), //PG px
    reader.hepeup.PUP.at (iPart).at (1), //PG py
    reader.hepeup.PUP.at (iPart).at (2), //PG pz
    reader.hepeup.PUP.at (iPart).at (3) //PG E
  ) ;
  return part ;
}



int main (int argc, char **argv) 
{

  // Open a stream connected to an event file:
  if (argc < 2) exit (1) ;
  ifstream ifs(argv[1]);

  // Create the Reader object:
  LHEF::Reader reader(ifs);

  // Print out the header information:
  cout << reader.headerBlock;

  // Print out the addinional comments in the init block:
  cout << reader.initComments;

  // Print out the beam energies:
  cout << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << endl;

  // Now loop over all events:
  long ieve = 0;
  
  TH1F mH ("mH", "mH", 500, 0, 1000) ; 
  TH1F pTH ("pTH", "pTH", 500, 0, 1000) ; 
  TH1F mH_VBF ("mH_VBF", "mH_VBF", 500, 0, 1000) ; 
  TH1F pTH_VBF ("pTH_VBF", "pTH_VBF", 500, 0, 1000) ; 
  while ( reader.readEvent() ) 
    {

      ++ieve;
      if (ieve % 1000 == 0) cout << "event " << ieve << "\n" ;
  
      vector<int> leptons ;      
      vector<int> bosons ;      
      //PG loop over particles in the event

      TLorentzVector higgs ;
      TLorentzVector j1 ;
      TLorentzVector j2 ;
          
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
          if (reader.hepeup.IDUP.at (iPart) == 25) higgs = getPart (reader, iPart) ;
          if (reader.hepeup.IDUP.at (iPart) == 1) j1 = getPart (reader, iPart) ;
          if (reader.hepeup.IDUP.at (iPart) == 2) j2 = getPart (reader, iPart) ;
        } //PG loop over particles in the event

      mH.Fill (higgs.M ()) ;
      pTH.Fill (higgs.Pt ()) ;

      // VBF cuts
      if (j1.Pt () < 30 || j2.Pt () < 30) continue ;
      if (fabs (j1.Eta () - j2.Eta ()) < 3.5) continue ; 
      TLorentzVector diJet = j1 + j2 ;
      if (diJet.M () < 450) continue ;

      mH_VBF.Fill (higgs.M ()) ;
      pTH_VBF.Fill (higgs.Pt ()) ;

    } // Now loop over all events

  TFile histosFile ("output.root","recreate") ;
  histosFile.cd () ;
  mH.Write () ;
  pTH.Write () ;
  mH_VBF.Write () ;
  pTH_VBF.Write () ;
  histosFile.Close () ;

  return 0 ;
}

/*
number of particles: 6
	 part type [0] 2	 status -1
	 part type [1] 1	 status -1
	 part type [2] 25	 status 1
	 part type [3] 1	 status 1
	 part type [4] 2	 status 1
	 part type [5] 21	 status 1

*/
