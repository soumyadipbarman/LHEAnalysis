// g++ LHEreader LHEreader.C LHEF.h -L/home/soumyadip/govoni/ `root-config --cflags --libs`
// ./LHEreader events.lhe
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <assert.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"

#include "conversions.h"
#include "vbf.h"
#include <cmath>
#include "LHEF.h"

using namespace std ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char** argv) 
{
  if (argc != 2)
    {
      std::cout << ">>> Usage:   " << argv[0] << " LHE.lhe" << std::endl;
      return -1;
    }

  ifstream ifs_wH (argv[1]) ;
  LHEF::Reader reader (ifs_wH) ;

  TH1F h_mH ("h_mH", "h_mH", 400, 0, 400) ;
  
  int ieve = 0 ;
  //PG loop over events
  while ( reader.readEvent ()) 
    {
      ++ieve;
      if (ieve % 1000 == 0) std::cerr << "event " << ieve << "\n" ;
  
      vector<lorentzVector> Hs ;
      
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
          lorentzVector particle 
            (
              reader.hepeup.PUP.at (iPart).at (0), //PG px
              reader.hepeup.PUP.at (iPart).at (1), //PG py
              reader.hepeup.PUP.at (iPart).at (2), //PG pz
              reader.hepeup.PUP.at (iPart).at (3) //PG E
            ) ;
          //PG leptons
          if (abs (reader.hepeup.IDUP.at (iPart)) == 6)  //PG W
            {
              Hs.push_back (particle) ;
//              cerr << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                   << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                   << "\n" ;
            }

          if (Hs.size () > 0) h_mH.Fill (sqrt (Hs.at (0).M2 ())) ;
	  //if (Hs.size () > 0) h_mH.Fill (Hs.at (0).perp2()) ;
        } //PG loop over particles in the event
    } //PG loop over events

  cout << "read " << ieve << " events\n" ;
  
  TFile outfile ("lhereader.root", "recreate") ;
  h_mH.Write () ;
  outfile.Close () ;
  
  return 0 ;
}

