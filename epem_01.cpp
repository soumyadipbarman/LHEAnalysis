// c++ -o epem_01 `root-config --glibs --cflags` -lm epem_01.cpp
// ./epem_01 LHE_withHiggs.lhe LHE_withoutHiggs.lhe
// ./epem_01 ../VBF_epem_withHiggs/Events/run_01/events.lhe ../VBF_epem/Events/run_01/events.lhe 

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

struct Histos
{
  TH1F h_deta_ee ;
  TH1F h_m_ee ;
  TH1F h_mWW ;
  Histos (TString prefix) :
    h_deta_ee (prefix + "_deta_ee", "deta_ee", 30, 0, 6) ,
    h_m_ee    (prefix + "_m_ee",    "m_ee",    250, 0, 1000) ,
    h_mWW     (prefix + "_mWW",     "mWW"  ,   200, 0, 400) 
    {
      cout << "Histos with prefix " << prefix << " created\n" ;
    }
  void fill (float deta_ee, float m_ee, float mWW)
    {
      h_deta_ee.Fill (deta_ee) ;
      h_m_ee.Fill (m_ee) ;
      h_mWW.Fill (mWW) ;
    }
  void save (TFile & output)
    {
      output.cd () ;
      h_deta_ee.Write () ;
      h_m_ee.Write () ;
      h_mWW.Write () ;
    }
  
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void loopAndFill (LHEF::Reader & reader, Histos & histos)
{
  int ieve = 0 ;
  //PG loop over events
  while ( reader.readEvent() ) 
    {
      ++ieve;
      if (ieve % 10000 == 0) std::cerr << "event " << ieve << "\n" ;
  
      vector<lorentzVector> eles ;
      vector<lorentzVector> Ws ;
      
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
//          std::cerr << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          if (reader.hepeup.ISTUP.at (iPart) < 0) 
            {
              continue ;
            }
          lorentzVector particle 
            (
              reader.hepeup.PUP.at (iPart).at (0), //PG px
              reader.hepeup.PUP.at (iPart).at (1), //PG py
              reader.hepeup.PUP.at (iPart).at (2), //PG pz
              reader.hepeup.PUP.at (iPart).at (3) //PG E
            ) ;
          //PG leptons
          if (abs (reader.hepeup.IDUP.at (iPart)) == 11)   //PG electron
            {
              eles.push_back (particle) ;
            }
          else if (abs (reader.hepeup.IDUP.at (iPart)) == 24)  //PG W
            {
              Ws.push_back (particle) ;
            }
        } //PG loop over particles in the event
      lorentzVector lv_ee = eles.at (0) + eles.at (1) ;
      lorentzVector lv_WW = Ws.at (0) + Ws.at (1) ;
      
      histos.fill (fabs (eles.at (0).eta () - eles.at (1).eta ()), sqrt (lv_ee.M2 ()), sqrt (lv_WW.M2 ())) ;
    } //PG loop over events


  return ;

}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char** argv) 
{
  if (argc != 3)
    {
      std::cout << ">>> Usage:   " << argv[0] << " LHE_withHiggs.lhe LHE_withoutHiggs.lhe" << std::endl;
      return -1;
    }

  ifstream ifs_wH (argv[1]) ;
  LHEF::Reader reader_wH (ifs_wH) ;
  Histos h_wH ("wH") ;
  loopAndFill (reader_wH, h_wH) ;

  ifstream ifs_wo (argv[1]) ;
  LHEF::Reader reader_wo (ifs_wo) ;
  Histos h_wo ("wo") ;
  loopAndFill (reader_wo, h_wo) ;

  TFile outfile ("epem_01.root", "recreate") ;
  h_wH.save (outfile) ;
  h_wo.save (outfile) ;
  outfile.Close () ;
  
  return 0 ;
}

