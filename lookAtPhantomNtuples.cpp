#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <cmath>

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

/*
  c++ -o lookAtPhantomNtuples `root-config --glibs --cflags` \
     -lm lookAtPhantomNtuples.cpp
*/


using namespace std ;


double 
deltaPhi (double phi1, double phi2)
{
  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  return deltaphi;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -


struct ptSort: 
public std::binary_function<TLorentzVector &, TLorentzVector &, bool>
{
  bool operator() (TLorentzVector & x, TLorentzVector & y)
    {
      return x.Pt () < y.Pt () ;
    }
} ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -


struct etaSort: 
public std::binary_function<TLorentzVector &, TLorentzVector &, bool>
{
  bool operator() (TLorentzVector & x, TLorentzVector & y)
    {
      return x.Eta () < y.Eta () ;
    }
} ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  std::ifstream ifs ("/Users/govoni/data/TP/phantom/testNtuples.lhe") ;
  LHEF::Reader reader (ifs) ;

  TH1F h_ptj1   ("ptj1  ", "ptj1  ", 100, 0., 500) ;
  TH1F h_etaj1  ("etaj1 ", "etaj1 ", 50, -5., 5.) ;
  TH1F h_phij1  ("phij1 ", "phij1 ", 50, -3.14, 3.14) ;
  TH1F h_ptj2   ("ptj2  ", "ptj2  ", 100, 0., 500) ;
  TH1F h_etaj2  ("etaj2 ", "etaj2 ", 50, -5., 5.) ;
  TH1F h_phij2  ("phij2 ", "phij2 ", 50, -3.14, 3.14) ;
  TH1F h_ptl1   ("ptl1  ", "ptl1  ", 100, 0., 500) ;
  TH1F h_etal1  ("etal1 ", "etal1 ", 50, -5., 5.) ;
  TH1F h_phil1  ("phil1 ", "phil1 ", 50, -3.14, 3.14) ;
  TH1F h_ptl2   ("ptl2  ", "ptl2  ", 100, 0., 500) ;
  TH1F h_etal2  ("etal2 ", "etal2 ", 50, -5., 5.) ;
  TH1F h_phil2  ("phil2 ", "phil2 ", 50, -3.14, 3.14) ;
  TH1F h_ptv1   ("ptv1  ", "ptv1  ", 100, 0., 500) ;
  TH1F h_etav1  ("etav1 ", "etav1 ", 50, -5., 5.) ;
  TH1F h_phiv1  ("phiv1 ", "phiv1 ", 50, -3.14, 3.14) ;
  TH1F h_ptv2   ("ptv2  ", "ptv2  ", 100, 0., 500) ;
  TH1F h_etav2  ("etav2 ", "etav2 ", 50, -5., 5.) ;
  TH1F h_phiv2  ("phiv2 ", "phiv2 ", 50, -3.14, 3.14) ;
  TH1F h_mjj    ("mjj   ", "mjj   ", 100, 0., 2000) ;
  TH1F h_detajj ("detajj", "detajj", 50 , 0., 10.) ;
  TH1F h_mll    ("mll   ", "mll   ", 100, 0., 2000) ;

  int ieve = 0 ;
  int maxevents = -1 ;
  // loop over events
  while ( reader.readEvent () ) 
    {
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
      if (maxevents > 0 && ieve >= maxevents) break ;
      ++ieve;
  
      vector<TLorentzVector> leptons ;      
      vector<TLorentzVector> neutrinos ;      
      vector<TLorentzVector> finalQuarks ;      
      vector<TLorentzVector> initialQuarks ;      
      
      //PG loop over particles in the event
      //PG and fill the variables of leptons and quarks
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
//          std::cout << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          //PG incoming particle          
          if (reader.hepeup.ISTUP.at (iPart) == -1)
            {
              initialQuarks.push_back (TLorentzVector
                (
                  reader.hepeup.PUP.at (iPart).at (0), //PG px
                  reader.hepeup.PUP.at (iPart).at (1), //PG py
                  reader.hepeup.PUP.at (iPart).at (2), //PG pz
                  reader.hepeup.PUP.at (iPart).at (3) //PG E
                )) ;
            } //PG incoming particle          

          //PG outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              //PG leptons
              if (abs (reader.hepeup.IDUP.at (iPart)) == 11 || //PG electron
                  abs (reader.hepeup.IDUP.at (iPart)) == 13 || //PG muon
                  abs (reader.hepeup.IDUP.at (iPart)) == 15)   //PG tau
                {
                  leptons.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3) //PG E
                    )) ;
                } //PG leptons
           else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 14 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 16)
                {
                  neutrinos.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } //PG neutrinos
           else if (abs (reader.hepeup.IDUP.at (iPart)) < 7) 
                {
                  finalQuarks.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3) //PG E
                    )) ;
                }
            } //PG outgoing particles
        } //PG loop over particles in the event

      if (leptons.size () != 2) {cout << "not 2 leptons\n" ; continue ; }
      if (neutrinos.size () != 2) {cout << "not 2 neutrinos\n" ; continue ; }
      if (finalQuarks.size () != 2) {cout << "not 2 finalQuarks\n" ; continue ; }

      sort (leptons.rbegin (), leptons.rend (), ptSort ()) ;
      sort (neutrinos.rbegin (), neutrinos.rend (), ptSort ()) ;
      sort (finalQuarks.rbegin (), finalQuarks.rend (), ptSort ()) ;

      h_ptj1  .Fill (finalQuarks.at (0).Pt ()) ;
      h_etaj1 .Fill (finalQuarks.at (0).Eta ()) ;
      h_phij1 .Fill (finalQuarks.at (0).Phi ()) ;
      h_ptj2  .Fill (finalQuarks.at (1).Pt ()) ;
      h_etaj2 .Fill (finalQuarks.at (1).Eta ()) ;
      h_phij2 .Fill (finalQuarks.at (1).Phi ()) ;
      h_ptl1  .Fill (leptons.at (0).Pt ()) ;
      h_etal1 .Fill (leptons.at (0).Eta ()) ;
      h_phil1 .Fill (leptons.at (0).Phi ()) ;
      h_ptl2  .Fill (leptons.at (1).Pt ()) ;
      h_etal2 .Fill (leptons.at (1).Eta ()) ;
      h_phil2 .Fill (leptons.at (1).Phi ()) ;
      h_ptv1  .Fill (neutrinos.at (0).Pt ()) ;
      h_etav1 .Fill (neutrinos.at (0).Eta ()) ;
      h_phiv1 .Fill (neutrinos.at (0).Phi ()) ;
      h_ptv2  .Fill (neutrinos.at (1).Pt ()) ;
      h_etav2 .Fill (neutrinos.at (1).Eta ()) ;
      h_phiv2 .Fill (neutrinos.at (1).Phi ()) ;

      float mjj = (finalQuarks.at (0) + finalQuarks.at (1)).M () ;
      h_mjj   .Fill (mjj) ;
      h_detajj.Fill (fabs (finalQuarks.at (1).Eta () - finalQuarks.at (0).Eta ())) ;
      float mll = (leptons.at (0) + leptons.at (1)).M () ;
      h_mll   .Fill (mll) ;
    } // loop over events

  cout << "RUN OVER " << ieve << " EVENTS\n" ;

  TFile outfile ("lookAtPhantomNtuples.root", "recreate") ;
  h_ptj1.Write () ;
  h_etaj1.Write () ;
  h_phij1.Write () ;
  h_ptj2.Write () ;
  h_etaj2.Write () ;
  h_phij2.Write () ;
  h_ptl1.Write () ;
  h_etal1.Write () ;
  h_phil1.Write () ;
  h_ptl2.Write () ;
  h_etal2.Write () ;
  h_phil2.Write () ;
  h_ptv1.Write () ;
  h_etav1.Write () ;
  h_phiv1.Write () ;
  h_ptv2.Write () ;
  h_etav2.Write () ;
  h_phiv2.Write () ;
  h_mjj.Write () ;
  h_detajj.Write () ;
  h_mll.Write () ;
  outfile.Close () ;
  
  return 0 ;
}

