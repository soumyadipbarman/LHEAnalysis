#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include <cmath>
#include "hFactory.h"
#include "h2Factory.h"
#include "hFunctions.h"

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

/*
  c++ -o testPowheg04 `root-config --glibs --cflags` \
     -lm testPowheg04.cpp
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


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


struct ptSort: public std::binary_function<TLorentzVector &, TLorentzVector &, bool>
{
  bool operator() (TLorentzVector & x, TLorentzVector & y)
    {
      return x.Pt () < y.Pt () ;
    }
} ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


struct etaSort: public std::binary_function<TLorentzVector &, TLorentzVector &, bool>
{
  bool operator() (TLorentzVector & x, TLorentzVector & y)
    {
      return x.Eta () < y.Eta () ;
    }
} ;




// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


template <class TLV>
pair<int, int> findPairWithWMass (const vector<TLV> & v_f_quarks)
{
  double ref_deltaM = 100000. ;
  int one = 0 ;
  int two = 0 ;
  
  for (int iJ = 0 ; iJ < 4 ; ++iJ)
    for (int iJ2 = iJ + 1 ; iJ2 < 4 ; ++iJ2)
      {
        double deltaM = fabs ((v_f_quarks.at (iJ) + v_f_quarks.at (iJ2)).M () - 80.4) ;
        if (deltaM < ref_deltaM)
          {
            ref_deltaM = deltaM ;
            one = iJ ;
            two = iJ2 ;
          }      
      }
  return pair<int, int> (one, two) ;
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TH1F *
addHistoToMap (map<string, TH1F *> & hmap, string name, int bin, float min, float max)
{
  TH1F * dummy = new TH1F (name.c_str (), name.c_str (), bin, min, max) ;
  dummy->Sumw2 () ;
  hmap[name] = dummy ;
  return dummy ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


void 
savemap (map<string, TH1F *> & hmap, TFile & outfile, float scale)
{
  outfile.cd () ;
  for (map<string, TH1F *>::iterator iMap = hmap.begin () ;
      iMap != hmap.end () ; ++ iMap)
    {
      iMap->second->Scale (scale) ;
      iMap->second->Write () ;
    }
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


TH1F * swap (TH1F * input)
{
  TString name = input->GetName () ;
  name += "_swap" ;
  TH1F * dummy = (TH1F *) input->Clone (name) ;
  dummy->Reset () ;
  for (int i = 0 ; i <= input->GetNbinsX () + 1 ; ++i)
    {
      dummy->SetBinContent (i, input->GetBinContent (input->GetNbinsX () + 1 - i)) ;
      dummy->SetBinError (i, input->GetBinError (input->GetNbinsX () + 1 - i)) ;
    }
  return dummy ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


map<string, TH1F *>
readSample (string sampleName, string radice, int maxevents = -1)
{
  cout << "reading " << sampleName << endl ;
  std::ifstream ifs (sampleName.c_str ()) ;
  LHEF::Reader reader (ifs) ;
  
  map<string, TH1F *> histos ;

  TH1F * h_higgs_eta   = addHistoToMap (histos, string ("higgs_eta_")    + radice, 40, -6, 6) ;
  TH1F * h_gluon_eta   = addHistoToMap (histos, string ("gluon_eta_")    + radice, 40, -6, 6) ;
  TH1F * h_vbf0_eta    = addHistoToMap (histos, string ("vbf0_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_vbf1_eta    = addHistoToMap (histos, string ("vbf1_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_quark0_eta  = addHistoToMap (histos, string ("quark0_eta_")   + radice, 40, -6, 6) ;
  TH1F * h_quark1_eta  = addHistoToMap (histos, string ("quark1_eta_")   + radice, 40, -6, 6) ;
    
  int ieve = 0 ;
  // loop over events
  while ( reader.readEvent () ) 
    {
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
      if (maxevents > 0 && ieve >= maxevents) break ;
      ++ieve;
  
      vector<TLorentzVector> higgs ;      
      vector<TLorentzVector> finalJets ;      
      vector<TLorentzVector> initialQuarks ;      
      vector<TLorentzVector> finalQuarks ;      
      vector<TLorentzVector> finalGluon ;      
      
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
              // jets  
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7) // quarks
                {
                  finalJets.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                  finalQuarks.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } // quarks
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 21 ) // gluons
                {
                  finalJets.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                  finalGluon.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } // gluons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 25) // higgs
                {
                  higgs.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } //PG higgs
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 9) // gluon?
                {
                  cout << "found gluon with pddgID == 9\n" ;
                  finalJets.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3) //PG E
                    )) ;
                } //PG gluon?
            } //PG outgoing particles
        } //PG loop over particles in the event

      // get the tag jets
      sort (finalJets.rbegin (), finalJets.rend (), ptSort ()) ;
      sort (finalQuarks.rbegin (), finalQuarks.rend (), ptSort ()) ;

//      cout << "Njs = " << finalJets.size () << "\n" ;
//      cout << "pts: " 
//           << finalJets.at (0).Pt () 
//           << "\t"
//           << finalJets.at (1).Pt () 
//           << "\n" ;

      if (finalGluon.size () != 0)
        {
          h_gluon_eta->Fill (finalGluon.at (0).Eta ()) ;            
        }
        
      h_vbf0_eta->Fill (finalJets.at (0).Eta ()) ;            
      h_vbf1_eta->Fill (finalJets.at (1).Eta ()) ;            
      h_higgs_eta->Fill (higgs.at (0).Eta ()) ;            
      h_quark0_eta->Fill (finalQuarks.at (0).Eta ()) ;            
      h_quark1_eta->Fill (finalQuarks.at (1).Eta ()) ;            

    } // loop over events
    
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

//  int N_V2_tot = 1000000 ;
  int N_V2_tot = 154000 ;
  map<string, TH1F *> hmap_V2 = 
//    readSample ("/Users/govoni/data/powheg_test/POWHEG-BOX-V2_test_VBF_H/LHE/VBF_H_V2_1M.lhe", "V2", N_V2_tot) ;
    readSample ("/Users/govoni/data/powheg_test/newRound/V2_fC.lhe", "V2", N_V2_tot) ;

  string outFolderName = "testPowheg04_plots/";
  system (Form ("mkdir -p %s", outFolderName.c_str ())) ;

  TFile outfile ((outFolderName + "testPowheg04.root").c_str (), "recreate") ;
  savemap (hmap_V2,  outfile,  1./N_V2_tot) ; 
  
  gStyle->SetStatStyle (0) ; 
  gStyle->SetTitleStyle (0) ; 
  gStyle->SetCanvasBorderSize (0) ; 
  gStyle->SetFrameBorderSize (0) ; 
  gStyle->SetLegendBorderSize (0) ; 
  gStyle->SetStatBorderSize (0) ; 
  gStyle->SetTitleBorderSize (0) ; 
  gStyle->SetTitleYOffset (2) ;

  TCanvas c1 ;
  c1.SetLeftMargin (0.17) ; 
  c1.SetTopMargin (0.1) ; 
  
  // plotting
  for (map<string, TH1F *>::iterator iMap_V2 = hmap_V2.begin () ;
       iMap_V2 != hmap_V2.end () ;
       ++iMap_V2)
    {
      iMap_V2->second->SetStats (0) ;
      
      iMap_V2->second->SetTitle ("") ;
      iMap_V2->second->SetLineWidth (2) ;
      iMap_V2->second->SetLineColor (9) ;
      
      TH1F * swapped = swap (iMap_V2->second) ;
      swapped->SetLineColor (2) ;
      swapped->SetFillColor (2) ;
//      swapped->SetFillStyle (3005) ;

      TLegend leg (0.3, 0.9, 0.8, 1) ;
      leg.SetNColumns (3) ;
      leg.SetLineStyle (0) ;
      leg.SetFillStyle (0) ;
      leg.AddEntry (iMap_V2->second, "V2", "fl") ;
      leg.AddEntry (swapped, "swapped", "fl") ;
      
      swapped->GetXaxis ()->SetTitle (iMap_V2->first.c_str ()) ;        
      iMap_V2->second->GetXaxis ()->SetTitle (iMap_V2->first.c_str ()) ;        
      swapped->Draw ("E2") ;   
      iMap_V2->second->SetFillStyle (3004) ;
      iMap_V2->second->SetFillColor (9) ;
      iMap_V2->second->DrawCopy ("E2same") ;           
      iMap_V2->second->SetFillStyle (0) ;
      iMap_V2->second->SetFillColor (0) ;
      iMap_V2->second->Draw ("E2same") ;           
      leg.Draw () ;
      
      c1.Print ((outFolderName + iMap_V2->first + ".png").c_str (), "png") ;
    }   

  return 0 ;
}


