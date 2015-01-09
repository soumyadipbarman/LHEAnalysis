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
  c++ -o testPowheg02 `root-config --glibs --cflags` \
     -lm testPowheg02.cpp
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


map<string, TH1F *>
readSample (string sampleName, string radice, int maxevents = -1)
{
  cout << "reading " << sampleName << endl ;
  std::ifstream ifs (sampleName.c_str ()) ;
  LHEF::Reader reader (ifs) ;
  
  map<string, TH1F *> histos ;

  TH1F * h_higgs_eta = addHistoToMap (histos, string ("higgs_eta_")  + radice, 40, -6, 6) ;
  TH1F * h_higgs_pt  = addHistoToMap (histos, string ("higgs_pt_")   + radice, 100, 0, 500) ;
  TH1F * h_higgs_phi = addHistoToMap (histos, string ("higgs_phi_")  + radice, 30, -3.14, 3.14) ;

  TH1F * h_vbf0_eta  = addHistoToMap (histos, string ("vbf0_eta_")   + radice, 40, -6, 6) ;
  TH1F * h_vbf0_pt   = addHistoToMap (histos, string ("vbf0_pt_")    + radice, 100, 0, 400) ;
  TH1F * h_vbf0_phi  = addHistoToMap (histos, string ("vbf0_phi_")   + radice, 30, -3.14, 3.14) ;
                                                                    
  TH1F * h_vbf1_eta  = addHistoToMap (histos, string ("vbf1_eta_")   + radice, 40, -6, 6) ;
  TH1F * h_vbf1_pt   = addHistoToMap (histos,  string ("vbf1_pt_")   + radice, 100, 0, 250) ;
  TH1F * h_vbf1_phi  = addHistoToMap (histos, string ("vbf1_phi_")   + radice, 30, -3.14, 3.14) ;
                                                                    
  TH1F * h_mjj_vbf   = addHistoToMap (histos, string ("mjj_vbf_")    + radice, 70, 0, 4000) ;
  TH1F * h_NJ        = addHistoToMap (histos, string ("NJ_")         + radice, 5, 0, 5) ;
    
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
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7 ||  // quarks
                  abs (reader.hepeup.IDUP.at (iPart)) == 21 ) // gluons
                {
                  finalJets.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } // jets
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
      h_NJ->Fill (finalJets.size ()) ;

//      cout << "Njs = " << finalJets.size () << "\n" ;
//      cout << "pts: " 
//           << finalJets.at (0).Pt () 
//           << "\t"
//           << finalJets.at (1).Pt () 
//           << "\n" ;
      
      h_higgs_eta->Fill (higgs.at (0).Eta ()) ;            
      h_higgs_phi->Fill (higgs.at (0).Phi ()) ;            
      h_higgs_pt-> Fill (higgs.at (0).Pt ()) ;        

      h_vbf0_eta->Fill (finalJets.at (0).Eta ()) ;            
      h_vbf0_phi->Fill (finalJets.at (0).Phi ()) ;            
      h_vbf0_pt-> Fill (finalJets.at (0).Pt ()) ;        

      h_vbf1_eta->Fill (finalJets.at (1).Eta ()) ;            
      h_vbf1_phi->Fill (finalJets.at (1).Phi ()) ;            
      h_vbf1_pt-> Fill (finalJets.at (1).Pt ()) ;        

      float mjj = (finalJets.at (0) + finalJets.at (1)).M () ;
      h_mjj_vbf->Fill (mjj) ;

    } // loop over events
    
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  int Ntot = 100000 ;
//  std::ifstream ifs (argv[1]) ;
  map<string, TH1F *> hmap_V1 = 
    readSample ("/Users/govoni/data/powheg_test/POWHEG-BOX-V2_test_VBF_H/LHE/VBF_H_V1.lhe", "V1", 100000) ;
  map<string, TH1F *> hmap_V2 = 
//    readSample ("/Users/govoni/data/powheg_test/POWHEG-BOX-V2_test_VBF_H/LHE/VBF_H_V2.lhe", "V2", Ntot) ;
    readSample ("/Users/govoni/data/powheg_test/VBF_H_V2_1M/VBF_H_V2.lhe", "V2", 1000000) ;

  TFile outfile ("testPowheg02.root", "recreate") ;
  savemap (hmap_V1,  outfile,  1./100000.) ; 
  savemap (hmap_V2,  outfile,  1./1000000.) ; 
  
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
  
  map<string, TH1F *>::iterator iMap_V2 = hmap_V2.begin () ;

  // plotting
  for (map<string, TH1F *>::iterator iMap_V1 = hmap_V1.begin () ;
       iMap_V1 != hmap_V1.end () ;
       ++iMap_V1)
    {
      iMap_V1->second->SetStats (0) ;
      iMap_V2->second->SetStats (0) ;
      
      iMap_V1->second->SetTitle ("") ;
      iMap_V2->second->SetTitle ("") ;
      
      iMap_V1->second->SetLineWidth (2) ;
      iMap_V2->second->SetLineWidth (2) ;
      
      iMap_V1->second->SetLineColor (30) ;
      iMap_V2->second->SetLineColor (9) ;
      
      iMap_V1->second->SetFillColor (30) ;

      TLegend leg (0.3, 0.9, 0.8, 1) ;
      leg.SetNColumns (3) ;
      leg.SetLineStyle (0) ;
      leg.SetFillStyle (0) ;
      leg.AddEntry (iMap_V1->second, "V1", "fl") ;
      leg.AddEntry (iMap_V2->second, "V2", "fl") ;
      
      iMap_V1->second->GetXaxis ()->SetTitle (iMap_V1->first.c_str ()) ;        
      iMap_V2->second->GetXaxis ()->SetTitle (iMap_V2->first.c_str ()) ;        
      iMap_V2->second->Draw () ;           
      iMap_V1->second->Draw ("histsame") ;   
      iMap_V2->second->SetFillStyle (3001) ;
      iMap_V2->second->SetFillColor (9) ;
      iMap_V2->second->DrawCopy ("E2same") ;           
      iMap_V2->second->SetFillStyle (0) ;
      iMap_V2->second->SetFillColor (0) ;
      iMap_V2->second->Draw ("same") ;           
      leg.Draw () ;
      c1.Print ((iMap_V1->first + ".png").c_str (), "png") ;

      ++iMap_V2 ;
    }   
  
  return 0 ;
}


