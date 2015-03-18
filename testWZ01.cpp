/*
  c++ -o testWZ01 `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags` \
     -lm testWZ01.cpp
*/



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

//#include "Math/Vector3D.h"
//#include "Math/Vector4D.h"
//using namespace ROOT::Math ;

#include "LHAPDF/LHAPDF.h"
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
readSample (string sampleName, string radice, float referenceScale = 0., int maxevents = -1)
{
  cout << "reading " << sampleName << endl ;
  std::ifstream ifs (sampleName.c_str ()) ;
  LHEF::Reader reader (ifs) ;
  
  map<string, TH1F *> histos ;

  TH1F * h_vbf0_eta    = addHistoToMap (histos, string ("vbf0_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_vbf0_pt     = addHistoToMap (histos, string ("vbf0_pt_")      + radice, 100, 0, 400) ;
  TH1F * h_vbf0_phi    = addHistoToMap (histos, string ("vbf0_phi_")     + radice, 30, -3.14, 3.14) ;
                                                                        
  TH1F * h_vbf1_eta    = addHistoToMap (histos, string ("vbf1_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_vbf1_pt     = addHistoToMap (histos, string ("vbf1_pt_")      + radice, 100, 0, 250) ;
  TH1F * h_vbf1_phi    = addHistoToMap (histos, string ("vbf1_phi_")     + radice, 30, -3.14, 3.14) ;
                                                                        
  TH1F * h_lep0_eta    = addHistoToMap (histos, string ("lep0_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_lep0_pt     = addHistoToMap (histos, string ("lep0_pt_")      + radice, 100, 0, 400) ;
  TH1F * h_lep0_phi    = addHistoToMap (histos, string ("lep0_phi_")     + radice, 30, -3.14, 3.14) ;
                                                                        
  TH1F * h_lep1_eta    = addHistoToMap (histos, string ("lep1_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_lep1_pt     = addHistoToMap (histos, string ("lep1_pt_")      + radice, 100, 0, 250) ;
  TH1F * h_lep1_phi    = addHistoToMap (histos, string ("lep1_phi_")     + radice, 30, -3.14, 3.14) ;
                                                                        
  TH1F * h_lep2_eta    = addHistoToMap (histos, string ("lep2_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_lep2_pt     = addHistoToMap (histos, string ("lep2_pt_")      + radice, 100, 0, 250) ;
  TH1F * h_lep2_phi    = addHistoToMap (histos, string ("lep2_phi_")     + radice, 30, -3.14, 3.14) ;
                                                                        
  TH1F * h_met_eta    = addHistoToMap (histos, string ("met_eta_")     + radice, 40, -6, 6) ;
  TH1F * h_met_pt     = addHistoToMap (histos, string ("met_pt_")      + radice, 100, 0, 250) ;
  TH1F * h_met_phi    = addHistoToMap (histos, string ("met_phi_")     + radice, 30, -3.14, 3.14) ;
                                                                        
  TH1F * h_mjj_vbf     = addHistoToMap (histos, string ("mjj_vbf_")      + radice, 70, 0, 4000) ;
  TH1F * h_deta_vbf    = addHistoToMap (histos, string ("deta_vbf_")     + radice, 70, 0, 10) ;
                                                                    
  TH1F * h_NJ          = addHistoToMap (histos, string ("NJ_")           + radice, 5, 0, 5) ;
  TH1F * h_NG          = addHistoToMap (histos, string ("NG_")           + radice, 5, 0, 5) ;

  TH1F * h_scale       = addHistoToMap (histos, string ("scale_")        + radice, 100, 0., 500.) ;
  TH1F * h_weight      = addHistoToMap (histos, string ("weight_")       + radice, 100, 0., 10.) ;
    
  int ieve = 0 ;
  // loop over events
  while ( reader.readEvent () ) 
    {
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
      if (maxevents > 0 && ieve >= maxevents) break ;
      ++ieve;
  
      vector<TLorentzVector> finalJets ;      
      vector<TLorentzVector> initialQuarks ;      
      vector<TLorentzVector> finalQuarks ;      
      vector<TLorentzVector> finalGluons ;      
      vector<TLorentzVector> leptons ;      
      vector<TLorentzVector> neutrinos ;      
      
      double x[2] = {0., 0.} ;
      int flavour[2] = {0, 0} ;

      int iquark = 0 ;
      //PG loop over particles in the event
      //PG and fill the variables of leptons and quarks
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
//          std::cout << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          TLorentzVector particle
                (
                  reader.hepeup.PUP.at (iPart).at (0), //PG px
                  reader.hepeup.PUP.at (iPart).at (1), //PG py
                  reader.hepeup.PUP.at (iPart).at (2), //PG pz
                  reader.hepeup.PUP.at (iPart).at (3) //PG E
                ) ;
          //PG incoming particle          
          if (reader.hepeup.ISTUP.at (iPart) == -1)
            {
               x[iquark] = particle.P () / 7000. ;
               flavour[iquark++] = reader.hepeup.IDUP.at (iPart) ;
               initialQuarks.push_back (particle) ;
            } //PG incoming particle          

          //PG outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              // jets  
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7) // quarks
                {
                  finalJets.push_back (particle) ;
                  finalQuarks.push_back (particle) ;
                } // quarks
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 21 ) // gluons
                {
                  finalJets.push_back (particle) ;
                  finalGluons.push_back (particle) ;
                } // gluons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 13) // charged leptons (e,u)
                {
                  leptons.push_back (particle) ;
                } //PG charged leptons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 14) // neutrinos
                {
                  neutrinos.push_back (particle) ;
                } //PG neutrinos
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 15 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 16) // charged leptons (tau)
                {
                  cout << "WARNING third family present!" << endl ; 
                } //PG charged leptons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 9) // gluon?
                {
                  cout << "WARNING found gluon with pddgID == 9\n" ;
                  finalJets.push_back (particle) ;
                } //PG gluon?
            } //PG outgoing particles
        } //PG loop over particles in the event

      double weight = 1. ;
      float scale = reader.hepeup.SCALUP ;
      if (referenceScale != 0 )
        weight = LHAPDF::xfx (x[0], referenceScale, flavour[0]) * LHAPDF::xfx (x[1], referenceScale, flavour[1]) /
                 (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;

      h_weight->Fill (weight) ;
      h_scale->Fill (scale) ;

      if (isnan (weight))
        {
          cout << "WARNING weight is not a number, setting to 1." << endl ;
          cout << "\t x0: " << x[0] << endl ; 
          cout << "\t x1: " << x[1] << endl ; 
          cout << "\t num, p0: " << LHAPDF::xfx (x[0], referenceScale, flavour[0]) << endl ; 
          cout << "\t num, p1: " << LHAPDF::xfx (x[1], referenceScale, flavour[1]) << endl ; 
          cout << "\t den, p0: " << LHAPDF::xfx (x[0], scale, flavour[0]) << endl ; 
          cout << "\t den, p1: " << LHAPDF::xfx (x[1], scale, flavour[1]) << endl ; 
          weight = 1. ;
        } 
      
      if (isinf (weight))
        {
          cout << "WARNING weight is infinite, setting to 1." << endl ;
          cout << "\t num: " << (LHAPDF::xfx (x[0], referenceScale, flavour[0]) * LHAPDF::xfx (x[1], referenceScale, flavour[1])) << endl ;
          cout << "\t den: " << (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) << endl ;
          weight = 1. ;
        } 
      
      sort (leptons.rbegin (), leptons.rend (), ptSort ()) ;
      sort (finalJets.rbegin (), finalJets.rend (), ptSort ()) ;

      float mjj = (finalJets.at (0) + finalJets.at (1)).M () ;
      float detajj = fabs (finalJets.at (0).Eta () - finalJets.at (1).Eta ()) ;

      //CUTS
      if (leptons.at (0).Pt () < 20) continue ;
      if (leptons.at (1).Pt () < 20) continue ;
      if (leptons.at (2).Pt () < 20) continue ;

      if (finalJets.at (0).Pt () < 20) continue ;
      if (finalJets.at (1).Pt () < 20) continue ;
      
      if (mjj < 300) continue ; 

      // get the tag jets
      h_NJ->Fill (finalJets.size (), weight) ;
      h_NG->Fill (finalGluons.size (), weight) ;

      h_vbf0_eta->Fill (finalJets.at (0).Eta (), weight) ;            
      h_vbf0_phi->Fill (finalJets.at (0).Phi (), weight) ;            
      h_vbf0_pt-> Fill (finalJets.at (0).Pt (), weight) ;        

      h_vbf1_eta->Fill (finalJets.at (1).Eta (), weight) ;            
      h_vbf1_phi->Fill (finalJets.at (1).Phi (), weight) ;            
      h_vbf1_pt-> Fill (finalJets.at (1).Pt (), weight) ;        

      h_mjj_vbf->Fill (mjj, weight) ;
      h_deta_vbf->Fill (detajj, weight) ;

      h_lep0_eta->Fill (leptons.at (0).Eta (), weight) ;            
      h_lep0_phi->Fill (leptons.at (0).Phi (), weight) ;            
      h_lep0_pt-> Fill (leptons.at (0).Pt (), weight) ;        

      h_lep1_eta->Fill (leptons.at (1).Eta (), weight) ;            
      h_lep1_phi->Fill (leptons.at (1).Phi (), weight) ;            
      h_lep1_pt-> Fill (leptons.at (1).Pt (), weight) ;        

      h_lep2_eta->Fill (leptons.at (2).Eta (), weight) ;            
      h_lep2_phi->Fill (leptons.at (2).Phi (), weight) ;            
      h_lep2_pt-> Fill (leptons.at (2).Pt (), weight) ;        

      h_met_eta->Fill (neutrinos.at (0).Eta (), weight) ;            
      h_met_phi->Fill (neutrinos.at (0).Phi (), weight) ;            
      h_met_pt-> Fill (neutrinos.at (0).Pt (), weight) ;        

    } // loop over events
    
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  const int SUBSET = 0 ;
  const string NAME = "cteq6ll" ; //"cteq6l1"

  LHAPDF::initPDFSet (NAME, LHAPDF::LHPDF, SUBSET) ;
  const int NUMBER = LHAPDF::numberPDF () ;

  LHAPDF::initPDF (0) ;

  gROOT->SetStyle ("Plain") ;

//  float commonScale = 125. ;
  float commonScale = 140. ;

//  int N_PH_tot = 480000 ;
  int N_PH_tot = 80000 ;
  float PH_xs = 32.58 ; // fb
  map<string, TH1F *> hmap_PH = 
  readSample ("/Users/govoni/data/TP/WZ/Phantom_QCD/total.lhe", "PH", commonScale, N_PH_tot) ;

//  int N_MG_tot = 400000 ;
  int N_MG_tot = 80000 ;
  float MG_xs = 287.59 ; // fb
  map<string, TH1F *> hmap_MG = 
  readSample ("/Users/govoni/data/TP/WZ/Madgraph_QCD/total.lhe", "MG", commonScale, N_MG_tot) ;

  // compare shapes only
//  PH_xs = 1. ;
//  MG_xs = 1. ; 
    
  TFile outfile ("testWZ01.root", "recreate") ;
  savemap (hmap_PH,  outfile,  PH_xs / N_PH_tot) ; 
  savemap (hmap_MG,  outfile,  MG_xs / N_MG_tot) ; 
  
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
  

  string outFolderName = "testWZ01_plots/";
  system (Form ("mkdir -p %s", outFolderName.c_str ())) ;

  // plotting
  map<string, TH1F *>::iterator iMap_MG = hmap_MG.begin () ;
  for (map<string, TH1F *>::iterator iMap_PH = hmap_PH.begin () ;
       iMap_PH != hmap_PH.end () ;
       ++iMap_PH)
    {
      iMap_PH->second->SetStats (0) ;
      iMap_MG->second->SetStats (0) ;
      
      iMap_PH->second->SetTitle ("") ;
      iMap_MG->second->SetTitle ("") ;
      
      iMap_PH->second->SetLineWidth (2) ;
      iMap_MG->second->SetLineWidth (2) ;
      
      iMap_PH->second->SetLineColor (30) ;
      iMap_MG->second->SetLineColor (9) ;
      
      iMap_PH->second->SetFillColor (30) ;

      TLegend leg (0.3, 0.9, 0.8, 1) ;
      leg.SetNColumns (3) ;
      leg.SetLineStyle (0) ;
      leg.SetFillStyle (0) ;
      leg.AddEntry (iMap_PH->second, "PH", "fl") ;
      leg.AddEntry (iMap_MG->second, "MG", "fl") ;
      
      iMap_PH->second->GetXaxis ()->SetTitle (iMap_PH->first.c_str ()) ;        
      iMap_MG->second->GetXaxis ()->SetTitle (iMap_MG->first.c_str ()) ;        
      iMap_MG->second->Draw () ;           
      iMap_PH->second->Draw ("histsame") ;   
      iMap_MG->second->SetFillStyle (3001) ;
      iMap_MG->second->SetFillColor (9) ;
      iMap_MG->second->DrawCopy ("E2same") ;           
      iMap_MG->second->SetFillStyle (0) ;
      iMap_MG->second->SetFillColor (0) ;
      iMap_MG->second->Draw ("same") ;           
      leg.Draw () ;
      
      c1.Print ((outFolderName + iMap_PH->first + ".png").c_str (), "png") ;

      ++iMap_MG ;
    }   

  return 0 ;
}


