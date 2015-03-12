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
  c++ -o lookAtMll `root-config --glibs --cflags` \
     -lm lookAtMll.cpp
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
addHistoToMap (map<string, TH1F *> & hmap, string name, string radice, int bin, float min, float max)
{
  TH1F * dummy = new TH1F ((name + "_" + radice).c_str (), (name + "_" + radice).c_str (), bin, min, max) ;
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


void fillVector (vector<TLorentzVector> & objects, LHEF::Reader & reader, int iPart)
  {
    objects.push_back (TLorentzVector
      (
        reader.hepeup.PUP.at (iPart).at (0), //PG px
        reader.hepeup.PUP.at (iPart).at (1), //PG py
        reader.hepeup.PUP.at (iPart).at (2), //PG pz
        reader.hepeup.PUP.at (iPart).at (3) //PG E
      )) ;
    return ;
  }


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


map<string, TH1F *>
readSample (string sampleName, string radice, int maxevents = -1, TH1F * weight_histo = 0)
{
  cout << "reading " << sampleName << endl ;
  std::ifstream ifs (sampleName.c_str ()) ;
  LHEF::Reader reader (ifs) ;
  
  map<string, TH1F *> histos ;

  TH1F * h_vbf0_eta    = addHistoToMap (histos, string ("vbf0_eta")     , radice, 10, 0, 6) ;
  TH1F * h_vbf0_pt     = addHistoToMap (histos, string ("vbf0_pt")      , radice, 10, 0, 400) ;
  TH1F * h_vbf0_phi    = addHistoToMap (histos, string ("vbf0_phi")     , radice, 10, -3.14, 3.14) ;
                                                                        
  TH1F * h_vbf1_eta    = addHistoToMap (histos, string ("vbf1_eta")     , radice, 10, 0, 6) ;
  TH1F * h_vbf1_pt     = addHistoToMap (histos, string ("vbf1_pt")      , radice, 10, 0, 250) ;
  TH1F * h_vbf1_phi    = addHistoToMap (histos, string ("vbf1_phi")     , radice, 10, -3.14, 3.14) ;
                                                                        
  TH1F * h_mjj_vbf     = addHistoToMap (histos, string ("mjj_vbf")      , radice, 15, 0, 4000) ;
  TH1F * h_deta_vbf    = addHistoToMap (histos, string ("deta_vbf")     , radice, 15, 0, 10) ;

  TH1F * h_mll         = addHistoToMap (histos, string ("mll")          , radice, 150, 0, 600) ;

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
      vector<TLorentzVector> finalLeptons ;      
      vector<TLorentzVector> finalNeutrinos ;      
      
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
             fillVector (initialQuarks, reader, iPart) ;
            } //PG incoming particle          

          //PG outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              // jets  
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7) // quarks
                {
                  fillVector (finalJets, reader, iPart) ;
                  fillVector (finalQuarks, reader, iPart) ;
                } // quarks
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 21 ) // gluons
                {
                  fillVector (finalJets, reader, iPart) ;
                  fillVector (finalGluon, reader, iPart) ;
                } // gluons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 25) // higgs
                {
                  fillVector (higgs, reader, iPart) ;
                } //PG higgs
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 13) // electron or muon
                {
                  fillVector (finalLeptons, reader, iPart) ;
                } // electron or muon
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 14) // electron or muon neutrinos
                {
                  fillVector (finalNeutrinos, reader, iPart) ;
                } // electron or muon neutrinos
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 9) // gluon?
                {
                  cout << "found gluon with pddgID == 9\n" ;
                  fillVector (finalJets, reader, iPart) ;
                } //PG gluon?
            } //PG outgoing particles
        } //PG loop over particles in the event
      
      float weight = 1. ;
      if (weight_histo != 0)
         weight = weight_histo->GetBinContent (weight_histo->GetXaxis ()->FindBin (higgs.at (0).Pt ())) ;

      if (finalLeptons.size () != 2) continue ;
      if (finalLeptons.at (0).Pt () < 20) continue ;
      if (finalLeptons.at (1).Pt () < 20) continue ;

      if (finalJets.size () < 2) continue ;
      sort (finalJets.rbegin (), finalJets.rend (), ptSort ()) ;      
      if (finalJets.at (0).Pt () < 30) continue ;
      if (finalJets.at (1).Pt () < 30) continue ;

      float mjj = (finalJets.at (0) + finalJets.at (1)).M () ;
      float detajj = fabs (finalJets.at (0).Eta () - finalJets.at (1).Eta ()) ;
      float detall = fabs (finalLeptons.at (0).Eta () - finalLeptons.at (1).Eta ()) ;

      if (mjj < 700) continue ;
      if (detajj < 3) continue ;
      if (detall > 2) continue ;
      
      float MET = (finalNeutrinos.at (0) + finalNeutrinos.at (1)).Pt () ;
      if (MET < 20) continue ;
      
      float mll = (finalLeptons.at (0) + finalLeptons.at (1)).M () ;

      h_mll->Fill (mll, weight) ;

      // get the tag jets
      sort (finalJets.rbegin (), finalJets.rend (), ptSort ()) ;
      sort (finalQuarks.rbegin (), finalQuarks.rend (), ptSort ()) ;

//      cout << "Njs = " << finalJets.size () << "\n" ;
//      cout << "pts: " 
//           << finalJets.at (0).Pt () 
//           << "\t"
//           << finalJets.at (1).Pt () 
//           << "\n" ;
        
      h_vbf0_eta->Fill (fabs (finalJets.at (0).Eta ()), weight) ;            
      h_vbf0_phi->Fill (finalJets.at (0).Phi (), weight) ;            
      h_vbf0_pt-> Fill (finalJets.at (0).Pt (), weight) ;        

      h_vbf1_eta->Fill (fabs (finalJets.at (1).Eta ()), weight) ;            
      h_vbf1_phi->Fill (finalJets.at (1).Phi (), weight) ;            
      h_vbf1_pt-> Fill (finalJets.at (1).Pt (), weight) ;        

      h_mjj_vbf->Fill (mjj, weight) ;
      h_deta_vbf->Fill (detajj, weight) ;

    } // loop over events
    
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  int NTOT = 100000 ;
  if (argc > 1) NTOT = atoi (argv[1]) ;

  int N_SMH_tot = NTOT ;
  map<string, TH1F *> hmap_SMH = 
    readSample ("/Users/govoni/data/TP/phantom/gen_TP_uvev_126/total.lhe", "SMH", N_SMH_tot) ;

  int N_noH_tot = NTOT ;
  map<string, TH1F *> hmap_noH = 
    readSample ("/Users/govoni/data/TP/phantom/gen_TP_uvev_noH/total.lhe", "noH", N_noH_tot) ;

  string outFolderName = "lookAtMll_plots/";
  system (Form ("mkdir -p %s", outFolderName.c_str ())) ;

  TFile outfile ((outFolderName + "/lookAtMll.root").c_str (), "recreate") ;
  savemap (hmap_SMH,  outfile,  4.1365 / N_SMH_tot) ; // XS in fb
  savemap (hmap_noH,  outfile,  4.4920 / N_noH_tot) ; // XS in fb
  
  // plotting
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

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

  for (map<string, TH1F *>::iterator iMap_noH = hmap_noH.begin () ;
       iMap_noH != hmap_noH.end () ;
       ++iMap_noH)
    {
      hmap_noH[iMap_noH->first]->SetLineWidth (2) ;
      hmap_SMH[iMap_noH->first]->SetLineWidth (2) ;

      hmap_SMH[iMap_noH->first]->SetStats (0) ;
      hmap_noH[iMap_noH->first]->SetStats (0) ;

      hmap_noH[iMap_noH->first]->SetLineColor (2) ;

      float min = 100. ;     
      if (hmap_noH[iMap_noH->first]->GetMinimum () < min) min = hmap_noH[iMap_noH->first]->GetMinimum () ;
      if (hmap_SMH[iMap_noH->first]->GetMinimum () < min) min = hmap_SMH[iMap_noH->first]->GetMinimum () ;
      if (min < 0) min *= 1.1 ;
      if (min > 0) min *= 0.9 ;

      float max = 0. ;     
      if (hmap_noH[iMap_noH->first]->GetMaximum () > max) max = hmap_noH[iMap_noH->first]->GetMaximum () ;
      if (hmap_SMH[iMap_noH->first]->GetMaximum () > max) max = hmap_SMH[iMap_noH->first]->GetMaximum () ;
      if (max < 0) max *= 0.9 ;
      if (max > 0) max *= 1.1 ;

      c1.DrawFrame (hmap_noH[iMap_noH->first]->GetXaxis ()->GetXmin (), min, 
                    hmap_noH[iMap_noH->first]->GetXaxis ()->GetXmax (), max) ;

      hmap_noH[iMap_noH->first]->Draw ("Ehisto same") ;
      hmap_SMH[iMap_noH->first]->Draw ("Ehisto same") ;

      TLegend leg (0.3, 0.9, 0.8, 1) ;
      leg.SetNColumns (3) ;
      leg.SetLineStyle (0) ;
      leg.SetFillStyle (0) ;
      leg.AddEntry (hmap_noH[iMap_noH->first], "noH", "fl") ;
      leg.AddEntry (hmap_SMH[iMap_noH->first], "SM", "fl") ;
      leg.Draw () ;

      c1.Print ((outFolderName + iMap_noH->first + ".png").c_str (), "png") ;
    }  

  return 0 ;
}


