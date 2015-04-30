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
//#include "hFactory.h"
//#include "h2Factory.h"
//#include "hFunctions.h"

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

/*
  c++ -o testWWInterf `root-config --glibs --cflags` \
     -lm testWWInterf.cpp
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

  TH1F * h_initPart  = addHistoToMap (histos, string ("initPart_")   + radice, 45, -22.5, 22.5) ;

  TH1F * h_vbf0_eta  = addHistoToMap (histos, string ("vbf0_eta_")   + radice, 40, -6, 6) ;
  TH1F * h_vbf0_pt   = addHistoToMap (histos, string ("vbf0_pt_")    + radice, 100, 0, 500) ;
  TH1F * h_vbf0_phi  = addHistoToMap (histos, string ("vbf0_phi_")   + radice, 30, -3.14, 3.14) ;
                                                                    
  TH1F * h_vbf1_eta  = addHistoToMap (histos, string ("vbf1_eta_")   + radice, 40, -6, 6) ;
  TH1F * h_vbf1_pt   = addHistoToMap (histos,  string ("vbf1_pt_")   + radice, 100, 0, 500) ;
  TH1F * h_vbf1_phi  = addHistoToMap (histos, string ("vbf1_phi_")   + radice, 30, -3.14, 3.14) ;
                                                                    
  TH1F * h_mjj_vbf   = addHistoToMap (histos, string ("mjj_vbf_")    + radice, 25, 0, 2500) ;
  TH1F * h_NJ        = addHistoToMap (histos, string ("NJ_")         + radice, 5, 0, 5) ;

  TH1F * h_m4l       = addHistoToMap (histos, string ("m4l_")        + radice, 35, 0, 2275) ;
    
  int ieve = 0 ;
  int btagged = 0 ;
  int selected = 0 ;
  // loop over events
  while ( reader.readEvent () ) 
    {
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
      if (maxevents > 0 && ieve >= maxevents) break ;
      ++ieve;
  
      vector<int> initialParticlesPDGId ;      
      vector<TLorentzVector> finalLeptons ;      
      vector<TLorentzVector> finalNeutrinos ;      
      vector<TLorentzVector> initialQuarks ;      
      vector<TLorentzVector> finalJets ;      
      vector<TLorentzVector> finalBs ;      
      
      int leptonsIDproduct = 1 ;

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
              initialParticlesPDGId.push_back (reader.hepeup.IDUP.at (iPart)) ;
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


                  if (abs (reader.hepeup.IDUP.at (iPart)) == 5)  // b's
                    {
                      finalBs.push_back (TLorentzVector
                        (
                          reader.hepeup.PUP.at (iPart).at (0), //PG px
                          reader.hepeup.PUP.at (iPart).at (1), //PG py
                          reader.hepeup.PUP.at (iPart).at (2), //PG pz
                          reader.hepeup.PUP.at (iPart).at (3)  //PG E
                        )) ;
                    } // b's
                } // jets
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 13 ) // charged leptons
                {
                  leptonsIDproduct *= reader.hepeup.IDUP.at (iPart) ;
                  finalLeptons.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } //PG charged leptons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 14 ) // neutrinos
                {
                  finalNeutrinos.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3)  //PG E
                    )) ;
                } // neutrinos

            } //PG outgoing particles
        } //PG loop over particles in the event

      if (initialParticlesPDGId.size () != 2) exit (1) ; 

      if (finalJets.size () != 2) continue ;
      if (finalNeutrinos.size () != 2) continue ;
      if (finalLeptons.size () != 2) continue ;

//      if (fabs (initialParticlesPDGId.at (0)) > 4) continue ;
//      if (fabs (initialParticlesPDGId.at (1)) > 4) continue ;

      // the following remove events with guons as initial particles
//      if (fabs (initialParticlesPDGId.at (0)) > 6) continue ;
//      if (fabs (initialParticlesPDGId.at (1)) > 6) continue ;
      
      sort (finalJets.rbegin (), finalJets.rend (), ptSort ()) ;
      if (finalJets.at (0).Pt () < 20) continue ;
      if (finalJets.at (1).Pt () < 20) continue ;
      if (finalJets.at (0).E () < 20) continue ;
      if (finalJets.at (1).E () < 20) continue ;
      if (fabs (finalJets.at (0).Eta ()) > 6.5) continue ;
      if (fabs (finalJets.at (1).Eta ()) > 6.5) continue ;

      sort (finalLeptons.rbegin (), finalLeptons.rend (), ptSort ()) ;
      if (finalLeptons.at (0).Pt () < 20) continue ;
      if (finalLeptons.at (1).Pt () < 20) continue ;
      if (finalLeptons.at (0).E () < 20) continue ;
      if (finalLeptons.at (1).E () < 20) continue ;
      if (fabs (finalLeptons.at (0).Eta ()) > 4) continue ;
      if (fabs (finalLeptons.at (1).Eta ()) > 4) continue ;

      float mjj = (finalJets.at (0) + finalJets.at (1)).M () ;
      if (mjj < 300) continue ;
      if (fabs (finalJets.at (0).Eta () - finalJets.at (1).Eta ()) < 2) continue ;
      float m4l = (finalNeutrinos.at (0) + finalLeptons.at (0) + finalLeptons.at (1) + finalNeutrinos.at (0)).M () ;
      if (m4l < 130) continue ;

      ++selected ;      
      h_NJ->Fill (finalJets.size ()) ;
//      cout << "Njs = " << finalJets.size () << "\n" ;
//      cout << "pts: " 
//           << finalJets.at (0).Pt () 
//           << "\t"
//           << finalJets.at (1).Pt () 
//           << "\n" ;

      for (int i = 0 ; i < initialParticlesPDGId.size () ; ++i)
        h_initPart->Fill (initialParticlesPDGId.at (i)) ;
      
      h_vbf0_eta->Fill (finalJets.at (0).Eta ()) ;            
      h_vbf0_phi->Fill (finalJets.at (0).Phi ()) ;            
      h_vbf0_pt-> Fill (finalJets.at (0).Pt ()) ;        

      h_vbf1_eta->Fill (finalJets.at (1).Eta ()) ;            
      h_vbf1_phi->Fill (finalJets.at (1).Phi ()) ;            
      h_vbf1_pt-> Fill (finalJets.at (1).Pt ()) ;        

      h_mjj_vbf->Fill (mjj) ;
      h_m4l->Fill (m4l) ;
      
    } // loop over events
  cout << "total eff  : " << selected * 1./ieve << endl ;
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  int N_EWK = 20000 ;
  float XS_EWK = 4.13649215685881443 ; /* fb */
  map<string, TH1F *> hmap_EWK = 
    readSample ("/Users/govoni/data/TP/forWWInterference/EWK.lhe", "EWK", N_EWK) ;

  int N_QCD = 20000 ;
  float XS_QCD = 1.06691296353271774 ; /* fb */
  map<string, TH1F *> hmap_QCD = 
    readSample ("/Users/govoni/data/TP/forWWInterference/QCD.lhe", "QCD", N_QCD) ;

  int N_ALL = 20000 ;
  float XS_ALL = 5.38030176830107346 ; /* fb */
  map<string, TH1F *> hmap_ALL = 
    readSample ("/Users/govoni/data/TP/forWWInterference/ALL.lhe", "ALL", N_ALL) ;


  TFile outfile ("testWWInterf.root", "recreate") ;
  savemap (hmap_EWK,  outfile,  XS_EWK * 1./N_EWK) ;  /*result in fb */
  savemap (hmap_QCD,  outfile,  XS_QCD * 1./N_QCD) ;  /*result in fb */
  savemap (hmap_ALL,  outfile,  XS_ALL * 1./N_ALL) ;  /*result in fb */
  
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
  
  map<string, TH1F *>::iterator iMap_QCD = hmap_QCD.begin () ;
  map<string, TH1F *>::iterator iMap_ALL = hmap_ALL.begin () ;

  // plotting
  for (map<string, TH1F *>::iterator iMap_EWK = hmap_EWK.begin () ;
       iMap_EWK != hmap_EWK.end () ;
       ++iMap_EWK)
    {
      iMap_EWK->second->SetStats (0) ;
      iMap_QCD->second->SetStats (0) ;
      iMap_ALL->second->SetStats (0) ;
      
      iMap_EWK->second->SetTitle ("") ;
      iMap_QCD->second->SetTitle ("") ;
      iMap_ALL->second->SetTitle ("") ;
      
      iMap_EWK->second->SetLineWidth (2) ;
      iMap_QCD->second->SetLineWidth (2) ;
      iMap_ALL->second->SetLineWidth (4) ;
      
      iMap_EWK->second->SetLineStyle (2) ;
      iMap_QCD->second->SetLineStyle (3) ;
      iMap_ALL->second->SetLineColor (50) ;
      
      TString name = "adding_" ;
      name += iMap_EWK->second->GetName () ;
      TH1F * dummy = (TH1F *) iMap_EWK->second->Clone (name) ;
      dummy->Add (iMap_QCD->second) ;
      dummy->SetFillColor (38) ;
      dummy->SetLineColor (38) ;
      dummy->SetLineStyle (1) ;

      TLegend leg (0.3, 0.9, 0.8, 1) ;
      leg.SetNColumns (3) ;
      leg.SetLineStyle (0) ;
      leg.SetFillStyle (0) ;
      leg.AddEntry (iMap_EWK->second, "EWK", "l") ;
      leg.AddEntry (iMap_QCD->second, "QCD", "l") ;
      leg.AddEntry (iMap_ALL->second, "ALL", "l") ;
      leg.AddEntry (dummy, "EWK+QCD", "fl") ;
      
      iMap_ALL->second->GetXaxis ()->SetTitle (iMap_EWK->first.c_str ()) ; 
      dummy->GetXaxis ()->SetTitle (iMap_EWK->first.c_str ()) ;        

      iMap_ALL->second->Draw ("histo") ;           
      dummy->Draw ("histo same") ;      
      iMap_ALL->second->Draw ("histo same") ;           
      iMap_EWK->second->Draw ("histo same") ;           
      iMap_QCD->second->Draw ("histo same") ;           

      leg.Draw () ;
      c1.Print (("WWInterf_" + iMap_EWK->first + ".png").c_str (), "png") ;

      ++iMap_QCD ;
      ++iMap_ALL ;
    }   
  
  return 0 ;
}


