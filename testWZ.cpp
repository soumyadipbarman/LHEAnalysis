/*
  c++ -o testWZ `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags` \
     -lm testWZ.cpp
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
//#include "hFactory.h"
//#include "h2Factory.h"
//#include "hFunctions.h"

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

#include "LHAPDF/LHAPDF.h"


using namespace std ;
using namespace LHAPDF ;


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
readSample (string sampleName, string radice, int maxevents = -1, 
            bool doPdfReweight = false)
{
  cout << "reading " << sampleName << endl ;
  std::ifstream ifs (sampleName.c_str ()) ;
  LHEF::Reader reader (ifs) ;
  
  map<string, TH1F *> histos ;

  TH1F * h_pdfWeight       = addHistoToMap (histos, string ("pdfWeight_")        + radice, 200, 0., 2.) ;

  TH1F * h_initPart        = addHistoToMap (histos, string ("initPart_")         + radice, 135, -22.5, 22.5) ;
  TH1F * h_finalJetsPDGId  = addHistoToMap (histos, string ("finalJetsPDGId_")   + radice, 135, -22.5, 22.5) ;

  TH1F * h_lep0_eta        = addHistoToMap (histos, string ("lep0_eta_")         + radice, 40, -6, 6) ;
  TH1F * h_lep0_pt         = addHistoToMap (histos, string ("lep0_pt_")          + radice, 40, 0, 400) ;
  TH1F * h_lep0_phi        = addHistoToMap (histos, string ("lep0_phi_")         + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_lep1_eta        = addHistoToMap (histos, string ("lep1_eta_")         + radice, 40, -6, 6) ;
  TH1F * h_lep1_pt         = addHistoToMap (histos, string ("lep1_pt_")          + radice, 40, 0, 400) ;
  TH1F * h_lep1_phi        = addHistoToMap (histos, string ("lep1_phi_")         + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_lep2_eta        = addHistoToMap (histos, string ("lep2_eta_")         + radice, 40, -6, 6) ;
  TH1F * h_lep2_pt         = addHistoToMap (histos, string ("lep2_pt_")          + radice, 40, 0, 400) ;
  TH1F * h_lep2_phi        = addHistoToMap (histos, string ("lep2_phi_")         + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_met_eta         = addHistoToMap (histos, string ("met_eta_")          + radice, 40, -6, 6) ;
  TH1F * h_met_pt          = addHistoToMap (histos, string ("met_pt_")           + radice, 40, 0, 400) ;
  TH1F * h_met_phi         = addHistoToMap (histos, string ("met_phi_")          + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_vbf0_eta        = addHistoToMap (histos, string ("vbf0_eta_")         + radice, 40, -6, 6) ;
  TH1F * h_vbf0_pt         = addHistoToMap (histos, string ("vbf0_pt_")          + radice, 40, 0, 400) ;
  TH1F * h_vbf0_phi        = addHistoToMap (histos, string ("vbf0_phi_")         + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_vbf1_eta        = addHistoToMap (histos, string ("vbf1_eta_")         + radice, 40, -6, 6) ;
  TH1F * h_vbf1_pt         = addHistoToMap (histos, string ("vbf1_pt_")          + radice, 40, 0, 400) ;
  TH1F * h_vbf1_phi        = addHistoToMap (histos, string ("vbf1_phi_")         + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_bjet0_eta       = addHistoToMap (histos, string ("bjet0_eta_")        + radice, 40, -6, 6) ;
  TH1F * h_bjet0_pt        = addHistoToMap (histos, string ("bjet0_pt_")         + radice, 40, 0, 400) ;
  TH1F * h_bjet0_phi       = addHistoToMap (histos, string ("bjet0_phi_")        + radice, 30, -3.14, 3.14) ;
                                                                                 
  TH1F * h_Ljet0_eta       = addHistoToMap (histos, string ("Ljet0_eta_")        + radice, 40, -6, 6) ;
  TH1F * h_Ljet0_pt        = addHistoToMap (histos, string ("Ljet0_pt_")         + radice, 40, 0, 400) ;
  TH1F * h_Ljet0_phi       = addHistoToMap (histos, string ("Ljet0_phi_")        + radice, 30, -3.14, 3.14) ;

  TH1F * h_detajj_vbf      = addHistoToMap (histos, string ("detajj_vbf_")       + radice, 25, 0, 10) ;
  TH1F * h_mjj_vbf         = addHistoToMap (histos, string ("mjj_vbf_")          + radice, 25, 0, 2500) ;
  TH1F * h_NJ              = addHistoToMap (histos, string ("NJ_")               + radice, 5, 0, 5) ;
                                                                                 
  TH1F * h_m4l             = addHistoToMap (histos, string ("m4l_")              + radice, 35, 0, 2275) ;
  TH1F * h_mtop            = addHistoToMap (histos, string ("mtop_")             + radice, 100, 0, 2000) ;
  TH1F * h_mtopSel         = addHistoToMap (histos, string ("mtopSel_")          + radice, 100, 0, 2000) ;
  TH1F * h_mllOS           = addHistoToMap (histos, string ("mllOS_")            + radice, 500, 4, 500) ;

  TH1F * h_NnotZleptons    = addHistoToMap (histos, string ("NnotZleptons_")     + radice, 15, -0.5, 4.5) ;
    
  int ieve = 0 ;
  int btagged = 0 ;
  int selected = 0 ;
  // loop over events

  double x[2] = {0., 0.} ;
  int flavour[2] = {0, 0} ;

  while ( reader.readEvent () ) 
    {
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
      if (maxevents > 0 && ieve >= maxevents) break ;
      ++ieve;
  
      vector<int> initialParticlesPDGId ;      
      vector<int> finalJetsPDGId ;      
      vector<TLorentzVector> finalLeptons ;      
      vector<int> finalLeptonsID ;      
      vector<TLorentzVector> finalLeptonsNOTfromZ ;      
      vector<TLorentzVector> finalNeutrinos ;      
      vector<TLorentzVector> initialQuarks ;      
      vector<TLorentzVector> finalJets ;      
      vector<TLorentzVector> finalBjets ;      
      vector<TLorentzVector> finalLjets ;      
      
      int leptonsIDproduct = 1 ;

      //PG loop over particles in the event
      //PG and fill the variables of leptons and quarks
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
          //PG incoming particle          
          TLorentzVector dummy (
              reader.hepeup.PUP.at (iPart).at (0), //PG px
              reader.hepeup.PUP.at (iPart).at (1), //PG py
              reader.hepeup.PUP.at (iPart).at (2), //PG pz
              reader.hepeup.PUP.at (iPart).at (3) //PG E
            ) ;

          if (reader.hepeup.ISTUP.at (iPart) == -1)
            {
              initialParticlesPDGId.push_back (reader.hepeup.IDUP.at (iPart)) ;
              initialQuarks.push_back (dummy) ;
              x[iPart] = initialQuarks.back ().P () / 7000. ;
              flavour[iPart] = reader.hepeup.IDUP.at (iPart) ;

            } //PG incoming particle          

          //PG outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              // jets  
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7 ||  // quarks
                  abs (reader.hepeup.IDUP.at (iPart)) == 21 ) // gluons
                {
                  finalJetsPDGId.push_back (reader.hepeup.IDUP.at (iPart)) ;
                  finalJets.push_back (dummy) ;
                  if (abs (reader.hepeup.IDUP.at (iPart)) == 5)  // b's
                    {
                      finalBjets.push_back (dummy) ;
                    } // b's
                  else
                    {
                      finalLjets.push_back (dummy) ;
                    }  
                } // jets
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 13 ) // charged leptons
                {
                  leptonsIDproduct *= reader.hepeup.IDUP.at (iPart) ;
                  finalLeptons.push_back (dummy) ;
                  finalLeptonsID.push_back (reader.hepeup.IDUP.at (iPart)) ;
                  int mother1 = reader.hepeup.MOTHUP.at (iPart).first ;
                  if (fabs (reader.hepeup.IDUP.at (mother1)) != 23)
                    {
                      finalLeptonsNOTfromZ.push_back (dummy) ;
                    }
                } //PG charged leptons
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 14 ) // neutrinos
                {
                  finalNeutrinos.push_back (dummy) ;
                } // neutrinos

            } //PG outgoing particles
        } //PG loop over particles in the event

      if (initialParticlesPDGId.size () != 2) exit (1) ; 

      if (finalJets.size () != 2) continue ;
      if (finalNeutrinos.size () != 1) continue ;
      if (finalLeptons.size () != 3) continue ;

      int totalCharge = 0 ;
      for (int i = 0 ; i < finalLeptonsID.size () ; ++i)
        totalCharge += finalLeptonsID.at (i) / abs (finalLeptonsID.at (i)) ;
        
      TLorentzVector lonelyLepton ;
      vector<TLorentzVector> otherLeptons ;
      for (int i = 0 ; i < finalLeptonsID.size () ; ++i)
        {
          if (totalCharge == finalLeptonsID.at (i) / abs (finalLeptonsID.at (i)))
            {
              otherLeptons.push_back (finalLeptons.at (i)) ;  
            }
          else
            { 
              lonelyLepton = finalLeptons.at (i) ;
            }  
        }

      // 11 * 11 * 11 = 1331
      // 11 * 11 * 13 = 1573
      // 11 * 13 * 13 = 1859
      // 13 * 13 * 13 = 2197
      int prodFlav = 1. ;
      for (int i = 0 ; i < finalLeptonsID.size () ; ++i) prodFlav *= finalLeptonsID.at (i) ;
      //PG discard events with three leptons of the same flavour
      if (abs (prodFlav) == 1331 ||
          abs (prodFlav) == 2197) continue ;

//      h_NnotZleptons->Fill (finalLeptonsNOTfromZ.size ()) ;

//      //PG no b's in the initial state
//      if (fabs (initialParticlesPDGId.at (0)) > 4) continue ;
//      if (fabs (initialParticlesPDGId.at (1)) > 4) continue ;

//      //PG only b's in the initial state   
//      if (fabs (initialParticlesPDGId.at (0)) < 5) continue ;
//      if (fabs (initialParticlesPDGId.at (1)) < 5) continue ;
// (no events in the sample)

//      //PG only one b in the initial state
//      if ((fabs (initialParticlesPDGId.at (0)) > 4 &&
//           fabs (initialParticlesPDGId.at (1)) > 4) ||
//          (fabs (initialParticlesPDGId.at (0)) < 5 &&
//           fabs (initialParticlesPDGId.at (1)) < 5)) continue ;

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

//      bool isSingleTop = false ;
//      if (isSingleTop) continue ;     

      sort (finalLeptons.rbegin (), finalLeptons.rend (), ptSort ()) ;
      if (finalLeptons.at (0).Pt () < 20) continue ;
      if (finalLeptons.at (1).Pt () < 20) continue ;
      if (finalLeptons.at (2).Pt () < 20) continue ;
      if (finalLeptons.at (0).E () < 20) continue ;
      if (finalLeptons.at (1).E () < 20) continue ;
      if (finalLeptons.at (2).E () < 20) continue ;
      if (fabs (finalLeptons.at (0).Eta ()) > 4) continue ;
      if (fabs (finalLeptons.at (1).Eta ()) > 4) continue ;
      if (fabs (finalLeptons.at (2).Eta ()) > 4) continue ;

      float mTop = 0. ;
      if (finalBjets.size () > 0)
        {
//           int iLep = 0 ;
//           for (int iLep = 0 ; iLep < 3 ; ++iLep)
//             {
//               if (isLeptonFromZ.at (iLep) == 0) break ;
//             }
//           TLorentzVector tl_top = finalLeptons.at (iLep) + finalNeutrinos.at (0) + finalBjets.at (0) ;
           TLorentzVector tl_top = finalLeptons.at (0) + finalNeutrinos.at (0) + finalBjets.at (0) ;
           mTop = tl_top.M () ;
//           h_mtop->Fill (tl_top.M ()) ;
//           if (tl_top.M () < (173 + 12) && tl_top.M () > (173 - 12)) 
//             {
//               isSingleTop = true ;
//             }   
        }

      float mjj = (finalJets.at (0) + finalJets.at (1)).M () ;
      if (mjj < 300) continue ;
      if (fabs (finalJets.at (0).Eta () - finalJets.at (1).Eta ()) < 2) continue ;
      float m4l = (finalNeutrinos.at (0) + finalLeptons.at (0) + finalLeptons.at (1) + finalLeptons.at (2)).M () ;
      if (m4l < 130) continue ;

      ++selected ;      
      //PG end of selections

      bool foundZ = false ;      
      for (int i = 0 ; i < otherLeptons.size () ; ++i)
        {
          float M = (otherLeptons.at (i) + lonelyLepton).M () ;
          if (M < 95 && M > 85) foundZ = true ;
        }
      if (!foundZ) continue ;
      
      float weight = 1. ;
      if (doPdfReweight)
        {
          //PG the scale:
          float scale = reader.hepeup.SCALUP ;
    
          //PG determine the scale according to the phantom recipe
          double phantomScale = 80.385 * 80.385 + 
              (finalJets.at (0).Pt () * finalJets.at (0).Pt () +
               finalJets.at (1).Pt () * finalJets.at (1).Pt () +
               finalLeptons.at (0).Pt () * finalLeptons.at (0).Pt () +
               finalLeptons.at (1).Pt () * finalLeptons.at (1).Pt () +
               finalLeptons.at (1).Pt () * finalLeptons.at (1).Pt () +
               finalLeptons.at (2).Pt () * finalLeptons.at (2).Pt () +
               finalNeutrinos.at (0).Pt () * finalNeutrinos.at (0).Pt ()) / 6. ;
          phantomScale = sqrt (phantomScale) ;
    
          //PG calculate the weight to be applied to the event
          weight = LHAPDF::xfx (x[0], phantomScale, flavour[0]) * LHAPDF::xfx (x[1], phantomScale, flavour[1]) /
                   (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;
//          cout << weight << endl ;
        }
      h_pdfWeight->Fill (weight) ;

      for (int i = 0 ; i < otherLeptons.size () ; ++i)
        {
          h_mllOS->Fill ((otherLeptons.at (i) + lonelyLepton).M (), weight) ;
        }

      h_NJ->Fill (finalJets.size (), weight) ;

      for (int i = 0 ; i < initialParticlesPDGId.size () ; ++i)
        h_initPart->Fill (initialParticlesPDGId.at (i), weight) ;
      for (int i = 0 ; i < finalJetsPDGId.size () ; ++i)
        h_finalJetsPDGId->Fill (finalJetsPDGId.at (i), weight) ;
      
      h_vbf0_eta->Fill (finalJets.at (0).Eta (), weight) ;            
      h_vbf0_phi->Fill (finalJets.at (0).Phi (), weight) ;            
      h_vbf0_pt-> Fill (finalJets.at (0).Pt (), weight) ;        

      h_vbf1_eta->Fill (finalJets.at (1).Eta (), weight) ;            
      h_vbf1_phi->Fill (finalJets.at (1).Phi (), weight) ;            
      h_vbf1_pt-> Fill (finalJets.at (1).Pt (), weight) ;        

      if (finalBjets.size () > 0)
        {
          h_bjet0_eta->Fill (finalBjets.at (0).Eta (), weight) ;            
          h_bjet0_phi->Fill (finalBjets.at (0).Phi (), weight) ;            
          h_bjet0_pt-> Fill (finalBjets.at (0).Pt (), weight) ;        
        }

      if (finalLjets.size () > 0)
        {
          h_Ljet0_eta->Fill (finalLjets.at (0).Eta (), weight) ;            
          h_Ljet0_phi->Fill (finalLjets.at (0).Phi (), weight) ;            
          h_Ljet0_pt-> Fill (finalLjets.at (0).Pt (), weight) ;        
        }

      h_lep0_eta->Fill (finalLeptons.at (0).Eta (), weight) ;            
      h_lep0_phi->Fill (finalLeptons.at (0).Phi (), weight) ;            
      h_lep0_pt-> Fill (finalLeptons.at (0).Pt (), weight) ;        

      h_lep1_eta->Fill (finalLeptons.at (1).Eta (), weight) ;            
      h_lep1_phi->Fill (finalLeptons.at (1).Phi (), weight) ;            
      h_lep1_pt-> Fill (finalLeptons.at (1).Pt (), weight) ;        

      h_lep2_eta->Fill (finalLeptons.at (2).Eta (), weight) ;            
      h_lep2_phi->Fill (finalLeptons.at (2).Phi (), weight) ;            
      h_lep2_pt-> Fill (finalLeptons.at (2).Pt (), weight) ;        

      h_met_eta->Fill (finalNeutrinos.at (0).Eta (), weight) ;            
      h_met_phi->Fill (finalNeutrinos.at (0).Phi (), weight) ;            
      h_met_pt-> Fill (finalNeutrinos.at (0).Pt (), weight) ;        

      h_detajj_vbf->Fill (fabs (finalJets.at (0).Eta () - finalJets.at (1).Eta ())) ;
      h_mjj_vbf->Fill (mjj, weight) ;
      h_m4l->Fill (m4l, weight) ;
      
    } // loop over events
  cout << "total eff  : " << selected * 1./ieve << endl ;
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  bool doPdfReweight = false ;
  if (argc > 1) doPdfReweight = true ;

  const int SUBSET = 0;
  const string NAME = "cteq6ll"; //"cteq6l1"

  LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
  const int NUMBER = LHAPDF::numberPDF();

  LHAPDF::initPDF (0) ;

  int N_MG = 10000 ;
  float XS_MG = 660.19 * 0.2 * 0.067 ; /* fb */ 
  map<string, TH1F *> hmap_MG = 
    readSample ("/Users/govoni/data/TP/compareEWK/WZ/madgraph_new.lhe", "MG", N_MG, doPdfReweight) ;

  int N_PH = 79938 ;
  float XS_PH = 7.854 ; /* fb */
  map<string, TH1F *> hmap_PH = 
    readSample ("/Users/govoni/data/TP/compareEWK/WZ/phantom.lhe", "PH", N_PH, doPdfReweight) ;

  TFile outfile ("testWZ.root", "recreate") ;
  savemap (hmap_MG,  outfile,  XS_MG * 1./N_MG) ;  /*result in fb */
  savemap (hmap_PH,  outfile,  XS_PH * 1./N_PH) ;  /*result in fb */
  
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
  
  map<string, TH1F *>::iterator iMap_PH = hmap_PH.begin () ;

  // plotting
  for (map<string, TH1F *>::iterator iMap_MG = hmap_MG.begin () ;
       iMap_MG != hmap_MG.end () ;
       ++iMap_MG)
    {
      iMap_MG->second->SetStats (0) ;
      iMap_PH->second->SetStats (0) ;
      
      iMap_MG->second->SetTitle ("") ;
      iMap_PH->second->SetTitle ("") ;
      
      iMap_MG->second->SetLineWidth (2) ;
      iMap_PH->second->SetLineWidth (2) ;
      
      iMap_MG->second->SetLineColor (30) ;
      iMap_PH->second->SetLineColor (9) ;
      
      iMap_MG->second->SetFillColor (30) ;

      TLegend leg (0.3, 0.9, 0.8, 1) ;
      leg.SetNColumns (3) ;
      leg.SetLineStyle (0) ;
      leg.SetFillStyle (0) ;
      leg.AddEntry (iMap_MG->second, "madgraph", "fl") ;
      leg.AddEntry (iMap_PH->second, "phantom", "fl") ;
      
      iMap_MG->second->GetXaxis ()->SetTitle (iMap_MG->first.c_str ()) ;        
      iMap_PH->second->GetXaxis ()->SetTitle (iMap_PH->first.c_str ()) ;        
      iMap_PH->second->Draw () ;           
      iMap_MG->second->Draw ("same") ;           
      iMap_MG->second->DrawCopy ("histsame") ;   
      iMap_PH->second->SetFillStyle (3001) ;
      iMap_PH->second->SetFillColor (9) ;
      iMap_PH->second->DrawCopy ("E2same") ;           
      iMap_PH->second->SetFillStyle (0) ;
      iMap_PH->second->SetFillColor (0) ;
      iMap_PH->second->Draw ("same") ;           
      leg.Draw () ;
      
      if (TString (iMap_MG->second->GetName ()).Contains ("mllOS_")) c1.SetLogy () ;
      
      c1.Print (("WZ_" + iMap_MG->first + ".png").c_str (), "png") ;
      iMap_PH->second->DrawNormalized ("E2") ;           
      iMap_MG->second->DrawNormalized ("histsame") ;   
      iMap_PH->second->DrawNormalized ("E2same") ;           
      leg.Draw () ;
      c1.Print (("WZ_" + iMap_MG->first + "_norm.png").c_str (), "png") ;

      if (TString (iMap_MG->second->GetName ()).Contains ("mllOS_")) c1.SetLogy (0) ;

      ++iMap_PH ;
    }   
  
  return 0 ;
}


