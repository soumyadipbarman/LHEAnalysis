// c++ -o Zeejets -lm `root-config --cflags --libs --glibs` Zeejets.cpp 

#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>


#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TText.h"
#include "TLegend.h"


//  ------------------------------------------------------------


double deltaPhi (double phi1, double phi2) ;
int selector (TChain * tree,  double weight) ;


TH1F * h_nJetsInside_gt20 ;


//  ========== M A I N    P R O G R A M =========================


int main (int argc, char *argv[])
{ 

  h_nJetsInside_gt20 = new TH1F ("nJetsInside_gt20","nJetsInside_gt20",20,0,20) ;

// ......................................................................

  
  TChain * chain_Zee_15_20 = new TChain ("VBFSimpleTree") ;
  chain_Zee_15_20->Add ("/Users/govoni/data/preselectedHWW/Zee_15_20/NTUPLE*_output_*.root");
  selector (chain_Zee_15_20, 353.1 * 0.293118) ;

  TChain * chain_Zee_20_30 = new TChain ("VBFSimpleTree") ;
  chain_Zee_20_30->Add ("/Users/govoni/data/preselectedHWW/Zee_20_30/NTUPLE*_output_*.root");
  selector (chain_Zee_20_30, 326.7 * 0.342509) ;
  
  TChain * chain_Zee_30_50 = new TChain ("VBFSimpleTree") ;
  chain_Zee_30_50->Add ("/Users/govoni/data/preselectedHWW/Zee_30_50/NTUPLE*_output_*.root");
  selector (chain_Zee_30_50, 227.0 * 0.418665) ;
  
  TChain * chain_Zee_50_80 = new TChain ("VBFSimpleTree") ;
  chain_Zee_50_80->Add ("/Users/govoni/data/preselectedHWW/Zee_50_80/NTUPLE*_output_*.root");
  selector (chain_Zee_50_80, 93.17 * 0.50172) ;
  
  TChain * chain_Zee_80_120 = new TChain ("VBFSimpleTree") ;
  chain_Zee_80_120->Add ("/Users/govoni/data/preselectedHWW/Zee_80_120/NTUPLE*_output_*.root");
  selector (chain_Zee_80_120, 31.48 * 0.57292) ;
  
  TChain * chain_Zee_120_170 = new TChain ("VBFSimpleTree") ;
  chain_Zee_120_170->Add ("/Users/govoni/data/preselectedHWW/Zee_120_170/NTUPLE*_output_*.root");
  selector (chain_Zee_120_170, 9.63 * 0.622112) ;
  
  TChain * chain_Zee_170_230 = new TChain ("VBFSimpleTree") ;
  chain_Zee_170_230->Add ("/Users/govoni/data/preselectedHWW/Zee_170_230/NTUPLE*_output_*.root");
  selector (chain_Zee_170_230, 2.92 * 0.67604) ;

  TFile output ("outputZeejets.root","recreate") ;
  output.cd () ;
  h_nJetsInside_gt20->Write () ;
  output.Close () ;

  return 0 ;
}


//  ------------------------------------------------------------


double 
deltaPhi (double phi1, double phi2)
  { 
    double deltaphi = fabs (phi1 - phi2) ;
    if (deltaphi > 6.283185308) deltaphi -= 6.283185308 ;
    if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi ;
    return deltaphi ;
  }


//  ------------------------------------------------------------


//!PG main function
int 
selector (TChain * tree, double weight)
{
  TClonesArray * tagJets = new TClonesArray ("TLorentzVector") ;
  tree->SetBranchAddress ("tagJets", &tagJets) ;
  TClonesArray * otherJets = new TClonesArray ("TLorentzVector") ;
  tree->SetBranchAddress ("otherJets", &otherJets) ;
  
  int nentries = (int) tree->GetEntries () ;

  //PG loop over the events
  //  int nentries = 100 ;
  for (int evt = 0 ; evt < nentries ; ++evt)
    {
      tree->GetEntry (evt) ;
      int cutId = 0 ;

      if (tagJets->GetEntries () != 2) continue ; 
      TLorentzVector * primoTAG = (TLorentzVector*) (tagJets->At (0)) ;
      TLorentzVector * secondoTAG = (TLorentzVector*) (tagJets->At (1)) ; 
      //PG get the first two in pt
      if (primoTAG->Eta () > secondoTAG->Eta ()) std::swap (primoTAG, secondoTAG) ;

      int ojetsNum = 0 ;
      for (int ojetIt = 0 ; ojetIt < otherJets->GetEntries () ; ++ojetIt)
        {
          if ( ((TLorentzVector*) (otherJets->At (ojetIt)))->Pt () < 20 /*GeV*/) continue ;
          if ( ((TLorentzVector*) (otherJets->At (ojetIt)))->Eta () < primoTAG->Eta () ||
               ((TLorentzVector*) (otherJets->At (ojetIt)))->Eta () > secondoTAG->Eta ()) continue ;
          ++ojetsNum ;        
        } //PG loop over other jets
      h_nJetsInside_gt20->Fill (ojetsNum,weight) ;

    } //PG loop over the events

  delete tagJets  ;  
  delete otherJets  ;

  return 0;
  
}
