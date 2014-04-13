#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
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
  c++ -o lookAtMWWSpectrum `root-config --glibs --cflags` \
     -lm lookAtMWWSpectrum.cpp
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
  TH1F * h_vbf0_eta  = addHistoToMap (histos, string ("vbf0_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_vbf1_eta  = addHistoToMap (histos, string ("vbf1_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_cjet0_eta = addHistoToMap (histos, string ("cjet0_eta_") + radice, 50, 0, 6) ;
  TH1F * h_cjet1_eta = addHistoToMap (histos, string ("cjet1_eta_") + radice, 50, 0, 6) ;
  TH1F * h_lep_eta   = addHistoToMap (histos, string ("lep_eta_")   + radice, 50, 0, 6) ;
  TH1F * h_mjj_cen   = addHistoToMap (histos, string ("mjj_cen_")   + radice, 200, 50, 150) ;
  TH1F * h_mjj_vbf   = addHistoToMap (histos, string ("mjj_vbf_")   + radice, 50, 0, 5000) ;

  TH1F * h_cut_vbf0_eta  = addHistoToMap (histos, string ("cut_vbf0_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_cut_vbf1_eta  = addHistoToMap (histos, string ("cut_vbf1_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_cut_cjet0_eta = addHistoToMap (histos, string ("cut_cjet0_eta_") + radice, 50, 0, 6) ;
  TH1F * h_cut_cjet1_eta = addHistoToMap (histos, string ("cut_cjet1_eta_") + radice, 50, 0, 6) ;
  TH1F * h_cut_lep_eta   = addHistoToMap (histos, string ("cut_lep_eta_")   + radice, 50, 0, 6) ;
  TH1F * h_cut_mjj_cen   = addHistoToMap (histos, string ("cut_mjj_cen_")   + radice, 200, 50, 150) ;
  TH1F * h_cut_mjj_vbf   = addHistoToMap (histos, string ("cut_mjj_vbf_")   + radice, 50, 0, 5000) ;
    
  int ieve = 0 ;
  // loop over events
  while ( reader.readEvent () ) 
    {
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
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

      // get the two jets with the W invariant mass
      pair<int, int> Wpair = findPairWithWMass (finalQuarks) ;

      // get the tag jets
      pair<int, int> tagpair (0, 1) ;
      for (int i = 0 ; i < 4 ; ++i)
        {
          if (i == Wpair.first || i == Wpair.second) continue ;
          else
            {
              tagpair.first = i ;
              for (int j = i + 1 ; j < 4 ; ++j)
                {
                  if (j == Wpair.first || j == Wpair.second) continue ;
                  else
                    {
                      tagpair.second = j ;
                      break ;
                    }
                }
              break ;   
            }
        }

      if (finalQuarks.at (tagpair.first).Pt () < finalQuarks.at (tagpair.second).Pt ())  
        {
          swap (tagpair.first, tagpair.second) ;
        }
            
      // get the WW system mass
      TLorentzVector total = (leptons.at (0) + neutrinos.at (0)) + 
                             (finalQuarks.at (Wpair.first) + finalQuarks.at (Wpair.second)) ;

      float mWjj = (finalQuarks.at (Wpair.first) + finalQuarks.at (Wpair.second)).M () ;
      h_mjj_cen->Fill (mWjj) ;

      float mjj = (finalQuarks.at (tagpair.first) + finalQuarks.at (tagpair.second)).M () ;
      h_mjj_vbf->Fill (mjj) ;

      if ( 6 != Wpair.first + Wpair.second + tagpair.first + tagpair.second)
        {
          cout << "WARNING strange jet assignment, skipping event : " ;
          cout << Wpair.first << " " << Wpair.second << " " 
               << tagpair.first << " " << tagpair.second << "\n" ;
          continue ;
        }  

      // put some cuts
      // ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
      
      if (leptons.at (0).Pt () < 30 /*GeV */) continue ; 
      if (neutrinos.at (0).Pt () < 30 /*GeV */) continue ;
      if (finalQuarks.at (Wpair.first).Pt () < 30 /*GeV */) continue ;
      if (finalQuarks.at (Wpair.second).Pt () < 30 /*GeV */) continue ;
      if (finalQuarks.at (tagpair.first).Pt () < 30 /*GeV */) continue ;
      if (finalQuarks.at (tagpair.second).Pt () < 30 /*GeV */) continue ;

      h_cut_mjj_cen->Fill (mWjj) ;
      h_cut_mjj_vbf->Fill (mjj) ;

    } // loop over events
    
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  int Ntot = 100000 ;
  std::ifstream ifs (argv[1]) ;
  map<string, TH1F *> hmap_8Tev = 
    readSample ("/Users/govoni/data/lvjj_samples/lookAtVBF/total.8TeV.lhe", "8TeV", Ntot) ;
  map<string, TH1F *> hmap_13Tev = 
    readSample ("/Users/govoni/data/lvjj_samples/lookAtVBF/total.13TeV.lhe", "13TeV", Ntot) ;

  TFile outfile ("output.root", "recreate") ;
  savemap (hmap_8Tev,  outfile,  76. * 2 /* fb */ / Ntot) ; // the multiplication accounts for muons and electrons
  savemap (hmap_13Tev, outfile, 262. * 2 /* fb */ / Ntot) ; // the multiplication accounts for muons and electrons
  
  // Now we are done.
  return 0 ;
}


/* 


TFile *_file0 = TFile::Open("output.root")

TCanvas c1 
c1.SetLogy ()

mjj_vbf_8TeV->SetLineColor (kBlue + 2)
mjj_vbf_13TeV->SetLineWidth (2)
mjj_vbf_13TeV->Draw ()
mjj_vbf_8TeV->Draw ("same")
cut_mjj_vbf_8TeV->SetLineColor (kBlue + 2)
cut_mjj_vbf_13TeV->SetLineWidth (2)
cut_mjj_vbf_8TeV->SetLineStyle (2)
cut_mjj_vbf_13TeV->SetLineStyle (2)
cut_mjj_vbf_13TeV->Draw ("same")
cut_mjj_vbf_8TeV->Draw ("same")

c1.Print ("ERC_XS.pdf", "pdf") 

mjj_vbf_13TeV->Scale (100)
mjj_vbf_8TeV->Scale (20)
mjj_vbf_13TeV->Draw ()
mjj_vbf_8TeV->Draw ("same")
cut_mjj_vbf_13TeV->Scale (100)
cut_mjj_vbf_8TeV->Scale (20)
cut_mjj_vbf_13TeV->Draw ("same")
cut_mjj_vbf_8TeV->Draw ("same")

c1.Print ("ERC_XS_lumi.pdf", "pdf") 




*/