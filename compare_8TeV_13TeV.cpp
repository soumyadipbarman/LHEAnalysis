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
  c++ -o compare_8TeV_13TeV `root-config --glibs --cflags` \
     -lm compare_8TeV_13TeV.cpp

FIXME e' bacato, considera anche i neutrini come getti!!!!
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


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


TH1F *
addHistoToMap (map<string, TH1F *> & hmap, string name, int bin, float min, float max)
{
  TH1F * dummy = new TH1F (name.c_str (), name.c_str (), bin, min, max) ;
  hmap[name] = dummy ;
  return dummy ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


void 
savemap (map<string, TH1F *> & hmap, TFile & outfile)
{
  outfile.cd () ;
  for (map<string, TH1F *>::iterator iMap = hmap.begin () ;
      iMap != hmap.end () ; ++ iMap)
    {
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
  TH1F * h_vbf0_eta = addHistoToMap (histos, string ("vbf0_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_vbf1_eta = addHistoToMap (histos, string ("vbf1_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_cjet0_eta = addHistoToMap (histos, string ("cjet0_eta_") + radice, 50, 0, 6) ;
  TH1F * h_cjet1_eta = addHistoToMap (histos, string ("cjet1_eta_") + radice, 50, 0, 6) ;
  TH1F * h_lep_eta = addHistoToMap (histos, string ("lep_eta_")   + radice, 50, 0, 6) ;
  TH1F * h_mjj_cen = addHistoToMap (histos, string ("mjj_cen_")   + radice, 200, 0, 200) ;

  TH1F * h_mjj_vbf = addHistoToMap (histos, string ("mjj_vbf_")   + radice, 200, 0, 200) ;
  
  TH1F * h_cut_vbf0_eta = addHistoToMap (histos, string ("cut_vbf0_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_cut_vbf1_eta = addHistoToMap (histos, string ("cut_vbf1_eta_")  + radice, 50, 0, 6) ;
  TH1F * h_cut_cjet0_eta = addHistoToMap (histos, string ("cut_cjet0_eta_") + radice, 50, 0, 6) ;
  TH1F * h_cut_cjet1_eta = addHistoToMap (histos, string ("cut_cjet1_eta_") + radice, 50, 0, 6) ;
  TH1F * h_cut_lep_eta = addHistoToMap (histos, string ("cut_lep_eta_")   + radice, 50, 0, 6) ;
  TH1F * h_cut_mjj_cen = addHistoToMap (histos, string ("cut_mjj_cen_")   + radice, 200, 0, 200) ;
  
  int ieve = 0 ;
  // loop over events
  while ( reader.readEvent () ) 
    {
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
      if (maxevents > 0 && ieve > maxevents) break ;
      ++ieve;
  
      vector<TLorentzVector> leptons ;      
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
                  abs (reader.hepeup.IDUP.at (iPart)) == 13)   //PG muon
                {
                  leptons.push_back (TLorentzVector
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3) //PG E
                    )) ;
                } //PG leptons
              else
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

      // sort jets in eta
      sort (finalQuarks.begin (), finalQuarks.end (), etaSort ()) ;
      
      // forward jets: the first number is the most forward jet
      float eta_tj_0 = finalQuarks.at (0).Eta () ;
      float eta_tj_1 = finalQuarks.back ().Eta () ;

      if (fabs (eta_tj_0) < fabs (eta_tj_1)) swap (eta_tj_0, eta_tj_1) ;

      h_vbf0_eta->Fill (eta_tj_0) ;
      h_vbf1_eta->Fill (eta_tj_1) ;
      
      // central jets: the first number is the one with the largest pT
      int index_0 = 1 ;
      int index_1 = 2 ;
      if (finalQuarks.at (index_0).Pt () < finalQuarks.at (index_1).Pt ()) swap (index_0, index_1) ;
      
      h_cjet0_eta->Fill (finalQuarks.at (index_0).Eta ()) ;
      h_cjet1_eta->Fill (finalQuarks.at (index_1).Eta ()) ;

      h_mjj_cen->Fill ((finalQuarks.at (index_0) + finalQuarks.at (index_1)).M ()) ;
      
      // lepton
      h_lep_eta->Fill (leptons.at (0).Eta ()) ;

      float mjj = (finalQuarks.at (0) + finalQuarks.at (3)).M () ;
      h_mjj_vbf->Fill (mjj) ;

      if (mjj < 500) continue ;

      h_cut_vbf0_eta->Fill (eta_tj_0) ;
      h_cut_vbf1_eta->Fill (eta_tj_1) ;
      h_cut_cjet0_eta->Fill (finalQuarks.at (index_0).Eta ()) ;
      h_cut_cjet1_eta->Fill (finalQuarks.at (index_1).Eta ()) ;
      h_cut_mjj_cen->Fill ((finalQuarks.at (index_0) + finalQuarks.at (index_1)).M ()) ;
      h_cut_lep_eta->Fill (leptons.at (0).Eta ()) ;

    } // loop over events
    
  return histos ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  gROOT->SetStyle ("Plain") ;

  std::ifstream ifs (argv[1]) ;
  map<string, TH1F *> hmap_8Tev = 
    readSample ("/Users/govoni/data/lvjj_samples/lookAtVBF/total.8TeV.lhe", "8TeV", 100000) ;
  map<string, TH1F *> hmap_13Tev = 
    readSample ("/Users/govoni/data/lvjj_samples/lookAtVBF/total.13TeV.lhe", "13TeV", 100000) ;

  TFile outfile ("output.root", "recreate") ;
  savemap (hmap_8Tev, outfile) ;
  savemap (hmap_13Tev, outfile) ;
  

  // Now we are done.
  return 0 ;
}


/* 


TCanvas c1 
c1.SetLogy ()

lep_eta_8TeV->SetFillColor (kGray)
lep_eta_8TeV->SetLineColor (kGray+1)

lep_eta_13TeV->SetLineColor (kBlue)
lep_eta_13TeV->SetLineWidth (2)

lep_eta_8TeV->SetStats (0)

lep_eta_8TeV->Draw ()
lep_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (lep_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (lep_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("lep_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

vbf0_eta_8TeV->SetFillColor (kGray)
vbf0_eta_8TeV->SetLineColor (kGray+1)

vbf0_eta_13TeV->SetLineColor (kBlue)
vbf0_eta_13TeV->SetLineWidth (2)

vbf0_eta_8TeV->SetStats (0)

vbf0_eta_8TeV->Draw ()
vbf0_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (vbf0_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (vbf0_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("vbf0_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

vbf1_eta_8TeV->SetFillColor (kGray)
vbf1_eta_8TeV->SetLineColor (kGray+1)

vbf1_eta_13TeV->SetLineColor (kBlue)
vbf1_eta_13TeV->SetLineWidth (2)

vbf1_eta_8TeV->SetStats (0)

vbf1_eta_8TeV->Draw ()
vbf1_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (vbf1_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (vbf1_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("vbf1_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

cjet0_eta_8TeV->SetFillColor (kGray)
cjet0_eta_8TeV->SetLineColor (kGray+1)

cjet0_eta_13TeV->SetLineColor (kBlue)
cjet0_eta_13TeV->SetLineWidth (2)

cjet0_eta_8TeV->SetStats (0)

cjet0_eta_8TeV->Draw ()
cjet0_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cjet0_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cjet0_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cjet0_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

cjet1_eta_8TeV->SetFillColor (kGray)
cjet1_eta_8TeV->SetLineColor (kGray+1)

cjet1_eta_13TeV->SetLineColor (kBlue)
cjet1_eta_13TeV->SetLineWidth (2)

cjet1_eta_8TeV->SetStats (0)

cjet1_eta_8TeV->Draw ()
cjet1_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cjet1_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cjet1_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cjet1_eta.pdf", "pdf")


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

mjj_cen_8TeV->SetFillColor (kGray)
mjj_cen_8TeV->SetLineColor (kGray+1)

mjj_cen_13TeV->SetLineColor (kBlue)
mjj_cen_13TeV->SetLineWidth (2)

mjj_cen_8TeV->SetStats (0)

mjj_cen_8TeV->Draw ()
mjj_cen_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (mjj_cen_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (mjj_cen_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("mjj_cen.pdf", "pdf")


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


mjj_vbf_8TeV->SetFillColor (kGray)
mjj_vbf_8TeV->SetLineColor (kGray+1)

mjj_vbf_13TeV->SetLineColor (kBlue)
mjj_vbf_13TeV->SetLineWidth (2)

mjj_vbf_8TeV->SetStats (0)

mjj_vbf_8TeV->Draw ()
mjj_vbf_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (mjj_vbf_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (mjj_vbf_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("mjj_vbf.pdf", "pdf")


// ====================================================================



TCanvas c1 
c1.SetLogy ()

cut_lep_eta_8TeV->SetFillColor (kGray)
cut_lep_eta_8TeV->SetLineColor (kGray+1)

cut_lep_eta_13TeV->SetLineColor (kBlue)
cut_lep_eta_13TeV->SetLineWidth (2)

cut_lep_eta_8TeV->SetStats (0)

cut_lep_eta_8TeV->Draw ()
cut_lep_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cut_lep_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cut_lep_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cut_lep_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

cut_vbf0_eta_8TeV->SetFillColor (kGray)
cut_vbf0_eta_8TeV->SetLineColor (kGray+1)

cut_vbf0_eta_13TeV->SetLineColor (kBlue)
cut_vbf0_eta_13TeV->SetLineWidth (2)

cut_vbf0_eta_8TeV->SetStats (0)

cut_vbf0_eta_8TeV->Draw ()
cut_vbf0_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cut_vbf0_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cut_vbf0_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cut_vbf0_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

cut_vbf1_eta_8TeV->SetFillColor (kGray)
cut_vbf1_eta_8TeV->SetLineColor (kGray+1)

cut_vbf1_eta_13TeV->SetLineColor (kBlue)
cut_vbf1_eta_13TeV->SetLineWidth (2)

cut_vbf1_eta_8TeV->SetStats (0)

cut_vbf1_eta_8TeV->Draw ()
cut_vbf1_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cut_vbf1_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cut_vbf1_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cut_vbf1_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

cut_cjet0_eta_8TeV->SetFillColor (kGray)
cut_cjet0_eta_8TeV->SetLineColor (kGray+1)

cut_cjet0_eta_13TeV->SetLineColor (kBlue)
cut_cjet0_eta_13TeV->SetLineWidth (2)

cut_cjet0_eta_8TeV->SetStats (0)

cut_cjet0_eta_8TeV->Draw ()
cut_cjet0_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cut_cjet0_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cut_cjet0_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cut_cjet0_eta.pdf", "pdf")

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

cut_cjet1_eta_8TeV->SetFillColor (kGray)
cut_cjet1_eta_8TeV->SetLineColor (kGray+1)

cut_cjet1_eta_13TeV->SetLineColor (kBlue)
cut_cjet1_eta_13TeV->SetLineWidth (2)

cut_cjet1_eta_8TeV->SetStats (0)

cut_cjet1_eta_8TeV->Draw ()
cut_cjet1_eta_13TeV->Draw ("same")

c1_leg = new TLegend (0.7, 0.7, 0.85, 0.85) ;
c1_leg->SetFillStyle (0) ;
c1_leg->SetBorderSize (0) ;
c1_leg->SetTextFont (42) ;
c1_leg->AddEntry (cut_cjet1_eta_8TeV, " 8 TeV","f") ;
c1_leg->AddEntry (cut_cjet1_eta_13TeV, "13 TeV","l") ;
c1_leg->Draw () ;

c1.Print ("cut_cjet1_eta.pdf", "pdf")




*/