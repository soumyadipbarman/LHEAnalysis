/*
c++ -o read_01 `root-config --glibs --cflags` -lm read_01.cpp 
./read_01 LHEfile.lhe
*/


#include <cmath>
#include <vector>
#include "LHEF.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <functional>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"



using namespace std ;

double 
deltaPhi (double phi1, double phi2)
{
  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  return deltaphi;
}

struct ptSort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator() (TLorentzVector x, TLorentzVector y)
    {
      return x.Pt () < y.Pt () ;
    }
} ;


double
zepp (double eta, double eta1, double eta2)
{
  return eta - 0.5 * (eta1 + eta2) ;
}


int main (int argc, char **argv) {

  // Open a stream connected to an event file:
  if (argc < 2) exit (1) ;
  std::ifstream ifs(argv[1]);

  // Create the Reader object:
  LHEF::Reader reader(ifs);

  // Copy header and init blocks and write them out.
  if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

  // Print out the header information:
  std::cerr << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cerr << reader.initComments;

  // Print out the beam energies:
  std::cerr << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;

  TH1F h_higgsZ  ("h_higgsZ" , "h_higgsZ" , 100, -6, 6   ) ;
  TH1F h_higgsPt ("h_higgsPt", "h_higgsPt", 50 , 0 , 200 ) ;
  TH1F h_j1Pt    ("h_j1Pt"   , "h_j1Pt"   , 200, 0 , 200 ) ;
  TH1F h_j2Pt    ("h_j2Pt"   , "h_j2Pt"   , 50 , 0 , 200 ) ;
  TH1F h_j3Pt    ("h_j3Pt"   , "h_j3Pt"   , 200, 0 , 50  ) ;
  TH1F h_j3Z     ("h_j3Z"    , "h_j3Z"    , 100, -6, 6   ) ;
  TH1F h_cjPt    ("h_cjPt"   , "h_cjPt"   , 200, 0 , 50  ) ;
  TH1F h_njets   ("h_njets"  , "h_njets"  , 4  , 0 , 4   ) ;
  TH1F h_mjj     ("h_mjj"    , "h_mjj"    , 100, 0 , 1000) ;
  TH1F h_detajj  ("h_detajj" , "h_detajj" , 100, 0 , 6   ) ;
  TH1F h_dphijj  ("h_dphijj" , "h_dphijj" , 50 , 0 , 3.15) ;
  TH2F h_detadphijj  ("h_detadphijj" , "h_detadphijj" , 20, 0 , 6, 10 , 0 , 3.15) ;
  TH2F h_ptj1_ptj2   ("h_ptj1_ptj2" , "h_ptj1_ptj2" , 200 , 0 , 200, 200 , 0 , 200) ;

  // Now loop over all events:
  long ieve = 0;
  while ( reader.readEvent() ) 
    {
      ++ieve;
      if (ieve % 100 == 0) std::cerr << "event " << ieve << "\n" ;
  
      vector<TLorentzVector> v_jets ;
      vector<TLorentzVector> v_quarks ;
      vector<TLorentzVector> v_gluons ;
      vector<TLorentzVector> v_higgs ;
      
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {  
 
          if (reader.hepeup.ISTUP.at (iPart) < 0) continue ;

//          std::cerr << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          int type = abs (reader.hepeup.IDUP.at (iPart)) ;
          TLorentzVector particle 
             (
               reader.hepeup.PUP.at (iPart).at (0), //PG px
               reader.hepeup.PUP.at (iPart).at (1), //PG py
               reader.hepeup.PUP.at (iPart).at (2), //PG pz
               reader.hepeup.PUP.at (iPart).at (3)  //PG E
             ) ;

          if (type == 25) //HIGGS
            v_higgs.push_back (particle) ;
          
          if (type < 9 && type > 0) //HIGGS
            {
              v_quarks.push_back (particle) ;
              v_jets.push_back (particle) ;
            }
          
          if (type == 21) //HIGGS
            {
              v_gluons.push_back (particle) ;
              v_jets.push_back (particle) ;
            }
        } //PG loop over particles in the event

      h_higgsPt.Fill (v_higgs.at (0).Pt (), reader.hepeup.XWGTUP) ;
      h_njets.Fill (v_jets.size (), reader.hepeup.XWGTUP) ;
      if (v_jets.size () > 0) 
        {
          sort (v_jets.rbegin (), v_jets.rend (), ptSort ()) ;
          h_j1Pt.Fill (v_jets.at (0).Pt (), reader.hepeup.XWGTUP) ;
          if (v_jets.size () == 1) continue ; 

          h_j2Pt.Fill (v_jets.at (1).Pt (), reader.hepeup.XWGTUP) ;
          h_ptj1_ptj2.Fill (v_jets.at (0).Pt (), v_jets.at (1).Pt (), reader.hepeup.XWGTUP) ;

          if (v_jets.size () == 2) continue ;


          h_j3Pt.Fill (v_jets.at (2).Eta (), reader.hepeup.XWGTUP) ;
          double mineta = v_jets.at (0).Eta () ;
          double maxeta = v_jets.at (1).Eta () ;
          TLorentzVector sumjj = v_jets.at (0) + v_jets.at (1) ;
          h_mjj.Fill (sumjj.M ()) ;
          h_dphijj.Fill (deltaPhi (v_jets.at (0).Phi (), v_jets.at (1).Phi ()), reader.hepeup.XWGTUP) ;
          if (mineta > maxeta) swap (mineta, maxeta) ;
          h_detajj.Fill (maxeta - mineta, reader.hepeup.XWGTUP) ;
          h_detadphijj.Fill (maxeta - mineta, deltaPhi (v_jets.at (0).Phi (), v_jets.at (1).Phi ()), reader.hepeup.XWGTUP) ;
          if (v_jets.at (2).Eta () > mineta && v_jets.at (2).Eta () < maxeta)
              h_cjPt.Fill (v_jets.at (2).Eta (), reader.hepeup.XWGTUP) ;
          h_j3Z.Fill (zepp (v_jets.at (2).Eta () , mineta, maxeta), reader.hepeup.XWGTUP) ;
          h_higgsZ.Fill (zepp (v_higgs.at (0).Eta () , mineta, maxeta), reader.hepeup.XWGTUP) ;
        }

    } // Now loop over all events


  TFile out ("read_01.root", "recreate") ;
  out.cd () ;
  h_higgsPt.Write () ;
  h_j1Pt.Write () ;
  h_j2Pt.Write () ;
  h_j3Pt.Write () ;
  h_cjPt.Write () ;
  h_njets.Write () ;
  h_mjj.Write () ;
  h_detajj.Write () ;
  h_dphijj.Write () ;
  h_j3Z.Write () ;
  h_higgsZ.Write () ;
  h_detadphijj.Write () ;
  h_ptj1_ptj2.Write () ;
  out.Close () ;

  // Now we are done.
  return 0 ;
}

/*

TCanvas c1
h_j1Pt->SetStats (0)
h_j1Pt->SetTitle ("")
h_j1Pt->SetLineWidth (2)
h_j1Pt->Draw ()
TCanvas c2
gStyle->SetPalette (1)
h_ptj1_ptj2->SetStats (0)
h_ptj1_ptj2->SetTitle ("")
h_ptj1_ptj2->Draw ("col")


h_j1Pt->SetStats (0)
h_j1Pt->SetTitle ("")
h_j1Pt->SetLineWidth (2)
h_j1Pt->DrawNormalized ()

h_j1Pt->ComputeIntegral () ;
Double_t * integral = h_j1Pt->GetIntegral () ;
TH1* h_j1Pt_int = (TH1*) h_j1Pt->Clone ("h_j1Pt_int") ;
h_j1Pt_int->SetContent (integral) ;
h_j1Pt_int->SetStats (0)
h_j1Pt_int->SetTitle ("")
h_j1Pt_int->SetLineWidth (2)
.L ReverseCumultive.C 
ReverseCumulative (h_j1Pt_int)

void ReverseCumulative (TH1* histo) 
{
  double integral = histo->GetBinContent (histo->GetNbinsX ()) ;
  for (int iBin = 0; iBin < histo->GetNbinsX (); iBin++) 
    {
     Êdouble value = histo->GetBinContent (iBin+1) ;
     Êhisto->SetBinContent (iBin+1, integral - value) ; 
    }
}

gStyle->SetPalette (1)
.L TH2FCumulative.C 
TH2F * h_ptj1_ptj2_int = TH2FCumulative (h_ptj1_ptj2)
h_ptj1_ptj2_int->Scale (1. / h_ptj1_ptj2->Integral ())
h_ptj1_ptj2_int->SetStats (0)
h_ptj1_ptj2_int->SetTitle (0)
TCanvas c1 ("c1", h_ptj1_ptj2_int->GetName ()) ;
h_ptj1_ptj2_int->Draw ("colz")

TH2F * h_ptj1_ptj2_eff = TH2FCumulative (h_ptj1_ptj2, 1)
h_ptj1_ptj2_eff->Scale (1. / h_ptj1_ptj2->Integral ())
h_ptj1_ptj2_eff->SetStats (0)
h_ptj1_ptj2_eff->SetTitle (0)
TCanvas c2 ("c2", h_ptj1_ptj2_eff->GetName ()) ;
h_ptj1_ptj2_eff->Draw ("colz")

TCanvas c3 ("c3", h_ptj1_ptj2->GetName ())
h_ptj1_ptj2->SetStats (0)
h_ptj1_ptj2->SetTitle ("")
h_ptj1_ptj2->Draw ("colz")


*/
