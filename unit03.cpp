#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <cmath>
#include "hChain.h"

// c++ -o unit03 `root-config --glibs --cflags` -lm hChain.cc unit03.cpp


double 
deltaPhi (double phi1, double phi2)
{

  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  return deltaphi;
}


//PG --------------------------------------------------------   


//! returns  0 : not a valid combination
//! returns  1 : combination compatible w/ a W+
//! returns -1 : combination compatible w/ W-
int testCombination (int inP, int outP, const LHEF::Reader & reader) 
{
  if (
      (abs (reader.hepeup.IDUP.at (inP)) < 5) && (abs (reader.hepeup.IDUP.at (outP)) < 5)   && //PG only u,d,c,s
      (abs (reader.hepeup.IDUP.at (inP)) != abs (reader.hepeup.IDUP.at (outP)))             && //PG not the same quark
      (reader.hepeup.IDUP.at (inP) * reader.hepeup.IDUP.at (outP) > 0)                      && //PG same sign
      ((abs (reader.hepeup.IDUP.at (inP)) > 2) == (abs (reader.hepeup.IDUP.at (outP)) > 2))    //PG same flavour family
    )
    {
      return reader.hepeup.IDUP.at (inP) / abs (reader.hepeup.IDUP.at (inP)) ; 
    }    
  return 0 ;
}


//! ========================================================================================


int main (int argc, char **argv) {

  gROOT->SetStyle ("Plain") ;

  hChain tags_etaProd ("tags_etaProd", "jet tags eta prod", 200, -10, 10, 3) ;
  hChain tags_eta ("tags_eta", "jet tags eta", 200, -7, 7, 3) ;
  hChain tags_Deta ("tags_Deta","jet tags Deta all", 200, 0, 7, 3) ;
  hChain tags_phi ("tags_phi", "tags_phi", 200, -7, 7, 3) ;
  hChain tags_Dphi ("tags_Dphi","jet tags Dphi all", 200, 0, 3.15, 3) ;
  hChain tags_pt ("tags_pt","jet tags pt", 300, 0, 300, 3) ;
  hChain total_pt ("total_pt","total system pt", 300, 0, 300, 3) ;
  hChain Z_pt ("Z_pt","Z pt", 300, 0, 300, 3) ;
  hChain tags_MaxPt ("tags_MaxPt","jet tags Max pt", 300, 0, 300, 3) ;
  hChain tags_MinPt ("tags_MinPt","jet tags Min pt", 300, 0, 300, 3) ;
  hChain Z_shape ("Z_shape","leptonic Z shape", 200, 0, 1000, 3) ;
  hChain quarks_shape ("quarks_shape","tag jets M shape", 200, 0, 1000, 3) ;
  hChain total_shape ("total_shape","total shape", 200, 0, 3000, 3) ;
  hChain tags_zep ("tags_zep","jet tags zeppenfeld var", 200, -7, 7, 3) ;

  hChain leps_pt ("leps_pt","jet leps pt", 300, 0, 300, 3) ;
  hChain leps_MaxPt ("leps_MaxPt","jet leps Max pt", 300, 0, 300, 3) ;
  hChain leps_MinPt ("leps_MinPt","jet leps Min pt", 300, 0, 300, 3) ;
  hChain leps_eta ("leps_eta", "leps_eta", 200, -7, 7, 3) ;
  hChain leps_Deta ("leps_Deta","jet leps Deta all", 200, 0, 7, 3) ;
  hChain leps_phi ("leps_phi", "leps_phi", 200, -7, 7, 3) ;
  hChain leps_Dphi ("leps_Dphi","jet leps Dphi all", 200, 0, 3.15, 3) ;

  TH2F h_0_tags_pt_vs_zep ("h_0_tags_pt_vs_zep","jet tags pt vs zeppenfeld var all",  25,-7,7,25,0,300) ;
  TH2F h_1_tags_pt_vs_zep ("h_1_tags_pt_vs_zep","jet tags pt vs zeppenfeld var sel 1",25,-7,7,25,0,300) ;
  TH2F h_2_tags_pt_vs_zep ("h_2_tags_pt_vs_zep","jet tags pt vs zeppenfeld var sel 2",25,-7,7,25,0,300) ;
  TH2F h_3_tags_pt_vs_zep ("h_3_tags_pt_vs_zep","jet tags pt vs zeppenfeld var sel 3",25,-7,7,25,0,300) ;
  TH2F h_4_tags_pt_vs_zep ("h_4_tags_pt_vs_zep","jet tags pt vs zeppenfeld var sel 4",25,-7,7,25,0,300) ;

  TH2F h_0_tags_maxPt_vs_Deta ("h_0_tags_maxPt_vs_Deta","jet tags max pt vs delta eta var all",  25,0,7,25,0,300) ;
  TH2F h_1_tags_maxPt_vs_Deta ("h_1_tags_maxPt_vs_Deta","jet tags max pt vs delta eta var sel 1",25,0,7,25,0,300) ;
  TH2F h_2_tags_maxPt_vs_Deta ("h_2_tags_maxPt_vs_Deta","jet tags max pt vs delta eta var sel 2",25,0,7,25,0,300) ;
  TH2F h_3_tags_maxPt_vs_Deta ("h_3_tags_maxPt_vs_Deta","jet tags max pt vs delta eta var sel 3",25,0,7,25,0,300) ;
  TH2F h_4_tags_maxPt_vs_Deta ("h_4_tags_maxPt_vs_Deta","jet tags max pt vs delta eta var sel 4",25,0,7,25,0,300) ;

  TH2F h_0_tags_Mtot_vs_Deta ("h_0_tags_Mtot_vs_Deta","final state total M vs delta eta all",  25,0,7,25,0,3000) ;
  TH2F h_1_tags_Mtot_vs_Deta ("h_1_tags_Mtot_vs_Deta","final state total M vs delta eta sel 1",25,0,7,25,0,3000) ;
  TH2F h_2_tags_Mtot_vs_Deta ("h_2_tags_Mtot_vs_Deta","final state total M vs delta eta sel 2",25,0,7,25,0,3000) ;
  TH2F h_3_tags_Mtot_vs_Deta ("h_3_tags_Mtot_vs_Deta","final state total M vs delta eta sel 3",25,0,7,25,0,3000) ;
  TH2F h_4_tags_Mtot_vs_Deta ("h_4_tags_Mtot_vs_Deta","final state total M vs delta eta sel 4",25,0,7,25,0,3000) ;

  TH2F h_0_tags_eta1_vs_eta2 ("h_0_tags_eta1_vs_eta2","jet tags eta 1 vs eta 2 all",  25,-7,7,25,-7,7) ;
  TH2F h_1_tags_eta1_vs_eta2 ("h_1_tags_eta1_vs_eta2","jet tags eta 1 vs eta 2 sel 1",25,-7,7,25,-7,7) ;
  TH2F h_2_tags_eta1_vs_eta2 ("h_2_tags_eta1_vs_eta2","jet tags eta 1 vs eta 2 sel 2",25,-7,7,25,-7,7) ;
  TH2F h_3_tags_eta1_vs_eta2 ("h_3_tags_eta1_vs_eta2","jet tags eta 1 vs eta 2 sel 3",25,-7,7,25,-7,7) ;
  TH2F h_4_tags_eta1_vs_eta2 ("h_4_tags_eta1_vs_eta2","jet tags eta 1 vs eta 2 sel 4",25,-7,7,25,-7,7) ;

  double vHmassWidth = 5. ;
 
  // Open a stream connected to an event file:
  if (argc < 2) exit (1) ;
  std::ifstream ifs(argv[1]);

  // Create the Reader object:
  LHEF::Reader reader(ifs);

  // Print out the header information:
  std::cout << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cout << reader.initComments;

  // Print out the beam energies:
  std::cout << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;

  // Now loop over all events:
  long ieve = 0;
  long VBFnumber = 0;
  long selVBFnumber = 0;
  long selNONVBFnumber = 0;
  while ( reader.readEvent() ) 
    {
      ++ieve;
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
  
      std::vector<int> leptons ;      
      std::vector<int> finalQuarks ;      
      std::vector<int> initialQuarks ;      
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
//          std::cout << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          //PG incoming particle          
          if (reader.hepeup.ISTUP.at (iPart) == -1)
            {
              initialQuarks.push_back (iPart) ;
            } //PG incoming particle          

          //PG outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              //PG leptons
              if (abs (reader.hepeup.IDUP.at (iPart)) == 11 || //PG electron
                  abs (reader.hepeup.IDUP.at (iPart)) == 13)   //PG muon
                {
                  leptons.push_back (iPart) ;
                } //PG leptons
              else
                {
                  finalQuarks.push_back (iPart) ;
                }

            } //PG outgoing particles
        } //PG loop over particles in the event

//PG --------------- S I G N A L   D E F I N I T I O N  -------------------

      if (leptons.size () != 2) std::cerr << "not two leptons in the final state\n" ;
      if (finalQuarks.size () != 2) std::cerr << "not two quarks in the final state\n" ;

      bool isVBF = true ;
            
      //PG controllo che i due quark NON vengano da una Z o W (4)     
      TLorentzVector fquark0
        (
          reader.hepeup.PUP.at (finalQuarks.at (0)).at (0), //PG px
          reader.hepeup.PUP.at (finalQuarks.at (0)).at (1), //PG py
          reader.hepeup.PUP.at (finalQuarks.at (0)).at (2), //PG pz
          reader.hepeup.PUP.at (finalQuarks.at (0)).at (3) //PG E
        ) ;
      
      TLorentzVector fquark1
        (
          reader.hepeup.PUP.at (finalQuarks.at (1)).at (0), //PG px
          reader.hepeup.PUP.at (finalQuarks.at (1)).at (1), //PG py
          reader.hepeup.PUP.at (finalQuarks.at (1)).at (2), //PG pz
          reader.hepeup.PUP.at (finalQuarks.at (1)).at (3) //PG E
        ) ;
      TLorentzVector vCand = fquark0 + fquark1 ;
   
      TLorentzVector fquarkMax, fquarkMin ;
      if (fquark0.Perp () > fquark1.Perp ()) 
        {
          fquarkMax = fquark0 ;
          fquarkMin = fquark1 ;
        }
      else
        { 
          fquarkMax = fquark1 ;
          fquarkMin = fquark0 ;
        }

      TLorentzVector lep0
        (
          reader.hepeup.PUP.at (leptons.at (0)).at (0), //PG px
          reader.hepeup.PUP.at (leptons.at (0)).at (1), //PG py
          reader.hepeup.PUP.at (leptons.at (0)).at (2), //PG pz
          reader.hepeup.PUP.at (leptons.at (0)).at (3) //PG E
        ) ;
      
      TLorentzVector lep1
        (
          reader.hepeup.PUP.at (leptons.at (1)).at (0), //PG px
          reader.hepeup.PUP.at (leptons.at (1)).at (1), //PG py
          reader.hepeup.PUP.at (leptons.at (1)).at (2), //PG pz
          reader.hepeup.PUP.at (leptons.at (1)).at (3) //PG E
        ) ;
      
      TLorentzVector lepMax, lepMin ;
      if (lep0.Perp () > lep1.Perp ()) 
        {
          lepMax = lep0 ;
          lepMin = lep1 ;
        }
      else
        { 
          lepMax = lep1 ;
          lepMin = lep0 ;
        }

      TLorentzVector zCand = lep0 + lep1 ;
      TLorentzVector total = zCand + vCand ;
      
      tags_etaProd.Fill (0, fquark0.Eta () * fquark1.Eta ()) ;
      tags_eta.Fill (0, fquark0.Eta ()) ;
      tags_eta.Fill (0, fquark1.Eta ()) ;
      tags_Deta.Fill (0, fabs (fquark0.Eta () - fquark1.Eta ())) ;
      tags_phi.Fill (0, fquark0.Phi ()) ;
      tags_phi.Fill (0, fquark1.Phi ()) ;
      tags_Dphi.Fill (0, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
      tags_pt.Fill (0, fquark0.Perp ()) ;
      tags_pt.Fill (0, fquark1.Perp ()) ;
      tags_MaxPt.Fill (0, fquarkMax.Perp ()) ;
      tags_MinPt.Fill (0, fquarkMin.Perp ()) ;
      tags_zep.Fill (0, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
      tags_zep.Fill (0, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
      quarks_shape.Fill (0, vCand.M ()) ;
      Z_shape.Fill (0, zCand.M ()) ;
      total_shape.Fill (0, total.M ()) ;
      Z_pt.Fill (0, zCand.Pt ()) ;
      total_pt.Fill (0, total.Pt ()) ;

      leps_eta.Fill (0, lep0.Eta ()) ;
      leps_eta.Fill (0, lep1.Eta ()) ;
      leps_Deta.Fill (0, fabs (lep0.Eta () - lep1.Eta ())) ;
      leps_phi.Fill (0, lep0.Phi ()) ;
      leps_phi.Fill (0, lep1.Phi ()) ;
      leps_Dphi.Fill (0, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
      leps_pt.Fill (0, lep0.Perp ()) ;
      leps_pt.Fill (0, lep1.Perp ()) ;
      leps_MaxPt.Fill (0, lepMax.Perp ()) ;
      leps_MinPt.Fill (0, lepMin.Perp ()) ;
 
      h_0_tags_pt_vs_zep.Fill (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
      h_0_tags_pt_vs_zep.Fill (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
      h_0_tags_maxPt_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
      h_0_tags_Mtot_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
      h_0_tags_eta1_vs_eta2.Fill (fquarkMax.Eta (), fquarkMin.Eta ()) ;
  
      //PG "continue" if NOT VBF
      //PG ---------------------
      
      //PG controllo che i due leptoni abbiano lo stesso sapore (1)
      if (abs (reader.hepeup.IDUP.at (leptons.at (0))) != 
          abs (reader.hepeup.IDUP.at (leptons.at (1)))) isVBF = false ;
          
      if (isVBF)
        {
          h_1_tags_pt_vs_zep.Fill (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          h_1_tags_pt_vs_zep.Fill (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          h_1_tags_maxPt_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          h_1_tags_Mtot_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          h_1_tags_eta1_vs_eta2.Fill (fquarkMax.Eta (), fquarkMin.Eta ()) ;
        }
        
      //PG check the flavour combination (2)
      if ((testCombination (initialQuarks.at (0), finalQuarks.at (1), reader) *
          testCombination (initialQuarks.at (1), finalQuarks.at (0), reader)) != -1 &&
          (testCombination (initialQuarks.at (0), finalQuarks.at (0), reader) *
          testCombination (initialQuarks.at (1), finalQuarks.at (1), reader)) != -1
          ) isVBF = false ;

      if (isVBF)
        {
          h_2_tags_pt_vs_zep.Fill (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          h_2_tags_pt_vs_zep.Fill (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          h_2_tags_maxPt_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          h_2_tags_Mtot_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          h_2_tags_eta1_vs_eta2.Fill (fquarkMax.Eta (), fquarkMin.Eta ()) ;
        }
      
      //PG controllo che i due leptoni vengano da una Z (3)
      if (zCand.M () > (91 + 5) || zCand.M () < (91 - 5)) isVBF = false ;

      if (isVBF)
        {
          h_3_tags_pt_vs_zep.Fill (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          h_3_tags_pt_vs_zep.Fill (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          h_3_tags_maxPt_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          h_3_tags_Mtot_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          h_3_tags_eta1_vs_eta2.Fill (fquarkMax.Eta (), fquarkMin.Eta ()) ;
        }
     
      //PG controllo che i due quark non vengano da una Z o W (4)
      if (vCand.M () < (91 + vHmassWidth) && vCand.M () > (91 - vHmassWidth)) isVBF = false ;
      if (vCand.M () < (80 + vHmassWidth) && vCand.M () > (80 - vHmassWidth)) isVBF = false ;

      if (isVBF)
        {
          tags_etaProd.Fill (1, fquark0.Eta () * fquark1.Eta ()) ;
          tags_eta.Fill (1, fquark0.Eta ()) ;
          tags_eta.Fill (1, fquark1.Eta ()) ;
          tags_Deta.Fill (1, fabs (fquark0.Eta () - fquark1.Eta ())) ;
          tags_phi.Fill (1, fquark0.Phi ()) ;
          tags_phi.Fill (1, fquark1.Phi ()) ;
          tags_Dphi.Fill (1, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          tags_pt.Fill (1, fquark0.Perp ()) ;
          tags_pt.Fill (1, fquark1.Perp ()) ;
          tags_MaxPt.Fill (1, fquarkMax.Perp ()) ;
          tags_MinPt.Fill (1, fquarkMin.Perp ()) ;
          tags_zep.Fill (1, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          tags_zep.Fill (1, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          Z_shape.Fill (1, zCand.M ()) ;
          quarks_shape.Fill (1, vCand.M ()) ;
          total_shape.Fill (1, total.M ()) ;
          Z_pt.Fill (1, zCand.Pt ()) ;
          total_pt.Fill (1, total.Perp ()) ;

          leps_eta.Fill (1, lep0.Eta ()) ;
          leps_eta.Fill (1, lep1.Eta ()) ;
          leps_Deta.Fill (1, fabs (lep0.Eta () - lep1.Eta ())) ;
          leps_phi.Fill (1, lep0.Phi ()) ;
          leps_phi.Fill (1, lep1.Phi ()) ;
          leps_Dphi.Fill (1, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
          leps_pt.Fill (1, lep0.Perp ()) ;
          leps_pt.Fill (1, lep1.Perp ()) ;
          leps_MaxPt.Fill (1, lepMax.Perp ()) ;
          leps_MinPt.Fill (1, lepMin.Perp ()) ;

          h_4_tags_pt_vs_zep.Fill (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          h_4_tags_pt_vs_zep.Fill (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          h_4_tags_maxPt_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          h_4_tags_Mtot_vs_Deta.Fill (fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          h_4_tags_eta1_vs_eta2.Fill (fquarkMax.Eta (), fquarkMin.Eta ()) ;
        }
      
      if (isVBF) { ++VBFnumber ; }

//PG --------------- S A M P L E   S E L E C T I O N  -------------------

      if ( //PG select VBF events
          fabs (fquark0.Eta () - fquark1.Eta ()) > 1.5 &&
          fquark0.Perp () > 30 && // GeV
          fquark1.Perp () > 30 && // GeV
          fquarkMax.Perp () > 50 && // GeV
          fabs (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) > 1 &&
          fabs (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) > 1 &&
          fabs (fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) < 3 &&
          fabs (fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) < 3 &&
          (vCand.M () > (91 + 20) || vCand.M () < (91 - 20)) && //PG low jet energy resol
          (vCand.M () > (80 + 20) || vCand.M () < (80 - 20)) && //PG low jet energy resol
          zCand.M () < (91 + 5) && zCand.M () > (91 - 5) &&
          zCand.Pt () > 30 && // GeV
          deltaPhi (fquark0.Phi (), fquark1.Phi ()) > 1 &&
          lepMin.Perp () > 20 // GeV
         )
        {
          tags_etaProd.Fill (2, fquark0.Eta () * fquark1.Eta ()) ;
          tags_eta.Fill (2, fquark0.Eta ()) ;
          tags_eta.Fill (2, fquark1.Eta ()) ;
          tags_Deta.Fill (2, fabs (fquark0.Eta () - fquark1.Eta ())) ;
          tags_phi.Fill (2, fquark0.Phi ()) ;
          tags_phi.Fill (2, fquark1.Phi ()) ;
          tags_Dphi.Fill (2, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          tags_pt.Fill (2, fquark0.Perp ()) ;
          tags_pt.Fill (2, fquark1.Perp ()) ;
          tags_MaxPt.Fill (2, fquarkMax.Perp ()) ;
          tags_MinPt.Fill (2, fquarkMin.Perp ()) ;
          tags_zep.Fill (2, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          tags_zep.Fill (2, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          Z_shape.Fill (2, zCand.M ()) ;
          quarks_shape.Fill (2, vCand.M ()) ;
          total_shape.Fill (2, total.M ()) ;
          Z_pt.Fill (2, zCand.Pt ()) ;
          total_pt.Fill (2, total.Perp ()) ;

          leps_eta.Fill (2, lep0.Eta ()) ;
          leps_eta.Fill (2, lep1.Eta ()) ;
          leps_Deta.Fill (2, fabs (lep0.Eta () - lep1.Eta ())) ;
          leps_phi.Fill (2, lep0.Phi ()) ;
          leps_phi.Fill (2, lep1.Phi ()) ;
          leps_Dphi.Fill (2, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
          leps_pt.Fill (2, lep0.Perp ()) ;
          leps_pt.Fill (2, lep1.Perp ()) ;
          leps_MaxPt.Fill (2, lepMax.Perp ()) ;
          leps_MinPt.Fill (2, lepMin.Perp ()) ;


          if (isVBF) ++selVBFnumber ;
          else ++selNONVBFnumber ;
        }

    } // Now loop over all events
    
    
  std::cout << "VBF NUMBER " << VBFnumber 
            << "\t fraction : " << VBFnumber / (double) ieve << "\n" ;
  std::cout << "SELECTIONS " << "\n" 
            << "\t VBF : " << selVBFnumber / (double) VBFnumber << "\n" 
            << "\t NONVBF : " << selNONVBFnumber / (double) (ieve - VBFnumber) << "\n" 
            << "\t VBF/NONVBF : " << selVBFnumber / (double) selNONVBFnumber << "\n" ;



//PG --------------- D R A W I N G ---------------------------


  int colors[3] = {92,98,0} ;
  tags_eta.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_Deta.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_pt.SetColors (std::vector<int> (colors, colors + 3)) ;
  total_pt.SetColors (std::vector<int> (colors, colors + 3)) ;
  Z_pt.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_MaxPt.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_MinPt.SetColors (std::vector<int> (colors, colors + 3)) ;
  Z_shape.SetColors (std::vector<int> (colors, colors + 3)) ;
  quarks_shape.SetColors (std::vector<int> (colors, colors + 3)) ;
  total_shape.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_zep.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_etaProd.SetColors (std::vector<int> (colors, colors + 3)) ;

  leps_pt.SetColors (std::vector<int> (colors, colors + 3)) ;
  leps_MaxPt.SetColors (std::vector<int> (colors, colors + 3)) ;
  leps_MinPt.SetColors (std::vector<int> (colors, colors + 3)) ;
  leps_eta.SetColors (std::vector<int> (colors, colors + 3)) ;
  leps_Deta.SetColors (std::vector<int> (colors, colors + 3)) ;
  leps_phi.SetColors (std::vector<int> (colors, colors + 3)) ;
  leps_Dphi.SetColors (std::vector<int> (colors, colors + 3)) ;
  
  tags_phi.SetColors (std::vector<int> (colors, colors + 3)) ;
  tags_Dphi.SetColors (std::vector<int> (colors, colors + 3)) ;

  tags_eta.Print (0,4) ;
  tags_Deta.Print (0,4) ;
  tags_pt.Print (0,4) ;
  total_pt.Print (0,4) ;
  Z_pt.Print (0,4) ;
  tags_MaxPt.Print (0,4) ;
  tags_MinPt.Print (0,4) ;
  Z_shape.Print (1,2) ;
  quarks_shape.Print (1,4) ;
  total_shape.Print (1,4) ;
  tags_zep.Print (0,4) ;

  leps_pt.Print (0,4) ;
  leps_MaxPt.Print (0,4) ;
  leps_MinPt.Print (0,4) ;
  leps_eta.Print (0,4) ;
  leps_Deta.Print (0,4) ;
  leps_phi.Print (0,4) ;
  leps_Dphi.Print (0,4) ;
  tags_phi.Print (0,4) ;
  tags_Dphi.Print (0,4) ;

  tags_etaProd.Print (0,4) ;

  TCanvas c7 ;
  c7.cd () ;
  h_0_tags_pt_vs_zep.SetMarkerColor (9) ;
  h_0_tags_pt_vs_zep.SetLineColor (9) ;
  h_0_tags_pt_vs_zep.SetFillColor (9) ;
  h_0_tags_pt_vs_zep.SetFillStyle (3001) ;
  h_4_tags_pt_vs_zep.SetMarkerColor (98) ;
  h_4_tags_pt_vs_zep.SetLineColor (98) ;
  h_4_tags_pt_vs_zep.SetFillColor (98) ;
  h_4_tags_pt_vs_zep.SetFillStyle (3001) ;
  h_0_tags_pt_vs_zep.Draw ("box") ;
  h_4_tags_pt_vs_zep.Draw ("boxsame") ;
  c7.Print ("h_tagJetsPtVSZep.gif","gif") ;
  h_0_tags_pt_vs_zep.DrawNormalized ("box") ;
  h_4_tags_pt_vs_zep.DrawNormalized ("boxsame") ;
  c7.Print ("hN_tagJetsPtVSZep.gif","gif") ;
   
  TCanvas c10 ;
  c10.cd () ;
  h_0_tags_maxPt_vs_Deta.SetMarkerColor (9) ;
  h_0_tags_maxPt_vs_Deta.SetLineColor (9) ;
  h_0_tags_maxPt_vs_Deta.SetFillColor (9) ;
  h_0_tags_maxPt_vs_Deta.SetFillStyle (3001) ;
  h_4_tags_maxPt_vs_Deta.SetMarkerColor (98) ;
  h_4_tags_maxPt_vs_Deta.SetLineColor (98) ;
  h_4_tags_maxPt_vs_Deta.SetFillColor (98) ;
  h_4_tags_maxPt_vs_Deta.SetFillStyle (3001) ;
  h_0_tags_maxPt_vs_Deta.Draw ("box") ;
  h_4_tags_maxPt_vs_Deta.Draw ("boxsame") ;
  c10.Print ("h_tagJetsMaxPtVSDeta.gif","gif") ;
  h_0_tags_maxPt_vs_Deta.DrawNormalized ("box") ;
  h_4_tags_maxPt_vs_Deta.DrawNormalized ("boxsame") ;
  c10.Print ("hN_tagJetsMaxPtVSDeta.gif","gif") ;
  
  TCanvas c12 ;
  c12.cd () ;
  h_0_tags_Mtot_vs_Deta.SetMarkerColor (9) ;
  h_0_tags_Mtot_vs_Deta.SetLineColor (9) ;
  h_0_tags_Mtot_vs_Deta.SetFillColor (9) ;
  h_0_tags_Mtot_vs_Deta.SetFillStyle (3001) ;
  h_4_tags_Mtot_vs_Deta.SetMarkerColor (98) ;
  h_4_tags_Mtot_vs_Deta.SetLineColor (98) ;
  h_4_tags_Mtot_vs_Deta.SetFillColor (98) ;
  h_4_tags_Mtot_vs_Deta.SetFillStyle (3001) ;
  h_0_tags_Mtot_vs_Deta.Draw ("box") ;
  h_4_tags_Mtot_vs_Deta.Draw ("boxsame") ;
  c12.Print ("h_totMVSDeta.gif","gif") ;
  h_0_tags_Mtot_vs_Deta.DrawNormalized ("box") ;
  h_4_tags_Mtot_vs_Deta.DrawNormalized ("boxsame") ;
  c12.Print ("hN_totMVSDeta.gif","gif") ;
    
  TCanvas c13 ;
  c13.cd () ;
  h_0_tags_eta1_vs_eta2.SetMarkerColor (9) ;
  h_0_tags_eta1_vs_eta2.SetLineColor (9) ;
  h_0_tags_eta1_vs_eta2.SetFillColor (9) ;
  h_0_tags_eta1_vs_eta2.SetFillStyle (3001) ;
  h_4_tags_eta1_vs_eta2.SetMarkerColor (98) ;
  h_4_tags_eta1_vs_eta2.SetLineColor (98) ;
  h_4_tags_eta1_vs_eta2.SetFillColor (98) ;
  h_4_tags_eta1_vs_eta2.SetFillStyle (3001) ;
  h_0_tags_eta1_vs_eta2.Draw ("box") ;
  h_4_tags_eta1_vs_eta2.Draw ("boxsame") ;
  c13.Print ("h_tagJetsPtVSZep.gif","gif") ;
  h_0_tags_eta1_vs_eta2.DrawNormalized ("box") ;
  h_4_tags_eta1_vs_eta2.DrawNormalized ("boxsame") ;
  c13.Print ("hN_tagJetsPtVSZep.gif","gif") ;

  TFile histosFile ("output.root","recreate") ;
  histosFile.cd () ;
  
  tags_eta.Write (histosFile) ;
  tags_Deta.Write (histosFile) ;
  tags_pt.Write (histosFile) ;
  total_pt.Write (histosFile) ;
  Z_pt.Write (histosFile) ;
  tags_MaxPt.Write (histosFile) ;
  tags_MinPt.Write (histosFile) ;
  Z_shape.Write (histosFile) ;
  quarks_shape.Write (histosFile) ;
  total_shape.Write (histosFile) ;
  tags_zep.Write (histosFile) ;
  leps_pt.Write (histosFile) ;
  leps_MaxPt.Write (histosFile) ;
  leps_MinPt.Write (histosFile) ;
  leps_eta.Write (histosFile) ;
  leps_Deta.Write (histosFile) ;
  leps_phi.Write (histosFile) ;
  leps_Dphi.Write (histosFile) ;
  tags_phi.Write (histosFile) ;
  tags_Dphi.Write (histosFile) ;
  tags_etaProd.Write (histosFile) ;
       
  h_0_tags_pt_vs_zep.Write () ;
  h_1_tags_pt_vs_zep.Write () ;
  h_2_tags_pt_vs_zep.Write () ;
  h_3_tags_pt_vs_zep.Write () ;
  h_4_tags_pt_vs_zep.Write () ;

  h_0_tags_maxPt_vs_Deta.Write () ;
  h_1_tags_maxPt_vs_Deta.Write () ;
  h_2_tags_maxPt_vs_Deta.Write () ;
  h_3_tags_maxPt_vs_Deta.Write () ;
  h_4_tags_maxPt_vs_Deta.Write () ;

  h_0_tags_Mtot_vs_Deta.Write () ;
  h_1_tags_Mtot_vs_Deta.Write () ;
  h_2_tags_Mtot_vs_Deta.Write () ;
  h_3_tags_Mtot_vs_Deta.Write () ;
  h_4_tags_Mtot_vs_Deta.Write () ;

  h_0_tags_eta1_vs_eta2.Write () ;
  h_1_tags_eta1_vs_eta2.Write () ;
  h_2_tags_eta1_vs_eta2.Write () ;
  h_3_tags_eta1_vs_eta2.Write () ;
  h_4_tags_eta1_vs_eta2.Write () ;

  histosFile.Close () ;

  // Now we are done.
  return 0 ;
}
