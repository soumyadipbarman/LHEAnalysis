// for the paper with Bho

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

/*

TO COMPILE:

  export ROOTSYS=~/Desktop/root
  export PATH=$ROOTSYS/bin:$PATH
  export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
  export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH 

  c++ -o analysis02 `root-config --glibs --cflags` \
     -lm hFactory.cc hChain.cc h2Factory.cc h2Chain.cc analysis02.cpp
     
TO RUN:     

  ./analysis02 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe 

numerazione dei plot:

  2 -- segnale selezionati
  4 -- segnale persi
  0 -- segnale totale
  3 -- bkg irr selezionati
  5 -- bkg irr persi
  1 -- bkg irr totale
  7 -- bkg qcd selezionati
  8 -- bkg qcd persi
  6 -- bkg qcd totale


*/



double 
deltaPhi (double phi1, double phi2)
{

  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  return deltaphi;
}


//PG --------------------------------------------------------   


//PG selections applied only to EWK samples
bool isRecoSignal (TLorentzVector & fquarkMin, TLorentzVector & fquarkMax, 
                   TLorentzVector & lepMin,    TLorentzVector & lepMax, 
                   TLorentzVector & vCand,     TLorentzVector & zCand)
{
     if ( //PG select VBF events
          fabs (fquarkMin.Eta () - fquarkMax.Eta ()) > 1.5 &&
          fabs (fquarkMin.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ())) > 1 &&
          fabs (fquarkMax.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ())) > 1 &&
          fabs (fquarkMin.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ())) < 3 &&
          fabs (fquarkMax.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ())) < 3 &&
          fquarkMin.Perp () > 30 && // GeV
          fquarkMax.Perp () > 30 && // GeV
          fquarkMax.Perp () > 50 && // GeV
          (vCand.M () > (91 + 20) || vCand.M () < (80 - 20)) && //PG low jet energy resol
          zCand.M () < (91 + 5) && zCand.M () > (91 - 5) &&
          zCand.Pt () > 30 // GeV
         ) return true ;
     return false ;
}



//PG --------------------------------------------------------   


bool isControl (TLorentzVector & fquarkMin, TLorentzVector & fquarkMax, 
                TLorentzVector & lepMin,    TLorentzVector & lepMax, 
                TLorentzVector & vCand,     TLorentzVector & zCand)
{
     if ( //PG select VBF events
//          fabs (fquarkMin.Eta () - fquarkMax.Eta ()) > 1.5 &&
// -------- this is the inverted selection to get the control region
          fquarkMin.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ()) > - 0.6 &&
          fquarkMin.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ()) < 0.6 &&
          fquarkMax.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ()) > - 0.6 &&
          fquarkMax.Eta () - 0.5 * (fquarkMin.Eta () + fquarkMax.Eta ()) < 0.6 &&
// --------
          fquarkMin.Perp () > 30 && // GeV
          fquarkMax.Perp () > 30 && // GeV
          fquarkMax.Perp () > 60 && // GeV
//          (vCand.M () > (91 + 20) || vCand.M () < (91 - 20)) && //PG low jet energy resol
//          (vCand.M () > (80 + 20) || vCand.M () < (80 - 20)) && //PG low jet energy resol
          (vCand.M () > 400) && //PG low jet energy resol
          zCand.M () < (91 + 5) && zCand.M () > (91 - 5) &&
          zCand.Pt () > 30 && // GeV
          deltaPhi (fquarkMin.Phi (), fquarkMax.Phi ()) > 1.5 &&
//          deltaPhi (lepMin.Phi (), lepMax.Phi ()) < 2.8 &&
          lepMin.Perp () > 15 && // GeV
          lepMax.Perp () > 30 && // GeV
          fabs (lepMin.Eta ()) < 2.5 && // detector effect
          fabs (lepMax.Eta ()) < 2.5 // detector effect
         ) return true ;
     return false ;
}


//PG --------------------------------------------------------   
//PG test the combination of the output and input quarks to ID a W+W- production


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

  hFactory H1fact ("output1D_02.root", true) ;

  int numOfHistos = 9 ;
  
  H1fact.add_h1 ("tags_etaProd", "jet tags eta prod", 200, -10, 10, numOfHistos) ;
  H1fact.add_h1 ("tags_eta", "jet tags eta", 200, -7, 7, numOfHistos) ;
  H1fact.add_h1 ("tags_Deta","jet tags Deta", 200, 0, 7, numOfHistos) ;
  H1fact.add_h1 ("tags_phi", "tags_phi", 200, -7, 7, numOfHistos) ;
  H1fact.add_h1 ("tags_Dphi","jet tags Dphi", 200, 0, 3.15, numOfHistos) ;
  H1fact.add_h1 ("tags_pt","jet tags pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("total_pt","total system pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("Z_pt","Z pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("tags_MaxPt","jet tags Max pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("tags_MinPt","jet tags Min pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("Z_shape","leptonic Z shape", 200, 0, 1000, numOfHistos) ;
  H1fact.add_h1 ("quarks_shape","tag jets M shape", 200, 0, 1000, numOfHistos) ;
  H1fact.add_h1 ("total_shape","total shape", 200, 0, 3000, numOfHistos) ;
  H1fact.add_h1 ("tags_zep","jet tags zeppenfeld var", 200, -7, 7, numOfHistos) ;
  H1fact.add_h1 ("leps_pt","jet leps pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("leps_MaxPt","jet leps Max pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("leps_MinPt","jet leps Min pt", 300, 0, 300, numOfHistos) ;
  H1fact.add_h1 ("leps_eta", "leps_eta", 200, -7, 7, numOfHistos) ;
  H1fact.add_h1 ("leps_Deta","leps Deta", 200, 0, 7, numOfHistos) ;
  H1fact.add_h1 ("leps_phi", "leps phi", 200, -7, 7, numOfHistos) ;
  H1fact.add_h1 ("leps_Dphi","jet leps Dphi", 200, 0, 3.15, numOfHistos) ;
  H1fact.add_h1 ("Vs_Deta","vector bosons Deta", 200, 0, 7, numOfHistos) ;

  h2Factory H2fact ("output2D_02.root", true) ;

  H2fact.add_h2 ("tags_pt_vs_zep","jet tags pt vs zeppenfeld var",  25,-7,7,25,0,300, numOfHistos) ;
  H2fact.add_h2 ("tags_maxPt_vs_Deta","jet tags max pt vs delta eta var", 25,0,7,25,0,300, numOfHistos) ;
  H2fact.add_h2 ("tags_Mtot_vs_Deta","final state total M vs delta eta", 25,0,7,25,0,3000, numOfHistos) ;
  H2fact.add_h2 ("tags_eta1_vs_eta2","jet tags eta 1 vs eta 2", 25,-7,7, 25,-7,7, numOfHistos) ;
  H2fact.add_h2 ("tags_z1_vs_Deta","jet tags zeppenfeld high p vs Deta", 25,0,7, 25,-7,7, numOfHistos) ;
  H2fact.add_h2 ("tags_z2_vs_Deta","jet tags zeppenfeld low p vs Deta", 25,0,7, 25,-7,7, numOfHistos) ;
  H2fact.add_h2 ("Peta_vs_Deta","jet tags eta product vs Deta", 25,0,7, 25,-10,10, numOfHistos) ;
  H2fact.add_h2 ("Dphi_vs_Deta","jet tags Dphi vs Deta", 25,0,7, 25,0,3.15, numOfHistos) ;
  H2fact.add_h2 ("DetaL_vs_DetaJ","leptons Deta vs jet tags Deta", 25,0,7, 25,0,7, numOfHistos) ;
  H2fact.add_h2 ("PhiZ_vs_PhiV","phi of lept Z vs phi of hadronic V", 25,-3.15,3.15, 25,-3.15,3.15, numOfHistos) ;
  H2fact.add_h2 ("EtaZ_vs_EtaV","eta of lept Z vs eta of hadronic V", 25,-7,7, 25,-7,7, numOfHistos) ;

  //PG width of the inv mass selections
  double vHmassWidth = 5. ;
 
  // Open a stream connected to an event file:
  if (argc < 3) exit (1) ;
  std::ifstream ifs (argv[1]) ;

  // Create the Reader object:
  LHEF::Reader reader (ifs) ;

  // Print out the header information:
  std::cout << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cout << reader.initComments;

  // Print out the beam energies:
  std::cout << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;

  // Now loop over all events:
  long ieve = 0 ;
  long VBFnumber = 0 ;
  long selVBFnumber = 0 ;
  long selNONVBFnumber = 0 ;
  long ctrl_40_Number = 0 ;
  long ctrl_22_Number = 0 ;
  long ctrl_sig_Number = 0 ;
  
  // loop over (4,0) events
  while ( reader.readEvent () ) 
    {
      ++ieve;
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
  
      std::vector<int> leptons ;      
      std::vector<int> finalQuarks ;      
      std::vector<int> initialQuarks ;      
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

      if (leptons.size () != 2) 
        {
          std::cerr << "not two leptons in the final state\n" ;
          continue ;
        }
      if (finalQuarks.size () != 2) 
        {
          std::cerr << "not two quarks in the final state\n" ;
          continue ;
        }
        
      bool isVBF = true ;
            
      //PG controllo che i due quark NON vengano da una Z o W (4)     

      //PG build the two quarks
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
      //PG the sum pf the two quarks
      TLorentzVector vCand = fquark0 + fquark1 ;
   
      //PG sort the two quarks in pt 
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

      //PG build the two leptons
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
      
      //PG sort the two leptons in pt 
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
 
      //PG check the MC definition of VBF signal
      //PG -------------------------------------
      
      //PG controllo che i due leptoni abbiano lo stesso sapore (1)
      if (abs (reader.hepeup.IDUP.at (leptons.at (0))) != 
          abs (reader.hepeup.IDUP.at (leptons.at (1)))) isVBF = false ;
          
      //PG check the flavour combination (2)
      if ((testCombination (initialQuarks.at (0), finalQuarks.at (1), reader) *
          testCombination (initialQuarks.at (1), finalQuarks.at (0), reader)) != -1 &&
          (testCombination (initialQuarks.at (0), finalQuarks.at (0), reader) *
          testCombination (initialQuarks.at (1), finalQuarks.at (1), reader)) != -1
          ) isVBF = false ;

      //PG controllo che i due leptoni vengano da una Z (3)
      if (zCand.M () > (91 + 5) || zCand.M () < (91 - 5)) isVBF = false ;

      //PG controllo che i due quark non vengano da una Z o W (4)
      if (vCand.M () < (91 + vHmassWidth) && vCand.M () > (91 - vHmassWidth)) isVBF = false ;
      if (vCand.M () < (80 + vHmassWidth) && vCand.M () > (80 - vHmassWidth)) isVBF = false ;

      //PG fill histograms for MC signal
      if (isVBF)
        {
           H1fact.Fill ("tags_etaProd", 0, fquark0.Eta () * fquark1.Eta ()) ;
           H1fact.Fill ("tags_eta", 0, fquark0.Eta ()) ;
           H1fact.Fill ("tags_eta", 0, fquark1.Eta ()) ;
           H1fact.Fill ("tags_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ())) ;
           H1fact.Fill ("tags_phi", 0, fquark0.Phi ()) ;
           H1fact.Fill ("tags_phi", 0, fquark1.Phi ()) ;
           H1fact.Fill ("tags_Dphi", 0, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
           H1fact.Fill ("tags_pt", 0, fquark0.Perp ()) ;
           H1fact.Fill ("tags_pt", 0, fquark1.Perp ()) ;
           H1fact.Fill ("tags_MaxPt", 0, fquarkMax.Perp ()) ;
           H1fact.Fill ("tags_MinPt", 0, fquarkMin.Perp ()) ;
           H1fact.Fill ("tags_zep", 0, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
           H1fact.Fill ("tags_zep", 0, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
           H1fact.Fill ("quarks_shape", 0, vCand.M ()) ;
           H1fact.Fill ("Z_shape", 0, zCand.M ()) ;
           H1fact.Fill ("total_shape", 0, total.M ()) ;
           H1fact.Fill ("Z_pt", 0, zCand.Pt ()) ;
           H1fact.Fill ("total_pt", 0, total.Pt ()) ;
           H1fact.Fill ("leps_eta", 0, lep0.Eta ()) ;
           H1fact.Fill ("leps_eta", 0, lep1.Eta ()) ;
           H1fact.Fill ("leps_Deta", 0, fabs (lep0.Eta () - lep1.Eta ())) ;
           H1fact.Fill ("leps_phi", 0, lep0.Phi ()) ;
           H1fact.Fill ("leps_phi", 0, lep1.Phi ()) ;
           H1fact.Fill ("leps_Dphi", 0, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
           H1fact.Fill ("leps_pt", 0, lep0.Perp ()) ;
           H1fact.Fill ("leps_pt", 0, lep1.Perp ()) ;
           H1fact.Fill ("leps_MaxPt", 0, lepMax.Perp ()) ;
           H1fact.Fill ("leps_MinPt", 0, lepMin.Perp ()) ;
           H1fact.Fill ("Vs_Deta", 0, fabs (vCand.Eta () - zCand.Eta ())) ;
       
           H2fact.Fill ("tags_pt_vs_zep", 0, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
           H2fact.Fill ("tags_pt_vs_zep", 0, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
           H2fact.Fill ("tags_maxPt_vs_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
           H2fact.Fill ("tags_Mtot_vs_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
           H2fact.Fill ("tags_eta1_vs_eta2", 0, fquarkMax.Eta (), fquarkMin.Eta ()) ;
           H2fact.Fill ("tags_z1_vs_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                              fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
           H2fact.Fill ("tags_z2_vs_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                              fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
           H2fact.Fill ("Peta_vs_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
           H2fact.Fill ("Dphi_vs_Deta", 0, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                           deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
           H2fact.Fill ("DetaL_vs_DetaJ", 0, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fabs (lep0.Eta () - lep1.Eta ())) ;
           H2fact.Fill ("PhiZ_vs_PhiV", 0, vCand.Phi (), zCand.Phi ()) ;
           H2fact.Fill ("EtaZ_vs_EtaV", 0, vCand.Eta (), zCand.Eta ()) ;
         } //PG fill histograms for MC signal
      else //PG fill histos for irreducible bkg
        {
          H1fact.Fill ("tags_etaProd", 1, fquark0.Eta () * fquark1.Eta ()) ;
          H1fact.Fill ("tags_eta", 1, fquark0.Eta ()) ;
          H1fact.Fill ("tags_eta", 1, fquark1.Eta ()) ;
          H1fact.Fill ("tags_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ())) ;
          H1fact.Fill ("tags_phi", 1, fquark0.Phi ()) ;
          H1fact.Fill ("tags_phi", 1, fquark1.Phi ()) ;
          H1fact.Fill ("tags_Dphi", 1, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          H1fact.Fill ("tags_pt", 1, fquark0.Perp ()) ;
          H1fact.Fill ("tags_pt", 1, fquark1.Perp ()) ;
          H1fact.Fill ("tags_MaxPt", 1, fquarkMax.Perp ()) ;
          H1fact.Fill ("tags_MinPt", 1, fquarkMin.Perp ()) ;
          H1fact.Fill ("tags_zep", 1, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H1fact.Fill ("tags_zep", 1, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H1fact.Fill ("Z_shape", 1, zCand.M ()) ;
          H1fact.Fill ("quarks_shape", 1, vCand.M ()) ;
          H1fact.Fill ("total_shape", 1, total.M ()) ;
          H1fact.Fill ("Z_pt", 1, zCand.Pt ()) ;
          H1fact.Fill ("total_pt", 1, total.Perp ()) ;
          H1fact.Fill ("leps_eta", 1, lep0.Eta ()) ;
          H1fact.Fill ("leps_eta", 1, lep1.Eta ()) ;
          H1fact.Fill ("leps_Deta", 1, fabs (lep0.Eta () - lep1.Eta ())) ;
          H1fact.Fill ("leps_phi", 1, lep0.Phi ()) ;
          H1fact.Fill ("leps_phi", 1, lep1.Phi ()) ;
          H1fact.Fill ("leps_Dphi", 1, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
          H1fact.Fill ("leps_pt", 1, lep0.Perp ()) ;
          H1fact.Fill ("leps_pt", 1, lep1.Perp ()) ;
          H1fact.Fill ("leps_MaxPt", 1, lepMax.Perp ()) ;
          H1fact.Fill ("leps_MinPt", 1, lepMin.Perp ()) ;
          H1fact.Fill ("Vs_Deta", 1, fabs (vCand.Eta () - zCand.Eta ())) ;

          H2fact.Fill ("tags_pt_vs_zep", 1, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          H2fact.Fill ("tags_pt_vs_zep", 1, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          H2fact.Fill ("tags_maxPt_vs_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          H2fact.Fill ("tags_Mtot_vs_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          H2fact.Fill ("tags_eta1_vs_eta2", 1, fquarkMax.Eta (), fquarkMin.Eta ()) ;
          H2fact.Fill ("tags_z1_vs_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H2fact.Fill ("tags_z2_vs_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H2fact.Fill ("Peta_vs_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
          H2fact.Fill ("Dphi_vs_Deta", 1, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                          deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          H2fact.Fill ("DetaL_vs_DetaJ", 1, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                            fabs (lep0.Eta () - lep1.Eta ())) ;
          H2fact.Fill ("PhiZ_vs_PhiV", 1, vCand.Phi (), zCand.Phi ()) ;
          H2fact.Fill ("EtaZ_vs_EtaV", 1, vCand.Eta (), zCand.Eta ()) ;
        } //PG fill histos for irreducible bkg
      
      if (isVBF) { ++VBFnumber ; }

//PG --------------- S A M P L E   S E L E C T I O N  -------------------

      //PG reco signal selection
      if (isRecoSignal (fquarkMin, fquarkMax, lepMin, lepMax, vCand, zCand))
        {
          //PG fill histograms for MC signal
          if (isVBF)
            {
              H1fact.Fill ("tags_etaProd", 2, fquark0.Eta () * fquark1.Eta ()) ;
              H1fact.Fill ("tags_eta", 2, fquark0.Eta ()) ;
              H1fact.Fill ("tags_eta", 2, fquark1.Eta ()) ;
              H1fact.Fill ("tags_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ())) ;
              H1fact.Fill ("tags_phi", 2, fquark0.Phi ()) ;
              H1fact.Fill ("tags_phi", 2, fquark1.Phi ()) ;
              H1fact.Fill ("tags_Dphi", 2, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H1fact.Fill ("tags_pt", 2, fquark0.Perp ()) ;
              H1fact.Fill ("tags_pt", 2, fquark1.Perp ()) ;
              H1fact.Fill ("tags_MaxPt", 2, fquarkMax.Perp ()) ;
              H1fact.Fill ("tags_MinPt", 2, fquarkMin.Perp ()) ;
              H1fact.Fill ("tags_zep", 2, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("tags_zep", 2, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("Z_shape", 2, zCand.M ()) ;
              H1fact.Fill ("quarks_shape", 2, vCand.M ()) ;
              H1fact.Fill ("total_shape", 2, total.M ()) ;
              H1fact.Fill ("Z_pt", 2, zCand.Pt ()) ;
              H1fact.Fill ("total_pt", 2, total.Perp ()) ;
              H1fact.Fill ("leps_eta", 2, lep0.Eta ()) ;
              H1fact.Fill ("leps_eta", 2, lep1.Eta ()) ;
              H1fact.Fill ("leps_Deta", 2, fabs (lep0.Eta () - lep1.Eta ())) ;
              H1fact.Fill ("leps_phi", 2, lep0.Phi ()) ;
              H1fact.Fill ("leps_phi", 2, lep1.Phi ()) ;
              H1fact.Fill ("leps_Dphi", 2, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
              H1fact.Fill ("leps_pt", 2, lep0.Perp ()) ;
              H1fact.Fill ("leps_pt", 2, lep1.Perp ()) ;
              H1fact.Fill ("leps_MaxPt", 2, lepMax.Perp ()) ;
              H1fact.Fill ("leps_MinPt", 2, lepMin.Perp ()) ;
              H1fact.Fill ("Vs_Deta", 2, fabs (vCand.Eta () - zCand.Eta ())) ;
    
              H2fact.Fill ("tags_pt_vs_zep", 2, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
              H2fact.Fill ("tags_pt_vs_zep", 2, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
              H2fact.Fill ("tags_maxPt_vs_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
              H2fact.Fill ("tags_Mtot_vs_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
              H2fact.Fill ("tags_eta1_vs_eta2", 2, fquarkMax.Eta (), fquarkMin.Eta ()) ;
              H2fact.Fill ("tags_z1_vs_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("tags_z2_vs_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("Peta_vs_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
              H2fact.Fill ("Dphi_vs_Deta", 2, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                              deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H2fact.Fill ("DetaL_vs_DetaJ", 2, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                fabs (lep0.Eta () - lep1.Eta ())) ;
              H2fact.Fill ("PhiZ_vs_PhiV", 2, vCand.Phi (), zCand.Phi ()) ;
              H2fact.Fill ("EtaZ_vs_EtaV", 2, vCand.Eta (), zCand.Eta ()) ;
            }
          else //PG fill histos for irreducible bkg
            {
              H1fact.Fill ("tags_etaProd", 3, fquark0.Eta () * fquark1.Eta ()) ;
              H1fact.Fill ("tags_eta", 3, fquark0.Eta ()) ;
              H1fact.Fill ("tags_eta", 3, fquark1.Eta ()) ;
              H1fact.Fill ("tags_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ())) ;
              H1fact.Fill ("tags_phi", 3, fquark0.Phi ()) ;
              H1fact.Fill ("tags_phi", 3, fquark1.Phi ()) ;
              H1fact.Fill ("tags_Dphi", 3, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H1fact.Fill ("tags_pt", 3, fquark0.Perp ()) ;
              H1fact.Fill ("tags_pt", 3, fquark1.Perp ()) ;
              H1fact.Fill ("tags_MaxPt", 3, fquarkMax.Perp ()) ;
              H1fact.Fill ("tags_MinPt", 3, fquarkMin.Perp ()) ;
              H1fact.Fill ("tags_zep", 3, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("tags_zep", 3, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("Z_shape", 3, zCand.M ()) ;
              H1fact.Fill ("quarks_shape", 3, vCand.M ()) ;
              H1fact.Fill ("total_shape", 3, total.M ()) ;
              H1fact.Fill ("Z_pt", 3, zCand.Pt ()) ;
              H1fact.Fill ("total_pt", 3, total.Perp ()) ;
              H1fact.Fill ("leps_eta", 3, lep0.Eta ()) ;
              H1fact.Fill ("leps_eta", 3, lep1.Eta ()) ;
              H1fact.Fill ("leps_Deta", 3, fabs (lep0.Eta () - lep1.Eta ())) ;
              H1fact.Fill ("leps_phi", 3, lep0.Phi ()) ;
              H1fact.Fill ("leps_phi", 3, lep1.Phi ()) ;
              H1fact.Fill ("leps_Dphi", 3, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
              H1fact.Fill ("leps_pt", 3, lep0.Perp ()) ;
              H1fact.Fill ("leps_pt", 3, lep1.Perp ()) ;
              H1fact.Fill ("leps_MaxPt", 3, lepMax.Perp ()) ;
              H1fact.Fill ("leps_MinPt", 3, lepMin.Perp ()) ;
              H1fact.Fill ("Vs_Deta", 3, fabs (vCand.Eta () - zCand.Eta ())) ;
    
              H2fact.Fill ("tags_pt_vs_zep", 3, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
              H2fact.Fill ("tags_pt_vs_zep", 3, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
              H2fact.Fill ("tags_maxPt_vs_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
              H2fact.Fill ("tags_Mtot_vs_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
              H2fact.Fill ("tags_eta1_vs_eta2", 3, fquarkMax.Eta (), fquarkMin.Eta ()) ;
              H2fact.Fill ("tags_z1_vs_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("tags_z2_vs_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("Peta_vs_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
              H2fact.Fill ("Dphi_vs_Deta", 3, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                              deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H2fact.Fill ("DetaL_vs_DetaJ", 3, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                fabs (lep0.Eta () - lep1.Eta ())) ;
              H2fact.Fill ("PhiZ_vs_PhiV", 3, vCand.Phi (), zCand.Phi ()) ;
              H2fact.Fill ("EtaZ_vs_EtaV", 3, vCand.Eta (), zCand.Eta ()) ;
            }    
          if (isVBF) ++selVBFnumber ;
          else ++selNONVBFnumber ;
        } //PG reco signal selection
      else
        { //PG reco bkg selection
          if (isVBF)
            {
              H1fact.Fill ("tags_etaProd", 4, fquark0.Eta () * fquark1.Eta ()) ;
              H1fact.Fill ("tags_eta", 4, fquark0.Eta ()) ;
              H1fact.Fill ("tags_eta", 4, fquark1.Eta ()) ;
              H1fact.Fill ("tags_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ())) ;
              H1fact.Fill ("tags_phi", 4, fquark0.Phi ()) ;
              H1fact.Fill ("tags_phi", 4, fquark1.Phi ()) ;
              H1fact.Fill ("tags_Dphi", 4, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H1fact.Fill ("tags_pt", 4, fquark0.Perp ()) ;
              H1fact.Fill ("tags_pt", 4, fquark1.Perp ()) ;
              H1fact.Fill ("tags_MaxPt", 4, fquarkMax.Perp ()) ;
              H1fact.Fill ("tags_MinPt", 4, fquarkMin.Perp ()) ;
              H1fact.Fill ("tags_zep", 4, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("tags_zep", 4, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("Z_shape", 4, zCand.M ()) ;
              H1fact.Fill ("quarks_shape", 4, vCand.M ()) ;
              H1fact.Fill ("total_shape", 4, total.M ()) ;
              H1fact.Fill ("Z_pt", 4, zCand.Pt ()) ;
              H1fact.Fill ("total_pt", 4, total.Perp ()) ;
              H1fact.Fill ("leps_eta", 4, lep0.Eta ()) ;
              H1fact.Fill ("leps_eta", 4, lep1.Eta ()) ;
              H1fact.Fill ("leps_Deta", 4, fabs (lep0.Eta () - lep1.Eta ())) ;
              H1fact.Fill ("leps_phi", 4, lep0.Phi ()) ;
              H1fact.Fill ("leps_phi", 4, lep1.Phi ()) ;
              H1fact.Fill ("leps_Dphi", 4, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
              H1fact.Fill ("leps_pt", 4, lep0.Perp ()) ;
              H1fact.Fill ("leps_pt", 4, lep1.Perp ()) ;
              H1fact.Fill ("leps_MaxPt", 4, lepMax.Perp ()) ;
              H1fact.Fill ("leps_MinPt", 4, lepMin.Perp ()) ;
              H1fact.Fill ("Vs_Deta", 4, fabs (vCand.Eta () - zCand.Eta ())) ;
    
              H2fact.Fill ("tags_pt_vs_zep", 4, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
              H2fact.Fill ("tags_pt_vs_zep", 4, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
              H2fact.Fill ("tags_maxPt_vs_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
              H2fact.Fill ("tags_Mtot_vs_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
              H2fact.Fill ("tags_eta1_vs_eta2", 4, fquarkMax.Eta (), fquarkMin.Eta ()) ;
              H2fact.Fill ("tags_z1_vs_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("tags_z2_vs_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("Peta_vs_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
              H2fact.Fill ("Dphi_vs_Deta", 4, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                              deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H2fact.Fill ("DetaL_vs_DetaJ", 4, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                fabs (lep0.Eta () - lep1.Eta ())) ;
              H2fact.Fill ("PhiZ_vs_PhiV", 4, vCand.Phi (), zCand.Phi ()) ;
              H2fact.Fill ("EtaZ_vs_EtaV", 4, vCand.Eta (), zCand.Eta ()) ;
            }
          else //PG fill histos for irreducible bkg
            {
              H1fact.Fill ("tags_etaProd", 5, fquark0.Eta () * fquark1.Eta ()) ;
              H1fact.Fill ("tags_eta", 5, fquark0.Eta ()) ;
              H1fact.Fill ("tags_eta", 5, fquark1.Eta ()) ;
              H1fact.Fill ("tags_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ())) ;
              H1fact.Fill ("tags_phi", 5, fquark0.Phi ()) ;
              H1fact.Fill ("tags_phi", 5, fquark1.Phi ()) ;
              H1fact.Fill ("tags_Dphi", 5, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H1fact.Fill ("tags_pt", 5, fquark0.Perp ()) ;
              H1fact.Fill ("tags_pt", 5, fquark1.Perp ()) ;
              H1fact.Fill ("tags_MaxPt", 5, fquarkMax.Perp ()) ;
              H1fact.Fill ("tags_MinPt", 5, fquarkMin.Perp ()) ;
              H1fact.Fill ("tags_zep", 5, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("tags_zep", 5, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H1fact.Fill ("Z_shape", 5, zCand.M ()) ;
              H1fact.Fill ("quarks_shape", 5, vCand.M ()) ;
              H1fact.Fill ("total_shape", 5, total.M ()) ;
              H1fact.Fill ("Z_pt", 5, zCand.Pt ()) ;
              H1fact.Fill ("total_pt", 5, total.Perp ()) ;
              H1fact.Fill ("leps_eta", 5, lep0.Eta ()) ;
              H1fact.Fill ("leps_eta", 5, lep1.Eta ()) ;
              H1fact.Fill ("leps_Deta", 5, fabs (lep0.Eta () - lep1.Eta ())) ;
              H1fact.Fill ("leps_phi", 5, lep0.Phi ()) ;
              H1fact.Fill ("leps_phi", 5, lep1.Phi ()) ;
              H1fact.Fill ("leps_Dphi", 5, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
              H1fact.Fill ("leps_pt", 5, lep0.Perp ()) ;
              H1fact.Fill ("leps_pt", 5, lep1.Perp ()) ;
              H1fact.Fill ("leps_MaxPt", 5, lepMax.Perp ()) ;
              H1fact.Fill ("leps_MinPt", 5, lepMin.Perp ()) ;
              H1fact.Fill ("Vs_Deta", 5, fabs (vCand.Eta () - zCand.Eta ())) ;
    
              H2fact.Fill ("tags_pt_vs_zep", 5, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
              H2fact.Fill ("tags_pt_vs_zep", 5, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
              H2fact.Fill ("tags_maxPt_vs_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
              H2fact.Fill ("tags_Mtot_vs_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
              H2fact.Fill ("tags_eta1_vs_eta2", 5, fquarkMax.Eta (), fquarkMin.Eta ()) ;
              H2fact.Fill ("tags_z1_vs_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("tags_z2_vs_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                 fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
              H2fact.Fill ("Peta_vs_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
              H2fact.Fill ("Dphi_vs_Deta", 5, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                              deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
              H2fact.Fill ("DetaL_vs_DetaJ", 5, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                                fabs (lep0.Eta () - lep1.Eta ())) ;
              H2fact.Fill ("PhiZ_vs_PhiV", 5, vCand.Phi (), zCand.Phi ()) ;
              H2fact.Fill ("EtaZ_vs_EtaV", 5, vCand.Eta (), zCand.Eta ()) ;
            }    
        
        } //PG reco bkg selection
      
      if (isControl (fquarkMin, fquarkMax, lepMin, lepMax, vCand, zCand))
        {
          ++ctrl_40_Number ;
          if (isVBF) ++ctrl_sig_Number ;
        }

    } // loop over (4,0) events
    
  //PG loop over bkg
  //PG -------------

  int BKGnumber = 0 ;
  int BKGnumberWithObj = 0 ;
  int VBFBKGnumber = 0 ;

  std::ifstream ifsbkg (argv[2]) ;
  // Create the Reader object
  LHEF::Reader bkgReader (ifsbkg) ;

  //PG loop over BKG
  while ( bkgReader.readEvent () ) 
    {
      ++BKGnumber;
      if (BKGnumber % 1000 == 0) std::cout << "BKG event " << BKGnumber << "\n" ;
  
      std::vector<int> leptons ;      
      std::vector<int> finalQuarks ;      
      std::vector<int> initialQuarks ;      
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < bkgReader.hepeup.IDUP.size (); ++iPart)
        {
//          std::cout << "\t part type [" << iPart << "] " << bkgReader.hepeup.IDUP.at (iPart)
//                    << "\t status " << bkgReader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          //PG incoming particle          
          if (bkgReader.hepeup.ISTUP.at (iPart) == -1)
            {
              initialQuarks.push_back (iPart) ;
            } //PG incoming particle          

          //PG outgoing particles          
          if (bkgReader.hepeup.ISTUP.at (iPart) == 1)
            {
              //PG leptons
              if (abs (bkgReader.hepeup.IDUP.at (iPart)) == 11 || //PG electron
                  abs (bkgReader.hepeup.IDUP.at (iPart)) == 13)   //PG muon
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

      if (leptons.size () != 2) 
        {
//          std::cerr << "bkg not two leptons in the final state\n" ;
          continue ;
        }
      if (finalQuarks.size () != 2) 
        {
//          std::cerr << "bkg not two quarks in the final state\n" ;
          continue ;
        }
        
      ++BKGnumberWithObj ;
      
      //PG controllo che i due quark NON vengano da una Z o W (4)     
      TLorentzVector fquark0
        (
          bkgReader.hepeup.PUP.at (finalQuarks.at (0)).at (0), //PG px
          bkgReader.hepeup.PUP.at (finalQuarks.at (0)).at (1), //PG py
          bkgReader.hepeup.PUP.at (finalQuarks.at (0)).at (2), //PG pz
          bkgReader.hepeup.PUP.at (finalQuarks.at (0)).at (3) //PG E
        ) ;
      
      TLorentzVector fquark1
        (
          bkgReader.hepeup.PUP.at (finalQuarks.at (1)).at (0), //PG px
          bkgReader.hepeup.PUP.at (finalQuarks.at (1)).at (1), //PG py
          bkgReader.hepeup.PUP.at (finalQuarks.at (1)).at (2), //PG pz
          bkgReader.hepeup.PUP.at (finalQuarks.at (1)).at (3) //PG E
        ) ;
      TLorentzVector vCand = fquark0 + fquark1 ;
   
      TLorentzVector fquarkMax = fquark0 ;
      TLorentzVector fquarkMin = fquark1 ;
      if (fquarkMax.Perp () < fquarkMin.Perp ()) std::swap (fquarkMax, fquarkMin) ;

      TLorentzVector lep0
        (
          bkgReader.hepeup.PUP.at (leptons.at (0)).at (0), //PG px
          bkgReader.hepeup.PUP.at (leptons.at (0)).at (1), //PG py
          bkgReader.hepeup.PUP.at (leptons.at (0)).at (2), //PG pz
          bkgReader.hepeup.PUP.at (leptons.at (0)).at (3) //PG E
        ) ;
      
      TLorentzVector lep1
        (
          bkgReader.hepeup.PUP.at (leptons.at (1)).at (0), //PG px
          bkgReader.hepeup.PUP.at (leptons.at (1)).at (1), //PG py
          bkgReader.hepeup.PUP.at (leptons.at (1)).at (2), //PG pz
          bkgReader.hepeup.PUP.at (leptons.at (1)).at (3) //PG E
        ) ;
      
      TLorentzVector lepMax = lep0 ;
      TLorentzVector lepMin = lep1 ;
      if (lepMax.Perp () < lepMin.Perp ()) std::swap (lepMax, lepMin) ;

      TLorentzVector zCand = lep0 + lep1 ;
      TLorentzVector total = zCand + vCand ;
      
      H1fact.Fill ("tags_etaProd", 6, fquark0.Eta () * fquark1.Eta ()) ;
      H1fact.Fill ("tags_eta", 6, fquark0.Eta ()) ;
      H1fact.Fill ("tags_eta", 6, fquark1.Eta ()) ;
      H1fact.Fill ("tags_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ())) ;
      H1fact.Fill ("tags_phi", 6, fquark0.Phi ()) ;
      H1fact.Fill ("tags_phi", 6, fquark1.Phi ()) ;
      H1fact.Fill ("tags_Dphi", 6, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
      H1fact.Fill ("tags_pt", 6, fquark0.Perp ()) ;
      H1fact.Fill ("tags_pt", 6, fquark1.Perp ()) ;
      H1fact.Fill ("tags_MaxPt", 6, fquarkMax.Perp ()) ;
      H1fact.Fill ("tags_MinPt", 6, fquarkMin.Perp ()) ;
      H1fact.Fill ("tags_zep", 6, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
      H1fact.Fill ("tags_zep", 6, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
      H1fact.Fill ("Z_shape", 6, zCand.M ()) ;
      H1fact.Fill ("quarks_shape", 6, vCand.M ()) ;
      H1fact.Fill ("total_shape", 6, total.M ()) ;
      H1fact.Fill ("Z_pt", 6, zCand.Pt ()) ;
      H1fact.Fill ("total_pt", 6, total.Perp ()) ;
      H1fact.Fill ("leps_eta", 6, lep0.Eta ()) ;
      H1fact.Fill ("leps_eta", 6, lep1.Eta ()) ;
      H1fact.Fill ("leps_Deta", 6, fabs (lep0.Eta () - lep1.Eta ())) ;
      H1fact.Fill ("leps_phi", 6, lep0.Phi ()) ;
      H1fact.Fill ("leps_phi", 6, lep1.Phi ()) ;
      H1fact.Fill ("leps_Dphi", 6, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
      H1fact.Fill ("leps_pt", 6, lep0.Perp ()) ;
      H1fact.Fill ("leps_pt", 6, lep1.Perp ()) ;
      H1fact.Fill ("leps_MaxPt", 6, lepMax.Perp ()) ;
      H1fact.Fill ("leps_MinPt", 6, lepMin.Perp ()) ;
      H1fact.Fill ("Vs_Deta", 6, fabs (vCand.Eta () - zCand.Eta ())) ;

      H2fact.Fill ("tags_pt_vs_zep", 6, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
      H2fact.Fill ("tags_pt_vs_zep", 6, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
      H2fact.Fill ("tags_maxPt_vs_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
      H2fact.Fill ("tags_Mtot_vs_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
      H2fact.Fill ("tags_eta1_vs_eta2", 6, fquarkMax.Eta (), fquarkMin.Eta ()) ;
      H2fact.Fill ("tags_z1_vs_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                         fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
      H2fact.Fill ("tags_z2_vs_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                         fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
      H2fact.Fill ("Peta_vs_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
      H2fact.Fill ("Dphi_vs_Deta", 6, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                      deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
      H2fact.Fill ("DetaL_vs_DetaJ", 6, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                        fabs (lep0.Eta () - lep1.Eta ())) ;
      H2fact.Fill ("PhiZ_vs_PhiV", 6, vCand.Phi (), zCand.Phi ()) ;
      H2fact.Fill ("EtaZ_vs_EtaV", 6, vCand.Eta (), zCand.Eta ()) ;

      //PG reco signal selection
      if (isRecoSignal (fquarkMin, fquarkMax, lepMin, lepMax, vCand, zCand))
        {
          H1fact.Fill ("tags_etaProd", 7, fquark0.Eta () * fquark1.Eta ()) ;
          H1fact.Fill ("tags_eta", 7, fquark0.Eta ()) ;
          H1fact.Fill ("tags_eta", 7, fquark1.Eta ()) ;
          H1fact.Fill ("tags_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ())) ;
          H1fact.Fill ("tags_phi", 7, fquark0.Phi ()) ;
          H1fact.Fill ("tags_phi", 7, fquark1.Phi ()) ;
          H1fact.Fill ("tags_Dphi", 7, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          H1fact.Fill ("tags_pt", 7, fquark0.Perp ()) ;
          H1fact.Fill ("tags_pt", 7, fquark1.Perp ()) ;
          H1fact.Fill ("tags_MaxPt", 7, fquarkMax.Perp ()) ;
          H1fact.Fill ("tags_MinPt", 7, fquarkMin.Perp ()) ;
          H1fact.Fill ("tags_zep", 7, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H1fact.Fill ("tags_zep", 7, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H1fact.Fill ("Z_shape", 7, zCand.M ()) ;
          H1fact.Fill ("quarks_shape", 7, vCand.M ()) ;
          H1fact.Fill ("total_shape", 7, total.M ()) ;
          H1fact.Fill ("Z_pt", 7, zCand.Pt ()) ;
          H1fact.Fill ("total_pt", 7, total.Perp ()) ;
          H1fact.Fill ("leps_eta", 7, lep0.Eta ()) ;
          H1fact.Fill ("leps_eta", 7, lep1.Eta ()) ;
          H1fact.Fill ("leps_Deta", 7, fabs (lep0.Eta () - lep1.Eta ())) ;
          H1fact.Fill ("leps_phi", 7, lep0.Phi ()) ;
          H1fact.Fill ("leps_phi", 7, lep1.Phi ()) ;
          H1fact.Fill ("leps_Dphi", 7, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
          H1fact.Fill ("leps_pt", 7, lep0.Perp ()) ;
          H1fact.Fill ("leps_pt", 7, lep1.Perp ()) ;
          H1fact.Fill ("leps_MaxPt", 7, lepMax.Perp ()) ;
          H1fact.Fill ("leps_MinPt", 7, lepMin.Perp ()) ;
          H1fact.Fill ("Vs_Deta", 7, fabs (vCand.Eta () - zCand.Eta ())) ;

          H2fact.Fill ("tags_pt_vs_zep", 7, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          H2fact.Fill ("tags_pt_vs_zep", 7, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          H2fact.Fill ("tags_maxPt_vs_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          H2fact.Fill ("tags_Mtot_vs_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          H2fact.Fill ("tags_eta1_vs_eta2", 7, fquarkMax.Eta (), fquarkMin.Eta ()) ;
          H2fact.Fill ("tags_z1_vs_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H2fact.Fill ("tags_z2_vs_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H2fact.Fill ("Peta_vs_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
          H2fact.Fill ("Dphi_vs_Deta", 7, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                          deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          H2fact.Fill ("DetaL_vs_DetaJ", 7, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                            fabs (lep0.Eta () - lep1.Eta ())) ;
          H2fact.Fill ("PhiZ_vs_PhiV", 7, vCand.Phi (), zCand.Phi ()) ;
          H2fact.Fill ("EtaZ_vs_EtaV", 7, vCand.Eta (), zCand.Eta ()) ;

          ++VBFBKGnumber ;
        } //PG reco signal selection
      else
        { //PG reco bkg selection
          H1fact.Fill ("tags_etaProd", 8, fquark0.Eta () * fquark1.Eta ()) ;
          H1fact.Fill ("tags_eta", 8, fquark0.Eta ()) ;
          H1fact.Fill ("tags_eta", 8, fquark1.Eta ()) ;
          H1fact.Fill ("tags_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ())) ;
          H1fact.Fill ("tags_phi", 8, fquark0.Phi ()) ;
          H1fact.Fill ("tags_phi", 8, fquark1.Phi ()) ;
          H1fact.Fill ("tags_Dphi", 8, deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          H1fact.Fill ("tags_pt", 8, fquark0.Perp ()) ;
          H1fact.Fill ("tags_pt", 8, fquark1.Perp ()) ;
          H1fact.Fill ("tags_MaxPt", 8, fquarkMax.Perp ()) ;
          H1fact.Fill ("tags_MinPt", 8, fquarkMin.Perp ()) ;
          H1fact.Fill ("tags_zep", 8, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H1fact.Fill ("tags_zep", 8, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H1fact.Fill ("Z_shape", 8, zCand.M ()) ;
          H1fact.Fill ("quarks_shape", 8, vCand.M ()) ;
          H1fact.Fill ("total_shape", 8, total.M ()) ;
          H1fact.Fill ("Z_pt", 8, zCand.Pt ()) ;
          H1fact.Fill ("total_pt", 8, total.Perp ()) ;
          H1fact.Fill ("leps_eta", 8, lep0.Eta ()) ;
          H1fact.Fill ("leps_eta", 8, lep1.Eta ()) ;
          H1fact.Fill ("leps_Deta", 8, fabs (lep0.Eta () - lep1.Eta ())) ;
          H1fact.Fill ("leps_phi", 8, lep0.Phi ()) ;
          H1fact.Fill ("leps_phi", 8, lep1.Phi ()) ;
          H1fact.Fill ("leps_Dphi", 8, deltaPhi (lep0.Phi (), lep1.Phi ())) ;
          H1fact.Fill ("leps_pt", 8, lep0.Perp ()) ;
          H1fact.Fill ("leps_pt", 8, lep1.Perp ()) ;
          H1fact.Fill ("leps_MaxPt", 8, lepMax.Perp ()) ;
          H1fact.Fill ("leps_MinPt", 8, lepMin.Perp ()) ;
          H1fact.Fill ("Vs_Deta", 8, fabs (vCand.Eta () - zCand.Eta ())) ;

          H2fact.Fill ("tags_pt_vs_zep", 8, fquark0.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark0.Perp ()) ;
          H2fact.Fill ("tags_pt_vs_zep", 8, fquark1.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ()), fquark1.Perp ()) ;
          H2fact.Fill ("tags_maxPt_vs_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ()), fquarkMax.Perp ()) ;
          H2fact.Fill ("tags_Mtot_vs_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ()), total.M ()) ;
          H2fact.Fill ("tags_eta1_vs_eta2", 8, fquarkMax.Eta (), fquarkMin.Eta ()) ;
          H2fact.Fill ("tags_z1_vs_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fquarkMax.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H2fact.Fill ("tags_z2_vs_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                             fquarkMin.Eta () - 0.5 * (fquark0.Eta () + fquark1.Eta ())) ;
          H2fact.Fill ("Peta_vs_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ()), fquark0.Eta () * fquark1.Eta ()) ;
          H2fact.Fill ("Dphi_vs_Deta", 8, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                          deltaPhi (fquark0.Phi (), fquark1.Phi ())) ;
          H2fact.Fill ("DetaL_vs_DetaJ", 8, fabs (fquark0.Eta () - fquark1.Eta ()), 
                                            fabs (lep0.Eta () - lep1.Eta ())) ;
          H2fact.Fill ("PhiZ_vs_PhiV", 8, vCand.Phi (), zCand.Phi ()) ;
          H2fact.Fill ("EtaZ_vs_EtaV", 8, vCand.Eta (), zCand.Eta ()) ;
        	
        } //PG reco bkg selection
      if (isControl (fquarkMin, fquarkMax, lepMin, lepMax, vCand, zCand))
        {
          ++ctrl_22_Number ;
        }
        
    } //PG loop over BKG
    
  double norm40 = 4.026 / 32500. ;
//  double norm22 = 2679.48 / 100000. ;
  double norm22 = 2079.47 / 100000. ; //PG pietro_test15
    
  std::cout << "ALL NUMBERS IN pb\n"
            << "SCALE FACTORS:\n"
            << "\t norm40 : " << norm40 << "\n"
            << "\t norm22 : " << norm22 << "\n"
            << "SIGNAL X-SECTION     " << VBFnumber * norm40 << " pb\n"
            << "EWK SIGNAL PURITY    " << selVBFnumber / static_cast<double>(selVBFnumber + selNONVBFnumber) << "\n"
            << "EWK BKG X-SECTION    " << (ieve -VBFnumber) * norm40 << " pb\n"
            << "SELECTED SIGNAL      " << selVBFnumber * norm40 << " pb"
            << "\t (eff : " << selVBFnumber * norm40 / (VBFnumber * norm40) << ")\n" 
            << "SELECTED EWK BKG     " << selNONVBFnumber * norm40 << " pb"
            << "\t (eff : " << selNONVBFnumber * norm40 / ((ieve -VBFnumber) * norm40) << ")\n" 
            << "SELECTED QCD BKG     " << VBFBKGnumber * norm22 << " pb"
            << "\t (eff : " << VBFBKGnumber * norm22 / (BKGnumber * norm22) << ")\n" 
            << "SELECTED TOT BKG     " << selNONVBFnumber * norm40 + VBFBKGnumber * norm22 << " pb"
            << "\t (eff : " << (selNONVBFnumber * norm40 + VBFBKGnumber * norm22) / 
                               (double) ((ieve - VBFnumber) * norm40 + BKGnumber * norm22) << ")\n"
            << "SIGNAL OVER BKG      " << selVBFnumber * norm40 / (selNONVBFnumber * norm40 + VBFBKGnumber * norm22) << "\n"

            << "TOTAL SAMPLE         " << (ieve * norm40 + BKGnumber * norm22) << " pb\n" 
            << "TOTAL SEL. SAMPLE    " << (selVBFnumber * norm40  + selNONVBFnumber * norm40 + VBFBKGnumber * norm22) 
                                       << " pb\n" 
            << "\t ctrl ctrl_40_Number  : " << ctrl_40_Number  * norm40 << " pb\n"
            << "\t ctrl ctrl_22_Number  : " << ctrl_22_Number  * norm22 << " pb\n"
            << "\t ctrl ctrl_sig_Number : " << ctrl_sig_Number * norm40 << " pb\n" ;

std::cout << "CHECK NUMBERS : "
          << "\tBKGnumber = " << BKGnumber
          << "\tBKGnumberWithObj = " << BKGnumberWithObj
          << "\tVBFBKGnumber = " << VBFBKGnumber
          << "\n" ;


//PG --------------- D R A W I N G ---------------------------


  H1fact.applyToAll (scaleH1 (0, norm40)) ;
  H1fact.applyToAll (scaleH1 (1, norm40)) ;
  H1fact.applyToAll (scaleH1 (2, norm40)) ;
  H1fact.applyToAll (scaleH1 (3, norm40)) ;
  H1fact.applyToAll (scaleH1 (4, norm40)) ;
  H1fact.applyToAll (scaleH1 (5, norm40)) ;
  H1fact.applyToAll (scaleH1 (6, norm22)) ;
  H1fact.applyToAll (scaleH1 (7, norm22)) ;
  H1fact.applyToAll (scaleH1 (8, norm22)) ;

  H2fact.applyToAll (scaleH2 (0, norm40)) ;
  H2fact.applyToAll (scaleH2 (1, norm40)) ;
  H2fact.applyToAll (scaleH2 (2, norm40)) ;
  H2fact.applyToAll (scaleH2 (3, norm40)) ;
  H2fact.applyToAll (scaleH2 (4, norm40)) ;
  H2fact.applyToAll (scaleH2 (5, norm40)) ;
  H2fact.applyToAll (scaleH2 (6, norm22)) ;
  H2fact.applyToAll (scaleH2 (7, norm22)) ;
  H2fact.applyToAll (scaleH2 (8, norm22)) ;

//  int colors[5] = {92,98,1,75,9} ;
  int colors[9] = {50,57,64,71,78, 85,92,99,106} ;
  H1fact.applyToAll (setH1Colors (std::vector<int> (colors, colors + 9))) ;

  H1fact["Z_shape"]->Print (1,2,"Z_shape_log") ;
  H1fact["tags_zep"]->Print (0,2,"tags_zep_bin") ;
  H1fact["quarks_shape"]->Print (1,2,"quarks_shape_log") ;
  H1fact["total_shape"]->Print (1,2,"total_shape_log") ;
  H1fact.Print (0,4) ;

//  H1fact.applyToAll (printEachH1 (0,1)) ;
  
  int colors2[9] = {50,57,64,71,78, 85,92,99,106} ;
  H2fact.applyToAll (setH2Colors (std::vector<int> (colors2, colors2 + 9))) ;
  H2fact.Print () ;
  H2fact.applyToAll (printEachH2 ()) ;

  // Now we are done.
  return 0 ;
}
