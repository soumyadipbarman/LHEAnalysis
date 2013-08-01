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

  c++ -o contamin01 `root-config --glibs --cflags` \
     -lm hFactory.cc hChain.cc h2Factory.cc h2Chain.cc contamin01.cpp
     
TO RUN:     

  ./contamin01 \
  /Users/govoni/private/job/cms/HIGGS/POWHEG_BoX/ggF/ggH170_WW_lnulnu_CTEQ6m_7TeV/pwgevents.lhe \
  /Users/govoni/private/job/cms/HIGGS/POWHEG_BoX/VBF/qqH170_WW_lnulnu_7TeV/pwgevents.lhe


*/



void fillHistos (std::string fileName, TNtuple & ntuple)
{
  std::ifstream ifs (fileName.c_str ()) ;
  LHEF::Reader reader (ifs) ;

  long ieve = 0 ;
  double jetsNum = 0. ;
  
  // loop over events
  while ( reader.readEvent () ) 
    {
      ++ieve;
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
  
      std::vector<int> finalJets ;      
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
              //PG quarks
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7 ||
                  reader.hepeup.IDUP.at (iPart) == 21)
                {
                  finalJets.push_back (iPart) ;
                }

            } //PG outgoing particles
        } //PG loop over particles in the event

      jetsNum += finalJets.size () ;
      if (finalJets.size () < 2) continue ;

      //PG build the two quarks
      TLorentzVector fquark0
        (
          reader.hepeup.PUP.at (finalJets.at (0)).at (0), //PG px
          reader.hepeup.PUP.at (finalJets.at (0)).at (1), //PG py
          reader.hepeup.PUP.at (finalJets.at (0)).at (2), //PG pz
          reader.hepeup.PUP.at (finalJets.at (0)).at (3) //PG E
        ) ;
      
      TLorentzVector fquark1
        (
          reader.hepeup.PUP.at (finalJets.at (1)).at (0), //PG px
          reader.hepeup.PUP.at (finalJets.at (1)).at (1), //PG py
          reader.hepeup.PUP.at (finalJets.at (1)).at (2), //PG pz
          reader.hepeup.PUP.at (finalJets.at (1)).at (3) //PG E
        ) ;
      //PG the sum pf the two quarks
      TLorentzVector vCand = fquark0 + fquark1 ;
       ntuple.Fill (
           vCand.M () ,
           fabs (fquark0.Eta () - fquark1.Eta ()) ,
           fquark0.DeltaR (fquark1) ,
           fabs (fquark0.DeltaPhi (fquark1)) 
         ) ;
   
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

    } // loop over events
  std::cout << jetsNum/ieve << " in " << fileName << std::endl ;

}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char **argv) 
{

  gROOT->SetStyle ("Plain") ;

  // Open a stream connected to an event file:
  if (argc < 3) exit (1) ;

  TNtuple ggFntuple ("ggFntuple", "ggFntuple", "Minv:Deta:DR:Dphi") ;
  fillHistos (argv[1], ggFntuple) ;
  TNtuple vbFntuple ("vbFntuple", "vbFntuple", "Minv:Deta:DR:Dphi") ;
  fillHistos (argv[2], vbFntuple) ;

  int bins = 50 ;
  TH2F Minv_vs_deta_cont ("Minv_vs_deta_cont", "contamination in VBF selections", 
                          bins, 0, 1500, bins, 0, 5) ;
  Minv_vs_deta_cont.SetStats (0) ;                          
  Minv_vs_deta_cont.GetXaxis ()->SetTitle ("M_{inv}") ;                          
  Minv_vs_deta_cont.GetYaxis ()->SetTitle ("#Delta#eta") ;                          

  double num_ggF = ggFntuple.GetEntries () ;
  double num_vbF = vbFntuple.GetEntries () ;
  double ggF_Xsec = 14.761 ; // pb NNLO
  double vbF_Xsec = 1.6864 ; // pb NLO
 
  for (double iMinv = 0.5 * (1500./bins) ; iMinv < 1500. ; iMinv += 1500./bins)
    for (double iDeta = 0.5 * (5./bins) ; iDeta < 5. ; iDeta += 5./bins)
      {
        char cut[20] ;
        sprintf (cut, "Minv>%f&&Deta>%f", iMinv, iDeta) ;
//        std::cout << cut << std::endl ;
        double num_ggF_sel = ggFntuple.GetEntries (cut) ;
        double num_vbF_sel = vbFntuple.GetEntries (cut) ;
        double val = (ggF_Xsec * num_ggF_sel / num_ggF) /
          (ggF_Xsec * num_ggF_sel / num_ggF + vbF_Xsec * num_vbF_sel / num_vbF) ;
        Minv_vs_deta_cont.Fill (iMinv, iDeta, val) ;
      }

  std::cout << "out of loop" << std::endl ;

  TCanvas c1 ;
  Minv_vs_deta_cont.Draw ("COLZ") ;
  c1.Print ("Minv_vs_deta_cont.gif","gif") ;

  TFile output ("ouptut.root","recreate") ;
  output.cd () ;
  Minv_vs_deta_cont.Write () ;
  output.Close () ;
    
  return 0 ;
}
