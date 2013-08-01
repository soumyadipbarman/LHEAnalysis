/*
 c++ -o ewk_ww_study `root-config --glibs --cflags` ewk_ww_study.cpp
*/

#include "LHEF.h"
#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <map>
#include <cmath>
#include <algorithm>


using namespace std ;



// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


void fillHistos (TH1F & h_H_mass_Phantom, ifstream & ifs) 
{
  LHEF::Reader reader (ifs) ;

  long ieve = 0 ;
  while ( reader.readEvent() ) 
    {
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
      ++ieve;
      
      std::vector<int> bosons ;      
      vector<pair<int, TLorentzVector> > fs_particles ;  //PG final state particles
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
          //PG outgoing particle          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              TLorentzVector particle 
                (
                  reader.hepeup.PUP.at (iPart).at (0), //PG px
                  reader.hepeup.PUP.at (iPart).at (1), //PG py
                  reader.hepeup.PUP.at (iPart).at (2), //PG pz
                  reader.hepeup.PUP.at (iPart).at (3) //PG E
                ) ;
//              fs_particles[reader.hepeup.IDUP.at (iPart)] = particle ;            
              fs_particles.push_back (pair<int, TLorentzVector> (reader.hepeup.IDUP.at (iPart), particle)) ;
            } //PG outgoing particle
        } //PG loop over particles in the event

      TLorentzVector Z_cand ;
      TLorentzVector H_cand ;
      vector <TLorentzVector *> tagJets ;
      vector <TLorentzVector *> leptons ;

      for (int iPart = 0 ; iPart < fs_particles.size () ; ++iPart)
        {
          if (abs (fs_particles.at (iPart).first) < 11) //PG tag jets
            {
              tagJets.push_back (&(fs_particles.at (iPart).second)) ;
              continue ;
            }
          if (abs (fs_particles.at (iPart).first) == 11 ||
              abs (fs_particles.at (iPart).first) == 13)
            {
               leptons.push_back (&(fs_particles.at (iPart).second)) ;
               Z_cand += fs_particles.at (iPart).second ; //PG electron or muon
            }
          H_cand += fs_particles.at (iPart).second ;
        }

      if (leptons.at (0)->Pt () < leptons.at (1)->Pt ()) swap (leptons.at (0), leptons.at (1)) ;
      if (leptons.at (0)->Pt () < 20 || leptons.at (1)->Pt () < 10) continue ;
      if (fabs (leptons.at (0)->Eta ()) > 2.5 || fabs (leptons.at (1)->Eta ()) > 2.5) continue ;
      if (tagJets.at (0)->Pt () < 30 || tagJets.at (1)->Pt () < 30) continue ;
      if (fabs (tagJets.at (0)->Eta ()) > 4.5 || fabs (tagJets.at (1)->Eta ()) > 4.5) continue ;

      //PG Z veto
      if (Z_cand.M () < 12 ||
          Z_cand.M () > 75) continue ;
          
      //PG VBF cuts
      double tj_dEta = fabs (tagJets.at (0)->Eta () - tagJets.at (1)->Eta ()) ;    
      if (tj_dEta < 3.5) continue ;

      TLorentzVector tj_pair = *(tagJets.at (0)) + *(tagJets.at (1)) ;
      double tj_mass = tj_pair.M () ;
      if (tj_mass < 450) continue ;
          
      h_H_mass_Phantom.Fill (H_cand.M ()) ;
    } // Now loop over all events

  return ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  TH1F h_H_mass_Phantom ("h_H_mass_Phantom", "higgs boson mass", 5000, 0, 1000) ;
  ifstream ifs ;

  ifs.open ("../genh170/gen1/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen2/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen3/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen4/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen5/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen6/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen7/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen8/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  ifs.open ("../genh170/gen9/phamom.dat") ; fillHistos (h_H_mass_Phantom, ifs) ; ifs.close () ;
  h_H_mass_Phantom.Scale (1./h_H_mass_Phantom.GetEntries ()) ;
  h_H_mass_Phantom.Scale (1.453042405755725E-002) ; //PG from the file res, is the average cross-section
  
  TFile histosFile ("ewk_ww_study.root","recreate") ;
  histosFile.cd () ;
  h_H_mass_Phantom.Write () ;
  histosFile.Close () ;

  // Now we are done.
  return 0 ;
}
