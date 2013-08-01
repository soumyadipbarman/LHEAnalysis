/*
 c++ -o readPowheg `root-config --glibs --cflags` readPowheg.cpp
superseeded by readNtuple in parent folder

*/

#include "LHEF.h"
#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <map>
#include <cmath>


using namespace std ;



// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


void fillHistosPowheg (TH1F & h_H_mass_Powheg, ifstream & ifs) 
{
  LHEF::Reader reader (ifs) ;

  long ieve = 0 ;
  while ( reader.readEvent() ) 
    {
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
      ++ieve;
      
      std::vector<int> leptons ;      
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

      TLorentzVector H_cand ;
      vector <TLorentzVector *> tagJets ;

      for (int iPart = 0 ; iPart < fs_particles.size () ; ++iPart)
        {
          if (abs (fs_particles.at (iPart).first) < 11) //PG tag jets
            {
              tagJets.push_back (&(fs_particles.at (iPart).second)) ;
              continue ;
            }
          if (abs (fs_particles.at (iPart).first) == 25) H_cand = fs_particles.at (iPart).second ;
        }

      //PG VBF cuts
      double tj_dEta = fabs (tagJets.at (0)->Eta () - tagJets.at (1)->Eta ()) ;    
      if (tj_dEta < 3.5) continue ;

      TLorentzVector tj_pair = *(tagJets.at (0)) + *(tagJets.at (1)) ;
      double tj_mass = tj_pair.M () ;
      if (tj_mass < 450) continue ;
      
      h_H_mass_Powheg.Fill (H_cand.M ()) ;
    } // Now loop over all events

  //PG FIXME i tagli applicati da alessandro vanno messi anche qui
  /*  
  pt(l,j,miss)>20 GeV, M(jj)>100 GeV, M(ll)>20 GeV,
  eta(j)<6.5, eta(l)<3, delta_eta(jj)>2
  */

  return ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  TH1F h_H_mass_Phantom ("h_H_mass_Phantom", "higgs boson mass", 5000, 0, 1000) ;
  ifstream ifs ;

  TH1F h_H_mass_Powheg ("h_H_mass_Powheg", "higgs boson mass", 5000, 0, 1000) ;
  ifs.open ("../h_qqH_WW_2L2NU_170.lhe") ; fillHistosPowheg (h_H_mass_Powheg, ifs) ; ifs.close () ;

  cout << h_H_mass_Powheg.GetEntries () << endl ;

  h_H_mass_Powheg.Scale (1./h_H_mass_Powheg.GetEntries ()) ;
  h_H_mass_Powheg.Scale (8.18901E-01 * 1.14E-02) ; //PG prod XS from powheg, H BR from XS WG
  //PG https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR

  TFile histosFile ("readPowheg.root","recreate") ;
  histosFile.cd () ;
  h_H_mass_Powheg.Write () ;
  histosFile.Close () ;

  // Now we are done.
  return 0 ;
}
