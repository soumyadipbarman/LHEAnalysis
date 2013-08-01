/*
 c++ -o vbf_test `root-config --glibs --cflags` vbf_test.cpp
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


struct mySort: public std::binary_function<TLorentzVector *, TLorentzVector *, bool>
{
  bool operator() (TLorentzVector * x, TLorentzVector * y)
    {
      return x->Pt () > y->Pt () ;
    }
} ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


void fillHistos (ifstream & ifs) 
{
  LHEF::Reader reader (ifs) ;
  TH1F h_deta ("h_deta", "h_deta", 100, 0, 10) ;
  TH1F h_jeta ("h_jeta", "h_jeta", 100, -5, 5) ;
  TH1F h_j1eta ("h_j1eta", "h_j1eta", 100, -5, 5) ;
  TH1F h_j2eta ("h_j2eta", "h_j2eta", 100, -5, 5) ;
  TH1F h_j1pt ("h_j1pt", "h_j1pt", 100, 0, 200) ;
  TH1F h_j2pt ("h_j2pt", "h_j2pt", 100, 0, 200) ;

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
      vector <TLorentzVector *> gluons ;

      for (int iPart = 0 ; iPart < fs_particles.size () ; ++iPart)
        {
          if (abs (fs_particles.at (iPart).first) < 11) //PG tag jets
            {
              tagJets.push_back (&(fs_particles.at (iPart).second)) ;
              continue ;
            }
          if (abs (fs_particles.at (iPart).first) == 21) //PG tag jets
            {
              gluons.push_back (&(fs_particles.at (iPart).second)) ;
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

      //PG VBF cuts
      sort (tagJets.begin (), tagJets.end (), mySort ()) ;
      double tj_dEta = fabs (tagJets.at (0)->Eta () - tagJets.at (1)->Eta ()) ; 
      h_deta.Fill (tj_dEta) ; 
      h_j1eta.Fill (tagJets.at (0)->Eta ()) ;
      h_j2eta.Fill (tagJets.at (1)->Eta ()) ;
      h_j1pt.Fill (tagJets.at (0)->Pt ()) ;
      h_j2pt.Fill (tagJets.at (1)->Pt ()) ;
      for (int k = 0 ; k < tagJets.size () ; ++k) h_jeta.Fill (tagJets.at (k)->Eta ()) ;
      if (gluons.size () > 0) h_jeta.Fill (gluons.at (0)->Eta ()) ;
        

    } // Now loop over all events

  TFile histosFile ("vbf_test.root","recreate") ;
  histosFile.cd () ;
  h_deta.Write () ;
  h_j1eta.Write () ;
  h_j2eta.Write () ;
  h_jeta.Write () ;
  h_j1pt.Write () ;
  h_j2pt.Write () ;

  histosFile.Close () ;

  return ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  ifstream ifs ;

  ifs.open ("../h_qqH_WW_2L2NU_170.lhe") ; fillHistos (ifs) ; ifs.close () ;
  
  // Now we are done.
  return 0 ;
}
