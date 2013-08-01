/*
 c++ -o madweight_test `root-config --glibs --cflags` madweight_test.cpp
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
  TH1F h_nlep ("h_nlep", "h_nlep", 10, 0, 10) ;
  TH1F h_nneu ("h_nneu", "h_nneu", 10, 0, 10) ;
  TH1F h_njets ("h_njets", "h_njets", 10, 0, 10) ;
  TH1F h_deta ("h_deta", "h_deta", 100, 0, 10) ;
  TH1F h_jeta ("h_jeta", "h_jeta", 100, -5, 5) ;
  TH1F h_leppt ("h_leppt", "h_leppt", 100, 0, 200) ;
  TH1F h_mH ("h_mH", "h_mH", 100, 0, 1000) ;
  TH1F h_lepeta ("h_lepeta", "h_lepeta", 100, -5, 5) ;
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
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size () ; ++iPart)
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

      vector <TLorentzVector *> quarks ;
      vector <TLorentzVector *> leptons ;
      vector <TLorentzVector *> neutrinos ;
      vector <TLorentzVector *> gluons ;

      for (int iPart = 0 ; iPart < fs_particles.size () ; ++iPart)
        {
          if (abs (fs_particles.at (iPart).first) < 11) //PG tag jets
            {
              quarks.push_back (&(fs_particles.at (iPart).second)) ;
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
            }
          if (abs (fs_particles.at (iPart).first) == 12 ||
              abs (fs_particles.at (iPart).first) == 14)
            {
               neutrinos.push_back (&(fs_particles.at (iPart).second)) ;
            }
        }

      //PG VBF cuts
      sort (quarks.begin (), quarks.end (), mySort ()) ;
      double tj_dEta = fabs (quarks.at (0)->Eta () - quarks.at (1)->Eta ()) ; 
      h_deta.Fill (tj_dEta) ; 
      h_j1eta.Fill (quarks.at (0)->Eta ()) ;
      h_j2eta.Fill (quarks.at (1)->Eta ()) ;
      h_j1pt.Fill (quarks.at (0)->Pt ()) ;
      h_j2pt.Fill (quarks.at (1)->Pt ()) ;
      for (int k = 0 ; k < quarks.size () ; ++k) h_jeta.Fill (quarks.at (k)->Eta ()) ;
      h_njets.Fill (quarks.size ()) ;
      h_nneu.Fill (neutrinos.size ()) ;
      h_nlep.Fill (leptons.size ()) ;
      for (int k = 0 ; k < leptons.size () ; ++k) 
        {
          h_lepeta.Fill (leptons.at (k)->Eta ()) ;
          h_leppt.Fill (leptons.at (k)->Pt ()) ;
        }
      if (gluons.size () > 0) h_jeta.Fill (gluons.at (0)->Eta ()) ;

      TLorentzVector H_cand ;
      H_cand += *quarks.at (0) ;
      H_cand += *quarks.at (1) ;
      H_cand += *leptons.at (0) ;
      H_cand += *neutrinos.at (0) ;

      h_mH.Fill (H_cand.M ()) ;  

    } // Now loop over all events

  TFile histosFile ("madweight_test.root","recreate") ;
  histosFile.cd () ;
  h_deta.Write () ;
  h_j1eta.Write () ;
  h_j2eta.Write () ;
  h_jeta.Write () ;
  h_j1pt.Write () ;
  h_j2pt.Write () ;
  h_lepeta.Write () ;
  h_leppt.Write () ;
  h_njets.Write () ;
  h_nneu.Write () ;
  h_nlep.Write () ;
  h_mH.Write () ;

  histosFile.Close () ;

  return ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  ifstream ifs ;

  ifs.open ("../test_004_unweighted_events.lhe") ; fillHistos (ifs) ; ifs.close () ;
  
  // Now we are done.
  return 0 ;
}
