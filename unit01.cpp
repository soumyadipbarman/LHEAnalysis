#include "LHEF.h"
#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"

// c++ `root-config --glibs --cflags` unit01.cpp 

int main (int argc, char **argv) {

  TH1F h_z_eta ("h_z_eta","outgoing Z",100,-5,5) ;

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
  while ( reader.readEvent() ) 
    {

      ++ieve;
      if (ieve % 10000 == 0) std::cout << "event " << ieve << "\n" ;
  
      std::vector<int> leptons ;      
      std::vector<int> bosons ;      
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
          std::cout << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
                    << "\n" ;

          if (abs (reader.hepeup.IDUP.at (iPart)) == 23)   //PG Z
            {
              TLorentzVector particle 
                (
                  reader.hepeup.PUP.at (iPart).at (0), //PG px
                  reader.hepeup.PUP.at (iPart).at (1), //PG py
                  reader.hepeup.PUP.at (iPart).at (2), //PG pz
                  reader.hepeup.PUP.at (iPart).at (3) //PG E
                ) ;
//              std::cout << "\t\t\teta " << particle.Eta () << "\n" ;
              h_z_eta.Fill (particle.Eta ()) ;
            }

          //PG leptons
          if (abs (reader.hepeup.IDUP.at (iPart)) == 23 || //PG electron
              abs (reader.hepeup.IDUP.at (iPart)) == 24)   //PG muon
            {
              bosons.push_back (iPart) ;
            }

          //PG outgoing particle          
          if (reader.hepeup.ISTUP.at (iPart) ==1)
            {
              //PG leptons
              if (abs (reader.hepeup.IDUP.at (iPart)) == 11 || //PG electron
                  abs (reader.hepeup.IDUP.at (iPart)) == 13)   //PG muon
                {
                  leptons.push_back (iPart) ;
                  TLorentzVector particle 
                    (
                      reader.hepeup.PUP.at (iPart).at (0), //PG px
                      reader.hepeup.PUP.at (iPart).at (1), //PG py
                      reader.hepeup.PUP.at (iPart).at (2), //PG pz
                      reader.hepeup.PUP.at (iPart).at (3) //PG E
                    ) ;
//                  reader.hepeup.PUP.at (iPart).at (4), //PG M

                } //PG leptons
            } //PG outcoming particle
        } //PG loop over particles in the event

      //PG search for VBF Z production
      bool isVBF = true ;
      int momId = reader.hepeup.MOTHUP.at (leptons.at (0)).first - 1 ;

std::cout << "test bosons num " << bosons.size () << "\n" ;
std::cout << "test leps num " << leptons.size () << "\n" ;
std::cout << "test leps " << abs (reader.hepeup.IDUP.at (leptons.at (0))) 
          << " " << abs (reader.hepeup.IDUP.at (leptons.at (1))) << "\n" ;
std::cout << "test lep mom " << reader.hepeup.MOTHUP.at (leptons.at (0)).first
          << " " << reader.hepeup.MOTHUP.at (leptons.at (0)).second << "\n" ;
std::cout << "test leps moms " << reader.hepeup.MOTHUP.at (leptons.at (0)).first
          << " " << reader.hepeup.MOTHUP.at (leptons.at (1)).first << "\n" ;
std::cout << "test mom [" << momId << "] " << reader.hepeup.IDUP.at (momId)  << "\n" ;
std::cout << "test granny " << reader.hepeup.MOTHUP.at (momId).first 
          << " " << reader.hepeup.MOTHUP.at (momId).second << "\n" ;

      //PG  - controllo che i due leptoni abbiano lo stesso sapore
      if (leptons.size () != 2) isVBF = false ;
      else if (abs (reader.hepeup.IDUP.at (leptons.at (0))) !=  
               abs (reader.hepeup.IDUP.at (leptons.at (1)))) isVBF = false ;
      //PG  - controllo che ogni leptone abbia una sola madre
      else if (reader.hepeup.MOTHUP.at (leptons.at (0)).first !=  
                 reader.hepeup.MOTHUP.at (leptons.at (0)).second ||
               reader.hepeup.MOTHUP.at (leptons.at (1)).first !=  
                 reader.hepeup.MOTHUP.at (leptons.at (1)).second)  isVBF = false ;
      //PG  - controllo che i due leptoni abbiano la stessa madre
      else if (reader.hepeup.MOTHUP.at (leptons.at (0)).first !=  
               reader.hepeup.MOTHUP.at (leptons.at (1)).first)  isVBF = false ;
      //PG  - controllo che questa madre sia una Z
      else if (reader.hepeup.IDUP.at (momId) != 23)  isVBF = false ;
      //PG  - controllo che questa Z abbia due madri
      else if (reader.hepeup.MOTHUP.at (momId).first == reader.hepeup.MOTHUP.at (momId).second)  isVBF = false ;
      //PG  - controllo che siano bosoni vettori
      else if (reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (momId).first + 1) != 23 &&
               abs (reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (momId).first + 1)) != 24)  isVBF = false ;
      else if (reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (momId).second + 1) != 23 &&
               abs (reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (momId).second + 1)) != 24)  isVBF = false ;
      //PG  - controllo che ciascuno abbia una madre entrante ed una uscente
      else if (
        reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (reader.hepeup.MOTHUP.at (momId).first).first + 1) *
        reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (reader.hepeup.MOTHUP.at (momId).first).second + 1) != -1 ||
        reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (reader.hepeup.MOTHUP.at (momId).second).first + 1) *
        reader.hepeup.ISTUP.at (reader.hepeup.MOTHUP.at (reader.hepeup.MOTHUP.at (momId).second).second + 1) != -1) isVBF = false ;
      if (isVBF)
        {
          std::cout << "VBF event\n" ;
          //PG do something to calculate exclusive x-section
          ++VBFnumber ;
        } 

    } // Now loop over all events

  std::cout << "VBF NUMBER " << VBFnumber << "\n" ;
  TFile histosFile ("output.root","recreate") ;
  histosFile.cd () ;
  h_z_eta.Write () ;
  histosFile.Close () ;

  // Now we are done.
  return 0 ;
}
