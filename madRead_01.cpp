// c++ -o madRead_01 `root-config --glibs --cflags` -lm madRead_01.cpp 

#include "LHEF.h"
#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>
#include <map>

using namespace std ;

int main () {

  TH1F njets ("njets", "njets", 5, 0, 5) ;
  TH1F njets_gt20 ("njets_gt20", "njets_gt20", 5, 0, 5) ;
  TH1F njets_gt30 ("njets_gt30", "njets_gt30", 5, 0, 5) ;

  // Open a stream connected to an event file:
  std::ifstream ifs("../events/3j/events.lhe");
//  std::ifstream ifs("../events/3j/test.lhe");

  // Create the Reader object:
  LHEF::Reader reader(ifs);

  // Print out the header information:
  std::cerr << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cerr << reader.initComments;

  // Print out the beam energies:
  cout << " --- --- --- --- --- --- --- --- --- --- --- ---\n" ;
  std::cout << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;
  cout << " --- --- --- --- --- --- --- --- --- --- --- ---\n" ;

  // Now loop over all events:
  long ieve = 0;
  while ( reader.readEvent() ) {

    ++ieve;
    if (ieve % 1000 == 0) cout << "reading " << ieve << " event" << endl ;
//    std::cout << "Event " << ieve << " contains "
//              << reader.hepeup.NUP << " particles." << std::endl;
              
    map <int, TLorentzVector> higgs ;      
    map <int, TLorentzVector> partons ;      
    //PG loop over particles in the event
    for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
      {
        if (reader.hepeup.ISTUP.at (iPart) == -1) continue ;
//        std::cout << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                  << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                  << "\n" ;

        TLorentzVector particle 
          (
            reader.hepeup.PUP.at (iPart).at (0), //PG px
            reader.hepeup.PUP.at (iPart).at (1), //PG py
            reader.hepeup.PUP.at (iPart).at (2), //PG pz
            reader.hepeup.PUP.at (iPart).at (3) //PG E
          ) ;
        if (abs (reader.hepeup.IDUP.at (iPart)) == 25)   //PG H
          higgs[iPart] = particle ;  
        else 
          partons[iPart] = particle ;  
      }
    int count20 = 0 ;
    int count30 = 0 ;
    for (map<int, TLorentzVector>::const_iterator iMap = partons.begin () ;
         iMap != partons.end () ; ++iMap)
       {
         if (iMap->second.Pt () > 20) ++count20 ;
         if (iMap->second.Pt () > 30) ++count30 ;
       }  
    njets_gt20.Fill (count20) ;
    njets_gt30.Fill (count30) ;
    njets.Fill (partons.size ()) ;
  }

  TFile out ("madRead_01.root","recreate") ;
  njets.Write () ;
  njets_gt20.Write () ;
  njets_gt30.Write () ;
  out.Close () ;

  return 0 ;
}
