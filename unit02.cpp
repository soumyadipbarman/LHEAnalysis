/*
per dare massa ai leptoni, assumendo che non cambi l'angolo del momento:

E^2 - (pc)^2 = (mc^2)^2

m = 0 :  E^2 - (pc)^2 = 0  =>  E^2 = (pc)^2

se aggiungo massa non nulla, il p scende:

(p'c)^2 = E^2 - (mc^2)^2

=> (p')^2 = (E/c)^2 - (mc)^2

quindi al momento attribuisco il raggio con il metodo :

   TLorentzVector::SetRho(10.); 

(rimane da capire che unita' di misura utilizza LHA)
(rimane da riscrivere l'output)

*/


#include <cmath>
#include "LHEF.h"
#include "TLorentzVector.h"

// c++ -o unit02 `root-config --glibs --cflags` -lm unit02.cpp 

int main (int argc, char **argv) {

  // Open a stream connected to an event file:
  if (argc < 2) exit (1) ;
  std::ifstream ifs(argv[1]);

  // Create the Reader object:
  LHEF::Reader reader(ifs);
  LHEF::Writer writer(std::cout);

  // Copy header and init blocks and write them out.
  if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
  writer.headerBlock() << reader.headerBlock;
  writer.initComments() << reader.initComments;
  writer.heprup = reader.heprup;
  writer.init();

  // Print out the header information:
  std::cerr << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cerr << reader.initComments;

  // Print out the beam energies:
  std::cerr << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;

  // Now loop over all events:
  long ieve = 0;
  long VBFnumber = 0;
  while ( reader.readEvent() ) 
    {
      ++ieve;
      if (ieve % 10000 == 0) std::cerr << "event " << ieve << "\n" ;
  
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
//          std::cerr << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                    << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                    << "\n" ;

          //PG leptons
          if (abs (reader.hepeup.IDUP.at (iPart)) == 11 || //PG electron
              abs (reader.hepeup.IDUP.at (iPart)) == 13)   //PG muon
            {
              TLorentzVector particle 
                (
                  reader.hepeup.PUP.at (iPart).at (0), //PG px
                  reader.hepeup.PUP.at (iPart).at (1), //PG py
                  reader.hepeup.PUP.at (iPart).at (2), //PG pz
                  reader.hepeup.PUP.at (iPart).at (3) //PG E
                ) ;
//                  reader.hepeup.PUP.at (iPart).at (4), //PG M
              double newMomentum = //PG in case c = 1
                reader.hepeup.PUP.at (iPart).at (3) * 
                    reader.hepeup.PUP.at (iPart).at (3) -
                reader.hepeup.PUP.at (iPart).at (4) * 
                    reader.hepeup.PUP.at (iPart).at (4) ;
                    
              particle.SetRho (sqrt (newMomentum)) ;

//                  std::cerr << " test resizing "
//                            << reader.hepeup.PUP.at (iPart).at (3)
//                            << " " << particle.E ()
//                            << "\n" ;

              reader.hepeup.PUP.at (iPart).at (0) = particle.X () ; //PG px
              reader.hepeup.PUP.at (iPart).at (1) = particle.Y () ; //PG py
              reader.hepeup.PUP.at (iPart).at (2) = particle.Z () ; //PG pz

            } //PG leptons

        } //PG loop over particles in the event

      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
      writer.eventComments() << reader.eventComments;
      writer.hepeup = reader.hepeup;
      writer.writeEvent();  

    } // Now loop over all events

  std::cerr << "VBF NUMBER " << VBFnumber << "\n" ;

  // Now we are done.
  return 0 ;
}
