#include "LHEF.h"

int main () {

  // Open a stream connected to an event file:
  std::ifstream ifs("../pietro_test4_events.lhe");

  // Create the Reader object:
  LHEF::Reader reader(ifs);

  // Create the Writer object:
  LHEF::Writer writer(std::cout);

  // Print out the header information:
  std::cerr << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cerr << reader.initComments;

  // Print out the beam energies:
  std::cerr << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;

  // Copy header and init blocks and write them out.
  if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
  writer.headerBlock() << reader.headerBlock;
  writer.initComments() << reader.initComments;
  writer.heprup = reader.heprup;
  writer.init();

  // Now loop over all events:
  long ieve = 0;
  while ( reader.readEvent() ) {

    ++ieve;

    // Some events may be specially tagged "# SPECIAL" in the comment
    // lines, in that case write out the number of particles in that
    // event:
    if ( reader.eventComments.find("# SPECIAL") != std::string::npos ) 
      std::cout << "Event " << ieve << " contained a special event with "
                << reader.hepeup.NUP << " particles." << std::endl;

  }

  // Now we are done.
  return 0 ;
}
