#include "LHEF.h"
#include <iomanip>

int main(int argc, char ** argv) {

  // Create Reader and Writer object
  LHEF::Reader reader(std::cin);
  LHEF::Writer writer(std::cout);

  // Copy header and init blocks and write them out.
  if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
  writer.headerBlock() << reader.headerBlock;
  writer.initComments() << reader.initComments;
  writer.heprup = reader.heprup;
  writer.init();

  // Read each event and write them out again.
  while ( reader.readEvent() ) {
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    writer.eventComments() << reader.eventComments;
    writer.hepeup = reader.hepeup;
    writer.writeEvent();
  }

}
