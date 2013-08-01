#include "LHEF.h"

void writerExample() {

  // Open a stream connected to an event file:
  std::ofstream ofs("some-event-file.dat");

  // Create the Writerer object:
  LHEF::Writer writer(ofs);

  // Add XML info to the header block:
  writer.headerBlock() << "<!-- Some XML code ... -->\n";

  // Add a comment to the init block (# marks will automatically be
  // inserted).
  writer.initComments() << "Produced by LHEF::Writer.";

  // Set the init information of the processes to be written.
  writer.heprup.NPRUP = 2;
  writer.heprup.resize();
  writer.heprup.LPRUP[0] = 81;
  writer.heprup.LPRUP[0] = 82;
  // etc...

  // Write out the header and init blocks.
  writer.init();

  // Loop over the events to be written.
  int neve = 100;
  for ( int ieve = 0; ieve < neve; ++ieve ) {

    // Generate an event here and fill information in writer.hepeup
    // equivalent of the /HEPEUP/ common block. Remember to resize its
    // vectors (using Write::resize()) to make room for the number of
    // particles needed. Approximately like this:
    //
    // writer.hepeup.NUP = npart;
    // writer.hepeup.resize();

    // Possibly add a comment to this event.
    if ( writer.hepeup.NUP >= 30 ) writer.eventComments() << "# SPECIAL";

    // Finally we write out the event:
    writer.writeEvent();

  }

}
