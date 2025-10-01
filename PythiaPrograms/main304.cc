// main302.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: colour reconnection; e+e- events;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Plot the colour reconnection rate in W+W- hadronic events
// as a function of collision energy.

// Pythia includes.

#include "Pythia8/Pythia.h"

// Preferably use HepMC3, but alternatively HepMC2.
#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif

using namespace Pythia8;

// Preferably use HepMC3, but alternatively HepMC2.


//==========================================================================

int main() {

  // Number of events.
  int nEvent = 4000;
  // Histograms.
  Pythia8ToHepMC toHepMC("main131.hepmc");
 
  // Loop over reconnection model and CM energy.
      // Quiet generator creation. Shorthand for the event record.
  Pythia pythia ("../share/Pythia8/xmldoc", false);

  // Incoming beams and energy. Optionally no ISR.
  pythia.readString("Beams:idA = -11");
  pythia.readString("Beams:idB = 11");
  pythia.settings.parm("Beams:eCM", 240.);
  pythia.readString("PDF:lepton = off");

  // Hard process, with hadronic W decays.
  pythia.readString("WeakDoubleBoson:ffbar2WW = on");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = 1 2 3 4 5");

  // Switch on CR.
  pythia.readString("ColourReconnection:reconnect = on");
  pythia.settings.mode("ColourReconnection:mode", 1);
     
  pythia.readString("ColourReconnection:forceResonance = on");
  // Reduce printout. Switch off hadronization and decays.
  pythia.readString("Print:quiet = on");
  pythia.readString("HadronLevel:Hadronize = off");
  pythia.readString("HadronLevel:Decay = off");

   // If Pythia fails to initialize, exit with error.
   if (!pythia.init()) return 1;


  // Begin event loop. Generate event. Skip if error.
   for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
	if (!pythia.next()) continue;

	toHepMC.writeNextEvent( pythia );
    // Find number of all final charged particles and fill histogram.
    // Construct new empty HepMC event, fill it and write it out.

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();


  // Done.
  return 0;

}
