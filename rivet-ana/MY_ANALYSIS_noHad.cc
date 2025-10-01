// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include <cmath>

namespace Rivet {

class QuarkPairDifferences : public Analysis {
public:
  RIVET_DEFAULT_ANALYSIS_CTOR(QuarkPairDifferences);

  void init() {
    declare(FinalState(), "FS");

    book(_h_deltaTheta,"quark_deltaTheta",50,0.0,M_PI);
    book(_h_deltaPhi,  "quark_deltaPhi",50,-M_PI,M_PI);
    book(_h_deltaR,    "quark_deltaR",50,0.0,6.0);
  }

  void analyze(const Event& event) {
    const Particles& fs = apply<FinalState>(event, "FS").particles();

    // Collect all final-state quarks
    Particles quarks;
    for (const Particle& p : fs)
      if (PID::isQuark(p.pid())) quarks.push_back(p);

    // Loop over all quark pairs
    for (size_t i = 0; i < quarks.size(); ++i) {
      for (size_t j = i+1; j < quarks.size(); ++j) {

        double theta1 = quarks[i].momentum().theta();
        double theta2 = quarks[j].momentum().theta();
        double deltaTheta = std::abs(theta1 - theta2);
        _h_deltaTheta->fill(deltaTheta);

        double phi1 = quarks[i].momentum().phi();
        double phi2 = quarks[j].momentum().phi();
        double deltaPhi = phi1 - phi2;
        while (deltaPhi > M_PI)  deltaPhi -= 2*M_PI;
        while (deltaPhi < -M_PI) deltaPhi += 2*M_PI;
        _h_deltaPhi->fill(deltaPhi);

        double deltaEta = quarks[i].momentum().eta() - quarks[j].momentum().eta();
        double deltaR = std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
        _h_deltaR->fill(deltaR);
      }
    }

    MSG_INFO("Event " << numEvents() << " has " << quarks.size() << " quarks");
  }

  void finalize() {
    normalize(_h_deltaTheta);
    normalize(_h_deltaPhi);
    normalize(_h_deltaR);
  }

private:
  Histo1DPtr _h_deltaTheta, _h_deltaPhi, _h_deltaR;
};

RIVET_DECLARE_PLUGIN(QuarkPairDifferences);

}
