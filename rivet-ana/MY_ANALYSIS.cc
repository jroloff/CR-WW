#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"


void init() {
  // Final hadrons (normal Rivet analysis)

  // Final partons (before hadronization)
}

namespace Rivet {

  class MY_ANALYSIS : public Analysis {
  public:

    /// Constructor
    MY_ANALYSIS() : Analysis("MY_ANALYSIS") {}

    /// Book histograms and projections
    void init() {
      // All final state particles
      VetoedFinalState hadrons(FinalState(Cuts::abseta < 4.4));
      declare(hadrons, "hadrons");

      declare(FinalState(), "FS");
      // Jet clustering with anti-kT, R=0.4
      // declare(FastJets(fs, JetAlg::ANTIKT, 0.4), "Jets");
      FastJets jets(hadrons, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jets, "Jets");
      // Book a histogram of jet pT
      book(_h_jetpt, "jet_pT", 50, 0.0, 200.0);
      book(_h_NJets, "nJets", 6, 1, 6);
      book(_h_deltaTheta,"jet_deltaTheta",50,0.0,M_PI);
      book(_h_deltaPhi,  "jet_deltaPhi",50,-M_PI,M_PI);
      book(_h_deltaR,    "jet_deltaR",50,0.0,6.0);

      book(_h_q_deltaTheta,"quark_deltaTheta",50,0.0,M_PI);
      book(_h_q_deltaPhi,  "quark_deltaPhi",50,-M_PI,M_PI);
      book(_h_q_deltaR,    "quark_deltaR",50,0.0,6.0);
    }

    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV);

      // Jet multiplicity
    size_t njets = jets.size();
    _h_NJets->fill(njets);
    for (size_t i = 0; i < jets.size(); ++i) {
      // Jet kinematics 
      _h_jetpt->fill(jets[i].pT()/GeV);
      for (size_t j = i+1; j < jets.size(); ++j) {
        // diJet kinematics 
        double theta1 = jets[i].momentum().theta();
        double theta2 = jets[j].momentum().theta();
        double deltaTheta = std::abs(theta1 - theta2);
        _h_deltaTheta->fill(deltaTheta);

        double phi1 = jets[i].momentum().phi();
        double phi2 = jets[j].momentum().phi();
        double deltaPhi = phi1 - phi2;
        while (deltaPhi > M_PI)  deltaPhi -= 2*M_PI;
        while (deltaPhi < -M_PI) deltaPhi += 2*M_PI;
        _h_deltaPhi->fill(deltaPhi);

        double deltaEta = jets[i].momentum().eta() - jets[j].momentum().eta();
        double deltaR = std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
        _h_deltaR->fill(deltaR);

      }
    }
    const Particles& partons = apply<FinalState>(event, "FS").particles();


    // Collect all final-state quarks
    Particles quarks;
    for (const Particle& p : partons)
      if (PID::isQuark(p.pid())) quarks.push_back(p);

    // Loop over all quark pairs
    for (size_t i = 0; i < quarks.size(); ++i) {
      for (size_t j = i+1; j < quarks.size(); ++j) {

        double theta1 = quarks[i].momentum().theta();
        double theta2 = quarks[j].momentum().theta();
        double deltaTheta = std::abs(theta1 - theta2);
        _h_q_deltaTheta->fill(deltaTheta);

        double phi1 = quarks[i].momentum().phi();
        double phi2 = quarks[j].momentum().phi();
        double deltaPhi = phi1 - phi2;
        while (deltaPhi > M_PI)  deltaPhi -= 2*M_PI;
        while (deltaPhi < -M_PI) deltaPhi += 2*M_PI;
        _h_q_deltaPhi->fill(deltaPhi);

        double deltaEta = quarks[i].momentum().eta() - quarks[j].momentum().eta();
        double deltaR = std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
        _h_q_deltaR->fill(deltaR);
      }
    }
  }
    /// Finalize (normalize etc.)
    void finalize() {
      normalize(_h_jetpt);
      normalize(_h_NJets);
      normalize(_h_deltaTheta);
      normalize(_h_deltaPhi);
      normalize(_h_deltaR);
      normalize(_h_q_deltaTheta);
      normalize(_h_q_deltaPhi);
      normalize(_h_q_deltaR);
    }

  private:
    Histo1DPtr _h_jetpt, _h_NJets, _h_deltaTheta, _h_deltaPhi, _h_deltaR, _h_q_deltaTheta, _h_q_deltaPhi, _h_q_deltaR;
  };

  // Correct plugin macro for Rivet ≥ 3.1
  RIVET_DECLARE_PLUGIN(MY_ANALYSIS);

}
