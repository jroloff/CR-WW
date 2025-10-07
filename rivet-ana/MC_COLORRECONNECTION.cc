// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#include <string.h>

/*

A RIVET routine for trying to understand what measurements would be useful for color reconnection
List of current included observables:

* Color reconnection studies

Written by Jennifer Roloff <jroloff2@gmail.com>

For more suggestions, please email Jennifer Roloff

*/

namespace Rivet {


  /// @brief MC_COLORRECONNECTION // placeholder number 1111111 for now
  class MC_COLORRECONNECTION : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_COLORRECONNECTION);
   
    // Given a list of bins, gives the bin which the value falls in
    // If it is below or above the last bin value, it will return -1
    int getBin(const double value, const std::vector<double> bins){
      for(unsigned int i=0; i<bins.size()-1; i++){
        if(value >= bins[i] && value < bins[i+1]) return i;
      }
      return -1;
    }



  
    //////////////////////////////////////////////////////////////////////////////////////////////////

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Note: if we want to get quarks and gluons before hadronization, can use FinalPartons

      // The basic final-state projection: all final-state particles within the given eta acceptance
      // 10 degrees < theta < 170 degrees
      const FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 0.1*GeV);
      declare(fs, "particles");

      ChargedFinalState tracks(Cuts::pT > 0.2*GeV && Cuts::abseta < m_maxEtaTracks);
      declare(tracks, "tracks");

      // TODO: I'm not sure why the Durham algorithm seems to be failing, but we should figure that out...
      //FastJets jetDurham(fs, FastJets::Algo::DURHAM, m_jetRadius, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      FastJets jetDurham(fs, FastJets::Algo::ANTIKT, m_jetRadius, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetDurham, "DurhamJets");



      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Book histograms
     
      // One histogram for each of the two jet pairs
      Histo1DPtr tmpDijetMass1;
      Histo1DPtr tmpDijetMass2;
      std::string nameDijetMass1 = "dijet_mass_1";
      std::string nameDijetMass2 = "dijet_mass_2";
      book(_h_dijetMass1, nameDijetMass1,  40, 0, 120);
      book(_h_dijetMass2, nameDijetMass2,  40, 0, 120);
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Most inspiration taken from https://arxiv.org/pdf/0704.0597
      // Some plots that would be good to add:
      // A cutflow of what event selections are leading to events being removed
      // A plot of the ycut of all durjet jets (before any cuts)
      // A plot of the number of selected jets (before any cuts)
      // A plot of the particle multiplicity (before any cuts)
      // A plot of the particle multiplicity (after all cuts)
      // A plot like Figure 6 (will probably be challenging), for all 4 regions
      


      // Things we don't need right now, but might want later
      //double collisionE = sqrtS();
      //const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();
      //const Particles& particles = apply<FinalState>(event, "particles").particles();
      

      const FastJets& durjetAlg = apply<FastJets>(event, "DurhamJets");
      const std::vector<fastjet::PseudoJet> durjet = durjetAlg.pseudojetsByPt();

      std::vector<fastjet::PseudoJet> selectedJets;
      for(unsigned int i=0; i<durjet.size(); i++){
        PseudoJet jj, j1, j2;
        jj = durjet[i];
        if(jj.has_parents(j1, j2)){
          FourMomentum fv1 = FourVector(j1.E(), j1.px() , j1.py(), j1.pz());
          FourMomentum fv2 = FourVector(j2.E(), j2.px() , j2.py(), j2.pz());
          double ycut = 2 * j2.E()*j2.E() * (1- std::cos(fv1.angle(fv2))) / (jj.E()*jj.E());
          if(ycut < 0.01) {
            selectedJets.push_back(jj);
          }
        }
      }


      // Need exactly 4 jets
      if (selectedJets.size() != 4) vetoEvent;

      // Two smallest angles < 100 degrees
      int minJetPair_1_jet1index = -1;
      int minJetPair_1_jet2index = -1;

      int minJetPair_2_jet1index = -1;
      int minJetPair_2_jet2index = -1;
 
      double minJetAngle1 = 1000;
      double minJetAngle2 = 1000;

      for(unsigned int i=0; i<selectedJets.size(); i++){
        for(unsigned int j=i; j<selectedJets.size(); j++){
          if(i==j) continue;
          FourMomentum fv1 = FourVector(selectedJets[i].E(), selectedJets[i].px() , selectedJets[i].py(), selectedJets[i].pz());
          FourMomentum fv2 = FourVector(selectedJets[j].E(), selectedJets[j].px() , selectedJets[j].py(), selectedJets[j].pz());
          double angleInDegrees = fv1.angle(fv2) * 180 / 3.14;

          if(angleInDegrees < minJetAngle1){
            minJetPair_2_jet1index = minJetPair_1_jet1index;
            minJetPair_2_jet2index = minJetPair_1_jet2index;
            minJetAngle2 = minJetAngle1;

            minJetPair_1_jet1index = i;
            minJetPair_1_jet2index = j;
            minJetAngle1 = angleInDegrees;
          }
          else if(angleInDegrees < minJetAngle2){
            minJetPair_2_jet1index = i;
            minJetPair_2_jet2index = j;
            minJetAngle2 = angleInDegrees;
          }
          
        }
      }

      // Find two other angles that are between 100 and 140 degrees, and aren't adjacent
      std::vector<double> jetAngles;
      std::vector<double> jetMasses;
      std::vector<int> jetIndices1;
      std::vector<int> jetIndices2;

      for(unsigned int i=0; i<selectedJets.size(); i++){
        for(unsigned int j=i+1; j<selectedJets.size(); j++){
          if( (i == minJetPair_1_jet1index && j == minJetPair_1_jet2index) || 
              (j == minJetPair_1_jet1index && i == minJetPair_1_jet2index) ||
              (i == minJetPair_2_jet1index && j == minJetPair_2_jet2index) || 
              (j == minJetPair_2_jet1index && i == minJetPair_2_jet2index) ) continue;

          FourMomentum fv1 = FourVector(selectedJets[i].E(), selectedJets[i].px() , selectedJets[i].py(), selectedJets[i].pz());
          FourMomentum fv2 = FourVector(selectedJets[j].E(), selectedJets[j].px() , selectedJets[j].py(), selectedJets[j].pz());
          double angleInDegrees = fv1.angle(fv2) * 180 / 3.14;
          fastjet::PseudoJet jetPair1 = selectedJets[i] + selectedJets[j];

          jetIndices1.push_back(i);
          jetIndices2.push_back(j);
          jetMasses.push_back(jetPair1.m());
          jetAngles.push_back( angleInDegrees);
        }
      }

      if(minJetAngle1 > 100) vetoEvent;
      if(minJetAngle2 > 100) vetoEvent;


      // These two angles don't share a jet
      // Not sure if we should veto the event, or just find another pairing
      if( minJetPair_1_jet1index ==  minJetPair_2_jet1index || minJetPair_1_jet1index ==  minJetPair_2_jet2index ||
          minJetPair_1_jet2index ==  minJetPair_2_jet1index || minJetPair_1_jet2index ==  minJetPair_2_jet2index){
         vetoEvent;
      }

      bool isPass = false;
      int index1 = -1;
      int index2 = -1;
      int index3 = -1;
      int index4 = -1;

      for(unsigned int i=0; i<jetAngles.size(); i++){
        if(jetAngles[i] < 100 || jetAngles[i] > 140) continue;
        for(unsigned int j=i+1; j<jetAngles.size(); j++){
          if(jetAngles[j] < 100 || jetAngles[j] > 140) continue;

          if( jetIndices1[i] ==  jetIndices1[j] || jetIndices2[i] ==  jetIndices2[j] ||
              jetIndices1[i] ==  jetIndices2[j] || jetIndices2[i] ==  jetIndices1[j]) continue;

          isPass = true;
          index1 = jetIndices1[i];
          index2 = jetIndices2[i];
          index3 = jetIndices1[j];
          index4 = jetIndices2[j];
        }
      }

      if(!isPass) vetoEvent;


      fastjet::PseudoJet correctJetPair1 = selectedJets[index1] + selectedJets[index2];
      fastjet::PseudoJet correctJetPair2 = selectedJets[index3] + selectedJets[index4];

      _h_dijetMass1->fill(correctJetPair1.m());
      _h_dijetMass2->fill(correctJetPair2.m());

  
    }

    /// Normalise histograms etc., after the run
    void finalize() 
    {      
    }


    bool m_debug = false;

    // If you want to match LEP, you should use this value for the max track eta
    //const double m_maxEtaTracks = 1.738; // Using cos(theta) < 0.94, which this should correspond to

    // For these studies, currently studying the same eta cutoff for all paricles and for tracks
    // This choice is somewhat arbitrary, and we could probably go higher in eta with future detectors
    const double m_maxEtaTracks = 2.1; 
    const double m_maxEta = 2.1;

    const double m_jetRadius = 0.7;

    //std::vector<Histo1DPtr> _h_dijetMasses;
    Histo1DPtr _h_dijetMass1;
    Histo1DPtr _h_dijetMass2;

  };


  RIVET_DECLARE_PLUGIN(MC_COLORRECONNECTION);

}


