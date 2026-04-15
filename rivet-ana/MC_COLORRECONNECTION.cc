// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include <fstream>
//#include "Rivet/Logging.hh"

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

      book(quarkMassWm, "Wm_Quark_Mass",  40, 0, 120);
      book(quarkMassWp, "Wp_Quark_Mass",  40, 0, 120);

      // One histogram for tracking what particles are left after cuts
      std::string cutflow = "Cutflow";
      book(_h_cutflow, cutflow, 5, 0.5, 5.5);

      //ycut before event vetoes
      book(_h_ycutAll, "ycutAll", 60, 1e-6, 0.5);

      //Selected jets before cut
      book(_h_selectedJets, "selectedJets", 13, -0.5, 12.5); 

      //particle multiplicity before cuts
      book(_h_mult_before, "PreMultiplicity", 151, -0.5, 150.5);   

      //particle multiplicity after cuts
      book(_h_mult_after, "PostMultiplicity", 151, -0.5, 150.5);

      //Figure 5
      book(phi_rescaled, "ParticleFlow_PhiRescaled", 84, -0.2, 4.2);

      //Figure 6
      book(insideTotal, "Inside_W_Region", 36, -0.2, 1.2);
      book(OutsideTotal, "Outside_W_Region", 36, -0.2, 1.2);
      book(region_ratio, "Region_Ratio", 36, -0.2, 1.2);

      //Diagnostic Plots
      book(_h_mult_nocut, "Multiplicity_before_cuts", 150, -0.5, 150.5);

      //ThetaRescaled for each region
      book(inside1, "inside1", 24, -0.2, 1.2);
      book(inside2, "inside2", 24, -0.2, 1.2);
      book(outside1, "outside1", 24, -0.2, 1.2);
      book(outside2, "outside2", 24, -0.2, 1.2);

      book(outside_ratio, "Outside_Ratio", 24, -0.2, 1.2);
      book(inside_ratio, "Inside_Ratio", 24, -0.2, 1.2);

      //Angle between each region
      book(regionA, "Region_A", 72, -0.5, 180.5);
      book(regionB, "Region_B", 72, -0.5, 180.5);
      book(regionC, "Region_C", 72, -0.5, 180.5);
      book(regionD, "Region_D", 72, -0.5, 180.5);

      book(quark1, "Quark1_deltaR", 50, 0, 1.5);
      book(quark2, "Quark2_deltaR", 50, 0, 1.5);
      book(quark3, "Quark3_deltaR", 50, 0, 1.5);
      book(quark4, "Quark4_deltaR", 50, 0, 1.5);

      book(mismatchDeltaR, "Matching_DeltaR", 50, 0, 3);

      book(check_index, "check_index", 6, -1, 5);

      book(jetToQuark1, "pT_Ratio_1", 50, 0, 7);
      book(jetToQuark2, "pT_Ratio_2", 50, 0, 7);
      book(jetToQuark3, "pT_Ratio_3", 50, 0, 7);
      book(jetToQuark4, "pT_Ratio_4", 50, 0, 7);

      book(quarkdR1, "quark_region_dR1", 50, 0, 1.5);
      book(quarkdR2, "quark_region_dR2", 50, 0, 1.5);
      book(quarkdR3, "quark_region_dR3", 50, 0, 1.5);
      book(quarkdR4, "quark_region_dR4", 50, 0, 1.5);

      //Using truth quarks
      book(regionQA, "Quark_Region_A", 72, -0.5, 180.5);
      book(regionQB, "Quark_Region_B", 72, -0.5, 180.5);
      book(regionQC, "Quark_Region_C", 72, -0.5, 180.5);
      book(regionQD, "Quark_Region_D", 72, -0.5, 180.5);

      book(sorted, "Sorted_particles", 13, -0.5, 12.5);
      

      out = new std::ofstream("outEventDisplay.csv");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Most inspiration taken from https://arxiv.org/pdf/0704.0597

      // Things we don't need right now, but might want later
      //double collisionE = sqrtS();
      //const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();
      const Particles& particles = apply<FinalState>(event, "particles").particles();

      //All events
      //Use event.weight() maybe to normalize
      _h_cutflow->fill(1);

      _h_mult_nocut->fill(particles.size());

      // Fastjet analysis - select algorithm and parameters
      fastjet::Strategy               strategy = fastjet::Best;
      fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
      fastjet::JetDefinition*         jetDef = new fastjet::JetDefinition(fastjet::ee_kt_algorithm,
                                      recombScheme, strategy);

      // Run Fastjet algorithm
      fastjet::ClusterSequence clustSeq(particles, *jetDef);
      int nJetMin = 4;
      //clustSeq.exclusive_dmerge_max(nJetMin-1);

      const std::vector<fastjet::PseudoJet> durjet = clustSeq.exclusive_jets(nJetMin);

      _h_mult_before->fill(particles.size());

      std::vector<fastjet::PseudoJet> selectedJets;

      //ycut calculation
      for(unsigned int i=0; i<durjet.size(); i++){
        PseudoJet jj, j1, j2;
        jj = durjet[i];
        if(jj.has_parents(j1, j2)){
          FourMomentum fv1 = FourVector(j1.E(), j1.px() , j1.py(), j1.pz());
          FourMomentum fv2 = FourVector(j2.E(), j2.px() , j2.py(), j2.pz());
          double ycut = 2 * j2.E()*j2.E() * (1- std::cos(fv1.angle(fv2))) / (jj.E()*jj.E());

          //Add ycut pre-veto
          _h_ycutAll->fill(ycut);

          if(ycut < 0.005) {
            selectedJets.push_back(jj);
          }
        }
      }

      //Histo of selected jets
      _h_selectedJets->fill(durjet.size());

      // Need exactly 4 jets
      if (selectedJets.size() != 4) vetoEvent;

      //Add event that survived 4-jet selection
      _h_cutflow->fill(2);

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
      //These become the jet pairs from the W
      std::vector<double> jetAngles;
      //std::vector<double> jetMasses;
      std::vector<int> jetIndices1;
      std::vector<int> jetIndices2;

      for(unsigned int i=0; i<selectedJets.size(); i++){
        for(unsigned int j=i+1; j<selectedJets.size(); j++){
          //The selected jet pairs aren't the smallest pairs
          if( (i == minJetPair_1_jet1index && j == minJetPair_1_jet2index) ||
              (j == minJetPair_1_jet1index && i == minJetPair_1_jet2index) ||
              (i == minJetPair_2_jet1index && j == minJetPair_2_jet2index) ||
              (j == minJetPair_2_jet1index && i == minJetPair_2_jet2index) ) continue;

          FourMomentum fv1 = FourVector(selectedJets[i].E(), selectedJets[i].px() , selectedJets[i].py(), selectedJets[i].pz());
          FourMomentum fv2 = FourVector(selectedJets[j].E(), selectedJets[j].px() , selectedJets[j].py(), selectedJets[j].pz());
          double angleInDegrees = fv1.angle(fv2) * 180 / 3.14;
          //fastjet::PseudoJet jetPair1 = selectedJets[i] + selectedJets[j];

          //Collect the 4 remaining jet combinations and their angles
          jetIndices1.push_back(i);
          jetIndices2.push_back(j);
          jetAngles.push_back(angleInDegrees);
        }
      }

      //Apply cut that the smallest angles must be less than 100 degrees
      if(minJetAngle1 > 100) vetoEvent;
      if(minJetAngle2 > 100) vetoEvent;

      //Cut if the pairs of jets are too wide
      _h_cutflow->fill(3);

      // These two angles don't share a jet
      // Not sure if we should veto the event, or just find another pairing
      if( minJetPair_1_jet1index ==  minJetPair_2_jet1index || minJetPair_1_jet1index ==  minJetPair_2_jet2index ||
          minJetPair_1_jet2index ==  minJetPair_2_jet1index || minJetPair_1_jet2index ==  minJetPair_2_jet2index){
         vetoEvent;
      }

      //Cut if jets are not distinct pairs
      _h_cutflow->fill(4);

      bool isPass = false;
      int index1 = -1;
      int index2 = -1;
      int index3 = -1;
      int index4 = -1;

      //Angles between W jets should be between 100 and 140 degrees
      for(unsigned int i=0; i<jetAngles.size(); i++){
        if(jetAngles[i] < 100 || jetAngles[i] > 140) continue;
        for(unsigned int j=i+1; j<jetAngles.size(); j++){
          if(jetAngles[j] < 100 || jetAngles[j] > 140) continue;

          //Loop through all pairings until we find a non adjacent one
          if( jetIndices1[i] ==  jetIndices1[j] || jetIndices2[i] ==  jetIndices2[j] ||
              jetIndices1[i] ==  jetIndices2[j] || jetIndices2[i] ==  jetIndices1[j]) continue;

          isPass = true;
          //Indexes 1 and 2 should be the first pair, indexes 3 and 4 should be the second
          index1 = jetIndices1[i];
          index2 = jetIndices2[i];
          index3 = jetIndices1[j];
          index4 = jetIndices2[j];
        }
      }

      if(!isPass) vetoEvent;

      //cut if the jets don't fit expected geometry
      _h_cutflow->fill(5);

      //We're gonna try to ensure the jets are listed in clockwise order
      //The paper says that the interjet regions have the smallest angles, which we found.
      
      std::array<int,4> indices = {index1,index2,index3,index4};
      std::vector<int> ordered_indices;

      //Start by putting the first of the smallest angle region
      ordered_indices.push_back(minJetPair_1_jet1index);

      //loop until we find which index matches it, then add its pair
      for(unsigned int i=0; i<indices.size();i++){
        if(minJetPair_1_jet1index == indices[i]){
          if(i%2 == 0){
            ordered_indices.push_back(indices[i+1]);
            break;
          }
          else{
            ordered_indices.push_back(indices[i-1]);
            break;
          }
        }
      }

      //Add the next jet based on the smallest jet pairs
      if(ordered_indices[1] == minJetPair_2_jet1index)
        ordered_indices.push_back(minJetPair_2_jet2index);
      else
        ordered_indices.push_back(minJetPair_2_jet1index);
      
      //Add the last remaining jet, and the first jet again for looping
      ordered_indices.push_back(minJetPair_1_jet2index);
      ordered_indices.push_back(minJetPair_1_jet1index);
      
      //Ok now ordered_indices should contain all jets in order clockwise
      //Note: First two are from the same W, next 2 are the other W

      fastjet::PseudoJet correctJetPair1 = selectedJets[ordered_indices[0]] + selectedJets[ordered_indices[1]];
      fastjet::PseudoJet correctJetPair2 = selectedJets[ordered_indices[2]] + selectedJets[ordered_indices[3]];

      _h_dijetMass1->fill(correctJetPair1.m());
      _h_dijetMass2->fill(correctJetPair2.m());
      _h_mult_after->fill(particles.size());

      //For checking proper ordering
      check_index->fill(ordered_indices[0]);
      check_index->fill(ordered_indices[1]);
      check_index->fill(ordered_indices[2]);
      check_index->fill(ordered_indices[3]);


      //Parton Comparison Graph/ Truth particle analysis from Iza

      const HepMC3::GenEvent &ge = *event.genEvent();
      //  Find hard-scatter W bosons
      HepMC3::ConstGenParticlePtr Wp = nullptr;
      HepMC3::ConstGenParticlePtr Wm = nullptr;

      for (const auto &p : ge.particles()) {
        if (abs(p->pid()) != 24)
          continue;
        if (!p->production_vertex())
          continue;
        bool prod_ee = true;
        for (const auto &parent : p->production_vertex()->particles_in()) {
          if (abs(parent->pid()) != 11) {
            prod_ee = false;
            break;
          }
        }
        if (!prod_ee)
          continue;

        if (p->pid() == 24)
          Wp = p;
        if (p->pid() == -24)
          Wm = p;
      }
      if (!Wp || !Wm)
        return;

      //  Get quarks from W decays
      std::vector<Particle> WpQuarks, WmQuarks;

      if (Wp->end_vertex()) {
        for (const auto &q : Wp->end_vertex()->particles_out())
          if (PID::isQuark(q->pid()))
            WpQuarks.push_back(Particle(q));
      }
      if (Wm->end_vertex()) {
        for (const auto &q : Wm->end_vertex()->particles_out())
          if (PID::isQuark(q->pid()))
            WmQuarks.push_back(Particle(q));
      }
      if (WpQuarks.size() != 2 || WmQuarks.size() != 2)
        return;

      //Iza's ordering scheme
      //  Order Wp quarks by pT
      if (WpQuarks[1].pT() > WpQuarks[0].pT())
        std::swap(WpQuarks[0], WpQuarks[1]);

      Particle q1 = WpQuarks[0];
      Particle q2 = WpQuarks[1];

      // Pick q3 based on min dR to q1
      double dR0 = deltaR(q1.momentum(), WmQuarks[0].momentum());
      double dR1 = deltaR(q1.momentum(), WmQuarks[1].momentum());

      Particle q3 = (dR0 < dR1 ? WmQuarks[0] : WmQuarks[1]);
      Particle q4 = (dR0 < dR1 ? WmQuarks[1] : WmQuarks[0]);

      //Checking invariant mass of quarks

      double massWp = (q1.momentum() + q2.momentum()).mass();
      double massWm = (q3.momentum() + q4.momentum()).mass();

      quarkMassWp -> fill(massWp);
      quarkMassWm -> fill(massWm);

      //Added to match region graphs

      // Setup 4-momentum
      FourMomentum fv_q1(q1.E(), q1.px(), q1.py(), q1.pz());
      FourMomentum fv_q2(q2.E(), q2.px(), q2.py(), q2.pz());
      FourMomentum fv_q3(q3.E(), q3.px(), q3.py(), q3.pz());
      FourMomentum fv_q4(q4.E(), q4.px(), q4.py(), q4.pz());
      FourMomentum fv_Wp(Wp->momentum().e(), Wp->momentum().px(), Wp->momentum().py(), Wp->momentum().pz());
      FourMomentum fv_Wm(Wm->momentum().e(), Wm->momentum().px(), Wm->momentum().py(), Wm->momentum().pz());

      std::vector<FourMomentum> momentaList = {fv_q1, fv_q2, fv_q3, fv_q4, fv_q1};
      double angleWpDeg = fv_q1.angle(fv_q2) * 180.0 / M_PI;
      double angleWmDeg = fv_q3.angle(fv_q4) * 180.0 / M_PI;

      if (angleWpDeg > 100.0 && angleWpDeg < 140.0 && angleWmDeg > 100.0 && angleWmDeg < 140.0) {
        for (int i = 0; i < 4; ++i) {
          double deg_theta = momentaList[i].angle(momentaList[i+1]) * 180/3.14;
          if(i==0)
            regionQA->fill(deg_theta);
          if(i==1)
            regionQB->fill(deg_theta);
          if(i==2)
            regionQC->fill(deg_theta);
          if(i==3)
            regionQD->fill(deg_theta);
        }
      }

      //Now we need 4 histograms, one for each quark, comparing its deltaR to the jets
      
      std::vector<Histo1DPtr> histos = {quark1,quark2,quark3,quark4};
      std::vector<Histo1DPtr> histosRatio = {jetToQuark1, jetToQuark2, jetToQuark3, jetToQuark4};
      std::vector<Particle> truth_quarks = {q1,q2,q3,q4,q1};
      std::vector<Histo1DPtr> quarkdR = {quarkdR1, quarkdR2, quarkdR3, quarkdR4};

      for(int i = 0; i <4; ++i){
        double deltaR2 = 1000;
        double matchedJetindex = -1;

        //Finds smallest deltaR between a truth quark and any jet
        for(int j = 0; j < 4; ++j){
          FourMomentum loop = FourVector(selectedJets[indices[j]].E(), selectedJets[indices[j]].px() , selectedJets[indices[j]].py(), selectedJets[indices[j]].pz());
          double deltaR1 = deltaR(truth_quarks[i].momentum(),loop);
          if(deltaR1 < deltaR2){
            deltaR2 = deltaR1;

            //identify matchedJetindex with the jet paired to selectedJets[j]
            if (j%2 == 0)
              matchedJetindex = indices[j+1];
            else
              matchedJetindex = indices[j-1];
          }
        }

        histos[i]->fill(deltaR2);

        //for the first quark, we found the jet cooresponding to minimum delta R
        //We then graph the deltaR of the other quark and jet
        if (i==0){
          FourMomentum matchedJet = FourVector(selectedJets[matchedJetindex].E(), selectedJets[matchedJetindex].px() , selectedJets[matchedJetindex].py(), selectedJets[matchedJetindex].pz());
          mismatchDeltaR -> fill(deltaR(matchedJet,truth_quarks[1].momentum()));
        }

        //Checking ratio of quark pT to jet pT
        for(int j = 0; j < 4; ++j){
          FourMomentum loop = FourVector(selectedJets[j].E(), selectedJets[j].px() , selectedJets[j].py(), selectedJets[j].pz());
          histosRatio[i]->fill(truth_quarks[i].pT()/loop.pT());
        }
        //Check dR between quarks
        quarkdR[i]->fill(deltaR(truth_quarks[i],truth_quarks[i+1]));
      }

      //
      (*out) << "Event" << std::endl;

      // Wp
      (*out) << correctJetPair1.theta() << " " << correctJetPair1.phi() << " " << correctJetPair1.E() << "\n";
      (*out) << fv_Wp.theta() << " " << fv_Wp.phi() << " " << fv_Wp.E() << "\n";
      (*out) << fv_q1.theta() << " " << fv_q1.phi() << " " << fv_q1.E() << "\n";
      (*out) << fv_q2.theta() << " " << fv_q2.phi() << " " << fv_q2.E() << "\n";

      // Wm
      (*out) << correctJetPair2.theta() << " " << correctJetPair2.phi() << " " << correctJetPair2.E() << "\n";
      (*out) << fv_Wm.theta() << " " << fv_Wm.phi() << " " << fv_Wm.E() << "\n";
      (*out) << fv_q3.theta() << " " << fv_q3.phi() << " " << fv_q3.E() << "\n";
      (*out) << fv_q4.theta() << " " << fv_q4.phi() << " " << fv_q4.E() << "\n";

      //invariant mass
      (*out) << massWp << " " << massWm << " " << correctJetPair1.m() << " " << correctJetPair2.m() << "\n";


      (*out) << "Jet1" << std::endl;
      (*out) << "Constit1" << std::endl;
      for(unsigned int i=0; i<selectedJets[ordered_indices[0]].constituents().size(); i++){
        (*out) << selectedJets[ordered_indices[0]].constituents()[i].theta() << " " << selectedJets[ordered_indices[0]].constituents()[i].phi() << " " << selectedJets[ordered_indices[0]].constituents()[i].e() << "\n";
      }

      (*out) << "Jet2" << std::endl;
      (*out) << "Constit2" << std::endl;
      for(unsigned int i=0; i<selectedJets[ordered_indices[1]].constituents().size(); i++){
        (*out) << selectedJets[ordered_indices[1]].constituents()[i].theta() << " " << selectedJets[ordered_indices[1]].constituents()[i].phi() << " " << selectedJets[ordered_indices[1]].constituents()[i].e() << "\n";
      }

      (*out) << "Jet3" << std::endl;
      (*out) << "Constit3" << std::endl;
      for(unsigned int i=0; i<selectedJets[ordered_indices[2]].constituents().size(); i++){
        (*out) << selectedJets[ordered_indices[2]].constituents()[i].theta() << " " << selectedJets[ordered_indices[2]].constituents()[i].phi() << " " << selectedJets[ordered_indices[2]].constituents()[i].e() << "\n";
      }

      (*out) << "Jet4" << std::endl;
      (*out) << "Constit4" << std::endl;
      for(unsigned int i=0; i<selectedJets[ordered_indices[3]].constituents().size(); i++){
        (*out) << selectedJets[ordered_indices[3]].constituents()[i].theta() << " " << selectedJets[ordered_indices[3]].constituents()[i].phi() << " " << selectedJets[ordered_indices[3]].constituents()[i].e() << "\n";
      }

      //A quest for Figure 6

      double angleA;
      double angleB;
      double angleC;
      double angleD;

      //So we see how much we are double counting or not counting particles in the sort
      std::vector<int> particleSort(particles.size(), 0);

      //Calculating thetaRescaled for each particle and ensuring each is placed in 1 or less regions
      for (unsigned int j = 0; j < particles.size(); j++){
        std::vector<double> sorter(4,0);
        std::vector<double> rescaledAngle(4,0);
        Vector3 p(particles[j].px(), particles[j].py(), particles[j].pz());

        for (int i = 0; i < 4; ++i) {
          //First we will define a normal vector from a couple of 3 Vectors:
          Vector3 a(selectedJets[ordered_indices[i]].px(), selectedJets[ordered_indices[i]].py(), selectedJets[ordered_indices[i]].pz());
          Vector3 b(selectedJets[ordered_indices[i+1]].px(), selectedJets[ordered_indices[i+1]].py(), selectedJets[ordered_indices[i+1]].pz());

          Vector3 n = a.cross(b);
          double nmag = n.mod();
          Vector3 unit_n = n / nmag;

          //Calculate refrence jet angle
          double thetaRef = atan2(n.dot(unit_n), a.dot(b));
          double deg_theta = a.angle(b) * 180/3.14;
          Vector3 proj_p = p - unit_n * p.dot(unit_n);

          if(j==0){
            if(i==0)
              angleA = deg_theta;
            if(i==1)
              angleB = deg_theta;
            if(i==2)
              angleC = deg_theta;
            if(i==3)
              angleD = deg_theta;
          }

          //Calculate angles with a as the reference jet
          double theta_p = atan2(a.cross(proj_p).dot(unit_n), a.dot(proj_p));
          //Find transverse momentum
          double pT_to_plane = std::abs(p.dot(unit_n));

          //if its in the region, have sorter write it down
          if ((0 < theta_p && theta_p < thetaRef) || (0 > theta_p && theta_p > thetaRef)){
            sorter[i] = pT_to_plane;
            rescaledAngle[i] = (theta_p/thetaRef);
          }
        }

        double minVal = 1000;
        int jetNum = -1;
        for(int k=0; k<4; ++k){
          if(sorter[k]!=0 && minVal>sorter[k]){
            minVal = sorter[k];
            jetNum = k;
          }
        }

        if(jetNum != -1){
          particleSort[j]++;
          phi_rescaled->fill(rescaledAngle[jetNum]+jetNum);

          //Check if inside or outside region and fill cooresponding histo
          if (jetNum == 0 || jetNum == 2){
            insideTotal->fill(rescaledAngle[jetNum]);
          }
          if (jetNum == 1 || jetNum == 3){
            OutsideTotal->fill(rescaledAngle[jetNum]);
          }
          //regionA
          if (jetNum==0)
            inside1->fill(rescaledAngle[jetNum]);
          //RegionB
          if (jetNum==1)
            outside1->fill(rescaledAngle[jetNum]);
          //RegionC
          if (jetNum==2)
            inside2->fill(rescaledAngle[jetNum]);
          //RegionD
          if (jetNum==3)
            outside2->fill(rescaledAngle[jetNum]);
        }
      }
      //Angle in degrees between each jet
      regionA->fill(angleA);
      regionB->fill(angleB);
      regionC->fill(angleC);
      regionD->fill(angleD);

      for (int i = 0; i<particles.size(); i++){
        sorted->fill(particleSort[i]);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      normalize(_h_dijetMass1);
      normalize(_h_dijetMass2);
      normalize(_h_ycutAll);
      normalize(_h_selectedJets);
      normalize(_h_mult_before);
      normalize(_h_mult_after);
      scale(phi_rescaled, 1.0 / sumOfWeights());
      scale(insideTotal, 1.0 / sumOfWeights());
      scale(OutsideTotal, 1.0 / sumOfWeights());

      //Find actual ratio for figure 6

      for (size_t i = 0; i < insideTotal->numBins(); ++i) {
          
          // Get the bin content
          double insideValue  = insideTotal->bin(i).sumW();
          double outsideValue = OutsideTotal->bin(i).sumW();

          // Compute the ratio safely
          double ratio = (outsideValue != 0) ? insideValue / outsideValue : 0;

          // Fill the ratio into the corresponding bin of region_ratio
          region_ratio->fillBin(i, ratio);
      } 

      for (unsigned int i = 0; i < inside1->numBins(); ++i) {
        double inside1Value  = inside1->bin(i).sumW();
        double inside2Value = inside2->bin(i).sumW();
        double outside1Value = outside1->bin(i).sumW();
        double outside2Value = outside2->bin(i).sumW();

        // Compute the ratio safely
        double in_ratio = (inside2Value != 0) ? inside1Value / inside2Value : 0;
        double out_ratio = (outside2Value != 0) ? outside1Value / outside2Value : 0;

        // Fill the ratio into the corresponding bin of region_ratio
        inside_ratio->fillBin(i, in_ratio);
        outside_ratio->fillBin(i, out_ratio);
      }
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

    Histo1DPtr quarkMassWm;
    Histo1DPtr quarkMassWp;

    //Same thing, std::vector<Histo1DPtr> _h_cutflow;
    Histo1DPtr _h_cutflow;
    //ycut before event cuts
    Histo1DPtr _h_ycutAll;
    //# of selected jets before cut
    Histo1DPtr _h_selectedJets;
    //particle multiplicity before cuts
    Histo1DPtr _h_mult_before;
    //particle multiplicity after cuts
    Histo1DPtr _h_mult_after;

    //Figure 6, particle flow per event in each region
    Histo1DPtr phi_rescaled;

    //Figure 5, ratio between particle count inside/outside
    Histo1DPtr insideTotal;
    Histo1DPtr OutsideTotal;
    Histo1DPtr region_ratio;

    //Diagnostic plots
    Histo1DPtr _h_mult_nocut;
    Histo1DPtr outside_ratio;
    Histo1DPtr inside_ratio;
    Histo1DPtr inside1;
    Histo1DPtr inside2;
    Histo1DPtr outside1;
    Histo1DPtr outside2;

    Histo1DPtr regionA;
    Histo1DPtr regionB;
    Histo1DPtr regionC;
    Histo1DPtr regionD;

    Histo1DPtr quark1;
    Histo1DPtr quark2;
    Histo1DPtr quark3;
    Histo1DPtr quark4;

    Histo1DPtr mismatchDeltaR;

    Histo1DPtr check_index;

    Histo1DPtr jetToQuark1;
    Histo1DPtr jetToQuark2;
    Histo1DPtr jetToQuark3;
    Histo1DPtr jetToQuark4;

    Histo1DPtr quarkdR1;
    Histo1DPtr quarkdR2;
    Histo1DPtr quarkdR3;
    Histo1DPtr quarkdR4;

    Histo1DPtr regionQA;
    Histo1DPtr regionQB;
    Histo1DPtr regionQC;
    Histo1DPtr regionQD;

    Histo1DPtr region_a_thetaRescaled;
    Histo1DPtr region_b_thetaRescaled;
    Histo1DPtr region_c_thetaRescaled;
    Histo1DPtr region_d_thetaRescaled;

    Histo1DPtr sorted;

    std::ofstream* out;
  };


  RIVET_DECLARE_PLUGIN(MC_COLORRECONNECTION);

}


