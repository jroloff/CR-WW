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


    double getAngleInDegrees(double radians){
      return radians*180/3.14159;
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

      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Book histograms

      // One histogram for each of the two jet pairs
      book(_h_dijetMass1, "dijet_mass_1",  40, 0, 120);
      book(_h_dijetMass2, "dijet_mass_2",  40, 0, 120);
      book(_h_dijetMass1_allMatched, "dijet_mass_1_allMatched",  40, 0, 120);
      book(_h_dijetMass2_allMatched, "dijet_mass_2_allMatched",  40, 0, 120);

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
      book(insideTotal, "Inside_W_Region", 20, 0, 1.);
      book(outsideTotal, "Outside_W_Region", 20, 0, 1.);
      book(region_ratio, "Region_Ratio");

      book(insideTotal_allMatched, "Inside_W_Region_allMatched", 20, 0, 1.);
      book(outsideTotal_allMatched, "Outside_W_Region_allMatched", 20, 0, 1.);
      book(region_ratio_allMatched, "Region_Ratio_allMatched");

      //Diagnostic Plots
      book(_h_mult_nocut, "Multiplicity_before_cuts", 151, -0.5, 150.5);

      //ThetaRescaled for each region
      book(inside1, "inside1", 24, 0.0, 1.0);
      book(inside2, "inside2", 24, 0.0, 1.0);
      book(outside1, "outside1", 24, 0.0, 1.0);
      book(outside2, "outside2", 24, 0.0, 1.0);

      book(inside1_allMatched, "inside1_allMatched", 24, 0.0, 1.0);
      book(inside2_allMatched, "inside2_allMatched", 24, 0.0, 1.0);
      book(outside1_allMatched, "outside1_allMatched", 24, 0.0, 1.0);
      book(outside2_allMatched, "outside2_allMatched", 24, 0.0, 1.0);

      book(outside_ratio, "Outside_Ratio");
      book(inside_ratio, "Inside_Ratio");

      book(outside_ratio_allMatched, "Outside_Ratio_allMatched");
      book(inside_ratio_allMatched, "Inside_Ratio_allMatched");

      book(orderingChoice, "orderingChoice", 4, -0.5, 3.5);
      book(orderingChoice_allMatched, "orderingChoice_allMatched", 4, -0.5, 3.5);
      book(isMatched, "isMatched", 2, -0.5, 1.5);

      //Angle between each region
      book(regionA, "Region_A", 72, -0.5, 180.5);
      book(regionB, "Region_B", 72, -0.5, 180.5);
      book(regionC, "Region_C", 72, -0.5, 180.5);
      book(regionD, "Region_D", 72, -0.5, 180.5);
      regions = {regionA, regionB, regionC, regionD};

      book(quark1, "Quark1_deltaR", 50, 0, 1.5);
      book(quark2, "Quark2_deltaR", 50, 0, 1.5);
      book(quark3, "Quark3_deltaR", 50, 0, 1.5);
      book(quark4, "Quark4_deltaR", 50, 0, 1.5);
      histos = {quark1,quark2,quark3,quark4};

      book(mismatchDeltaR, "Matching_DeltaR", 50, 0, 3);

      book(check_index, "check_index", 6, -1, 5);

      book(jetToQuark1, "pT_Ratio_1", 100, 0, 3);
      book(jetToQuark2, "pT_Ratio_2", 100, 0, 3);
      book(jetToQuark3, "pT_Ratio_3", 100, 0, 3);
      book(jetToQuark4, "pT_Ratio_4", 100, 0, 3);
      histosRatio = {jetToQuark1, jetToQuark2, jetToQuark3, jetToQuark4};

      book(quarkdR1, "quark_region_dR1", 50, 0, 1.5);
      book(quarkdR2, "quark_region_dR2", 50, 0, 1.5);
      book(quarkdR3, "quark_region_dR3", 50, 0, 1.5);
      book(quarkdR4, "quark_region_dR4", 50, 0, 1.5);
      quarkdR = {quarkdR1, quarkdR2, quarkdR3, quarkdR4};

      //Using truth quarks
      book(regionQA, "Quark_Region_A", 72, -0.5, 180.5);
      book(regionQB, "Quark_Region_B", 72, -0.5, 180.5);
      book(regionQC, "Quark_Region_C", 72, -0.5, 180.5);
      book(regionQD, "Quark_Region_D", 72, -0.5, 180.5);
      regionsQ = {regionQA, regionQB, regionQC, regionQD};

      book(sorted, "Sorted_particles", 13, -0.5, 12.5);

      out = new std::ofstream("outEventDisplay.csv");

    }

    bool tryMatch(int u,
              const std::vector<std::vector<int>>& adj,
              std::vector<int>& matchTo,
              std::vector<bool>& visited) {
      for (int v : adj[u]) {
        if (visited[v]) continue;
        visited[v] = true;

        // If v is free OR we can re-match its partner
        if (matchTo[v] == -1 || tryMatch(matchTo[v], adj, matchTo, visited)) {
          matchTo[v] = u;
          return true;
        }
      }
      return false;
    }

    void writeEvent(fastjet::PseudoJet correctJetPair1, fastjet::PseudoJet correctJetPair2, FourMomentum fv_Wp, FourMomentum fv_Wm, FourMomentum fv_q1, FourMomentum fv_q2, FourMomentum fv_q3, FourMomentum fv_q4,
                    double massWp, double massWm, std::vector<fastjet::PseudoJet> selectedJets, std::vector<int> ordered_indices){
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
      // Just so we don't have to keep doing this conversion over and over
      std::vector<FourMomentum> jet4mom;

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
            FourMomentum fv1 = FourVector(selectedJets[i].E(), selectedJets[i].px() , selectedJets[i].py(), selectedJets[i].pz());
            jet4mom.push_back(fv1);
          }
        }
      }

      //Histo of selected jets
      _h_selectedJets->fill(durjet.size());

      // Need exactly 4 jets
      if (selectedJets.size() != 4) vetoEvent;

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
      if (!Wp || !Wm){
        std::cout << "Did not find both Ws" << std::endl;
        return;
      }

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
      if (WpQuarks.size() != 2 || WmQuarks.size() != 2){
        //std::cout << "Did not find enough quarks" << std::endl;
        //std::cout << WpQuarks.size() << "\t" << WmQuarks.size() << std::endl;
        return;
      }

      //Add event that survived 4-jet selection
      _h_cutflow->fill(2);

      // Two smallest angles < 100 degrees
      int minJetPair_1_jet1index = -1;
      int minJetPair_1_jet2index = -1;

      int minJetPair_2_jet1index = -1;
      int minJetPair_2_jet2index = -1;

      double minJetAngle1 = 1000;
      double minJetAngle2 = 1000;

      for(unsigned int i=0; i<jet4mom.size(); i++){
        for(unsigned int j=i+1; j<jet4mom.size(); j++){
          if(i==j) continue;
          double angleInDegrees = getAngleInDegrees(jet4mom[i].angle(jet4mom[j]));

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

      //Apply cut that the smallest angles must be less than 100 degrees
      if(minJetAngle1 > 100) vetoEvent;
      if(minJetAngle2 > 100) vetoEvent;


      //Cut if the pairs of jets are too wide
      _h_cutflow->fill(3);

      // These two angles don't share a jet
      // Not sure if we should veto the event, or just find another pairing,
      // but we are assuming we veto the event
      if( minJetPair_1_jet1index ==  minJetPair_2_jet1index || minJetPair_1_jet1index ==  minJetPair_2_jet2index ||
          minJetPair_1_jet2index ==  minJetPair_2_jet1index || minJetPair_1_jet2index ==  minJetPair_2_jet2index){
         vetoEvent;
      }

      //Cut if jets are not distinct pairs
      _h_cutflow->fill(4);


      
      //Parton Comparison Graph/ Truth particle analysis from Iza
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
      std::vector<Particle> truth_quarks = {q1,q2,q3,q4};

      //Checking invariant mass of quarks
      double massWp = (q1.momentum() + q2.momentum()).mass();
      double massWm = (q3.momentum() + q4.momentum()).mass();

      quarkMassWp -> fill(massWp);
      quarkMassWm -> fill(massWm);

      // Note that this region matching is different than what is implemented for reco level
      // I'm assuming part of this motivation is because we know which W is which?
      //Added to match region graphs
      // Setup 4-momentum
      FourMomentum fv_q1(q1.E(), q1.px(), q1.py(), q1.pz());
      FourMomentum fv_q2(q2.E(), q2.px(), q2.py(), q2.pz());
      FourMomentum fv_q3(q3.E(), q3.px(), q3.py(), q3.pz());
      FourMomentum fv_q4(q4.E(), q4.px(), q4.py(), q4.pz());
      FourMomentum fv_Wp(Wp->momentum().e(), Wp->momentum().px(), Wp->momentum().py(), Wp->momentum().pz());
      FourMomentum fv_Wm(Wm->momentum().e(), Wm->momentum().px(), Wm->momentum().py(), Wm->momentum().pz());

      std::vector<FourMomentum> momentaList = {fv_q1, fv_q2, fv_q3, fv_q4};
      double angleWpDeg = getAngleInDegrees(fv_q1.angle(fv_q2));
      double angleWmDeg = getAngleInDegrees(fv_q3.angle(fv_q4));

      if (angleWpDeg > 100.0 && angleWpDeg < 140.0 && angleWmDeg > 100.0 && angleWmDeg < 140.0) {
        for (int i = 0; i < 4; ++i) {
          double deg_theta = getAngleInDegrees(momentaList[i].angle(momentaList[(i+1)%(truth_quarks.size())]));
          regionsQ[i]->fill(deg_theta);
        }
      }

      //Now we need 4 histograms, one for each quark, comparing its deltaR to the jets
      for(int i = 0; i <truth_quarks.size(); ++i){
        double minDR = 1000;
        double matchedJetIndex = -1;

        //Finds smallest deltaR between a truth quark and any jet
        for(int j = 0; j < selectedJets.size(); j++){
          double deltaR1 = deltaR(truth_quarks[i].momentum(),jet4mom[j]);
          if(deltaR1 < minDR){
            minDR = deltaR1;
            matchedJetIndex = j;
          }
        }

        histos[i]->fill(minDR);

        //for the first quark, we found the jet cooresponding to minimum delta R
        //We then graph the deltaR of the other quark and jet
        // JKR: I might have messed this up, so sorry about that. I was quite confused about what the code was doing...
        for(int j = 0; j < selectedJets.size(); ++j){
          if(j==matchedJetIndex) continue;
          mismatchDeltaR -> fill(deltaR(jet4mom[matchedJetIndex],truth_quarks[i].momentum()));
        }

        //Checking ratio of quark pT to jet pT
        // TODO shouldn't we just do this for the case where it is matched?
        histosRatio[i]->fill(truth_quarks[i].pT()/jet4mom[matchedJetIndex].pT());

        //Check dR between quarks
        quarkdR[i]->fill(deltaR(truth_quarks[i],truth_quarks[(i+1)%truth_quarks.size()]));
      }

      bool allJetsMatched = true;

      // This comes from chatGPT, so I'm not 100% sure it's all correct
      std::vector<std::vector<int>> adj(4);

      double maxDeltaR = 0.1;
      // Build all pairwise distances
      for (size_t j = 0; j < selectedJets.size(); j++) {
        for (size_t i = 0; i < truth_quarks.size(); i++) {
          double dR = deltaR(truth_quarks[i].momentum(), jet4mom[j]); 
          if (dR < maxDeltaR && std::abs(truth_quarks[i].pT()/jet4mom[j].pT()-1)<0.2) {
            adj[i].push_back(j);
          }
        }
      }

      int N = 4;
      std::vector<int> matchTo(N, -1); // which v1 index each v2 node is matched to

      int matchCount = 0;
      for (int u = 0; u < N; ++u) {
        std::vector<bool> visited(N, false);
        if (tryMatch(u, adj, matchTo, visited)) {
          matchCount++;
        }
      }


      // Two options for non-adjacent jet pairings, once we have selected the jets that make regions B and D
      // j1_1 j2_1, j1_2 j2_2
      double angle1_1 = getAngleInDegrees(jet4mom[minJetPair_1_jet1index].angle(jet4mom[minJetPair_2_jet1index]));
      double angle1_2 = getAngleInDegrees(jet4mom[minJetPair_1_jet2index].angle(jet4mom[minJetPair_2_jet2index]));
      double sumAngles1 = (angle1_1 > 100 && angle1_1 < 140 && angle1_2 > 100 && angle1_2 < 140)? angle1_1 + angle1_2:0;

      // j1_1 j2_2, j1_2 j2_1
      double angle2_1 = getAngleInDegrees(jet4mom[minJetPair_1_jet1index].angle(jet4mom[minJetPair_2_jet2index]));
      double angle2_2 = getAngleInDegrees(jet4mom[minJetPair_1_jet2index].angle(jet4mom[minJetPair_2_jet1index]));
      double sumAngles2 = (angle2_1 > 100 && angle2_1 < 140 && angle2_2 > 100 && angle2_2 < 140)?angle2_1 + angle2_2:0;

      // One of the pairings must pass the selection
      if(sumAngles1==0 && sumAngles2==0) vetoEvent;
      //cut if the jets don't fit expected geometry
      _h_cutflow->fill(5);

      int index1, index2, index3, index4;
      // Use the pairing based on the larger sum of angles
      // Jets 1-2 form region A (largest angle)
      // Jets 2-3 form region B (smallest angle, which must be from minJetPair_1 by definition)
      // Jets 3-4 form region C (second-largest angle)
      // Jets 4-1 form region D (second-smallest angle, which must be from minJetPair_2 by definition)
      //std::cout << "\t\t\t\t\t\t " << angle2_1 << "\t" << angle2_2 <<"\t" << angle1_1 << "\t" << angle1_2 << "\t" << minJetAngle1  << "\t" << minJetAngle2 << std::endl;
      if(sumAngles1 > sumAngles2){
        if(angle1_2 > angle1_1){
          index2 = minJetPair_1_jet2index;
          index3 = minJetPair_1_jet1index;
          index4 = minJetPair_2_jet1index;
          index1 = minJetPair_2_jet2index;
          orderingChoice->fill(0);
          //if(allJetsMatched){
          //  orderingChoice_allMatched->fill(0);
          //}
        }
        else{
          index2 = minJetPair_1_jet1index;
          index3 = minJetPair_1_jet2index;
          index4 = minJetPair_2_jet2index;
          index1 = minJetPair_2_jet1index;
          orderingChoice->fill(1);
          //if(allJetsMatched){
          //  orderingChoice_allMatched->fill(1);
          //}
        }
      }
      else{
        // Both of these versions look asymmetric, but in opposite ways,
        // so I think they basically cancel each other out.
        // I have not been able to figure out why though.

        if(angle2_2 > angle2_1){
          index2 = minJetPair_1_jet2index;
          index3 = minJetPair_1_jet1index;
          index4 = minJetPair_2_jet2index;
          index1 = minJetPair_2_jet1index;
          orderingChoice->fill(2);
          //if(allJetsMatched){
          //  orderingChoice_allMatched->fill(2);
          //}
        }
        else{
          //std::cout << "D" << std::endl;
          //std::cout << "\t\t\t\t\t\t " << angle2_1 << "\t" << angle2_2 <<"\t" << angle1_1 << "\t" << angle1_2 << "\t" << minJetAngle1  << "\t" << minJetAngle2 << std::endl;
          index2 = minJetPair_1_jet1index;
          index3 = minJetPair_1_jet2index;
          index4 = minJetPair_2_jet1index;
          index1 = minJetPair_2_jet2index;
          orderingChoice->fill(3);
          //if(allJetsMatched){
          //  orderingChoice_allMatched->fill(3);
          //}
        }
      }
      // We have passed the full cutflow now

      // We have already ordered them, so we don't need to do this again
      std::vector<int> ordered_indices = {index1,index2,index3,index4};
      
      for(unsigned int i=0; i<ordered_indices.size(); i++){
        check_index->fill(ordered_indices[i]);
      }


      //Note: First two are from the same W, next 2 are the other W
      fastjet::PseudoJet correctJetPair1 = selectedJets[index1] + selectedJets[index2];
      fastjet::PseudoJet correctJetPair2 = selectedJets[index3] + selectedJets[index4];
      if(matchCount !=4) allJetsMatched = false;
      // Check if the inside regions correspond to a W, or to mismatched W's
      // We only want the correct matching
      if(allJetsMatched){
          if((matchTo[ordered_indices[0]]==0 && matchTo[ordered_indices[1]] !=1 )|| (matchTo[ordered_indices[0]]==2 && matchTo[ordered_indices[1]] !=3 )) allJetsMatched = false;
          if((matchTo[ordered_indices[1]]==0 && matchTo[ordered_indices[0]] !=1 )|| (matchTo[ordered_indices[1]]==2 && matchTo[ordered_indices[0]] !=3 )) allJetsMatched = false;
      }

      if(sumAngles1 > sumAngles2){
        if(angle1_2 > angle1_1){
          if(allJetsMatched){
            orderingChoice_allMatched->fill(0);
          }
        }
        else{
          if(allJetsMatched){
            orderingChoice_allMatched->fill(1);
          }
        }
      }
      else{
        if(angle2_2 > angle2_1){
          if(allJetsMatched){
            orderingChoice_allMatched->fill(2);
          }
        }
        else{
          if(allJetsMatched){
            orderingChoice_allMatched->fill(3);
          }
        }
      }



      _h_dijetMass1->fill(correctJetPair1.m());
      _h_dijetMass2->fill(correctJetPair2.m());
      if(allJetsMatched){
        _h_dijetMass1_allMatched->fill(correctJetPair1.m());
        _h_dijetMass2_allMatched->fill(correctJetPair2.m());
      }
      _h_mult_after->fill(particles.size());


      isMatched->fill(allJetsMatched);

      // Pre-calculating stuff so we don't have to keep doing this every time
      std::vector<Vector3> a_vec;
      std::vector<Vector3> unit_n_vec;
      std::vector<double> thetaRef_vec;
      for (int i = 0; i < jet4mom.size(); i++) {
          Vector3 a(selectedJets[ordered_indices[i]].px(), selectedJets[ordered_indices[i]].py(), selectedJets[ordered_indices[i]].pz());
          // Need to make sure it wraps around cleanly
          Vector3 b(selectedJets[ordered_indices[(i+1)%jet4mom.size()]].px(), selectedJets[ordered_indices[(i+1)%jet4mom.size()]].py(), selectedJets[ordered_indices[(i+1)%jet4mom.size()]].pz());

          //First we will define a normal vector from a couple of 3-Vectors:
          Vector3 n = a.cross(b);
          double nmag = n.mod();
          Vector3 unit_n = n / nmag;

          //Calculate reference jet angle
          double thetaRef = atan2(n.dot(unit_n), a.dot(b));
          //double thetaRef = a.angle(b);
          //std::cout << thetaRef << "\t" << a.angle(b) << std::endl;

          a_vec.push_back(a);
          unit_n_vec.push_back(unit_n); 
          thetaRef_vec.push_back(thetaRef);

          double deg_theta = getAngleInDegrees(a.angle(b));
          regions[i]->fill(deg_theta);
      }  

      //A quest for Figure 6
      //So we see how much we are double counting or not counting particles in the sort
      std::vector<int> particleSort(particles.size(), 0);

      //Calculating thetaRescaled for each particle and ensuring each is placed in 1 or less regions
      for (unsigned int j = 0; j < particles.size(); j++){
        Vector3 p(particles[j].px(), particles[j].py(), particles[j].pz());

        // sorter has the pt_to_plane for a given particle
        std::vector<double> sorter(4,0);
        // The rescaled angle for each region for this particle
        std::vector<double> rescaledAngle(4,0);

        for (int i = 0; i < unit_n_vec.size(); i++) {
          //Calculate refrence jet angle by subtracting off the component of the vector in the direction of the unit vector
          Vector3 proj_p = p - unit_n_vec[i] * p.dot(unit_n_vec[i]);

          //Calculate angles with a as the reference jet
          double theta_p = atan2(a_vec[i].cross(proj_p).dot(unit_n_vec[i]), a_vec[i].dot(proj_p));
          //double theta_p = a_vec[i].angle(proj_p);
          //std::cout << theta_p << "\t" << atan2(a_vec[i].cross(proj_p).dot(unit_n_vec[i]), a_vec[i].dot(proj_p)) << std::endl;

          //Find transverse momentum 
          double pT_to_plane = std::abs(p.dot(unit_n_vec[i]));

          //if its in the region, have sorter write it down
          if ( (0 < theta_p && theta_p < thetaRef_vec[i])){
            sorter[i] = pT_to_plane;
            rescaledAngle[i] = (theta_p/thetaRef_vec[i]);
          }
        }

        double minVal = 1000;
        int jetNum = -1;
        double bestRescaledAngle = -1000;
        for(int k=0; k<4; ++k){
          if(sorter[k]!=0 && minVal>sorter[k]){
            minVal = sorter[k];
            bestRescaledAngle = rescaledAngle[k];
            jetNum = k;
          }
        }
        if(jetNum ==-1) continue;

        particleSort[j]++;
        phi_rescaled->fill(bestRescaledAngle+jetNum);

        //Check if inside or outside region and fill corresponding histo
        //regionA
        if (jetNum==0){
          inside1->fill(bestRescaledAngle);
          insideTotal->fill(bestRescaledAngle);
          if(allJetsMatched){
            inside1_allMatched->fill(bestRescaledAngle);
            insideTotal_allMatched->fill(bestRescaledAngle);
          }
        }
        //RegionB
        if (jetNum==1){
          outside1->fill(bestRescaledAngle);
          outsideTotal->fill(bestRescaledAngle);
          if(allJetsMatched){
            outside1_allMatched->fill(bestRescaledAngle);
            outsideTotal_allMatched->fill(bestRescaledAngle);
          }
        }
        //RegionC
        if (jetNum==2){
          inside2->fill(bestRescaledAngle);
          insideTotal->fill(bestRescaledAngle);
          if(allJetsMatched){
            inside2_allMatched->fill(bestRescaledAngle);
            insideTotal_allMatched->fill(bestRescaledAngle);
          }
        }
        //RegionD
        if (jetNum==3){
          outside2->fill(bestRescaledAngle);
          outsideTotal->fill(bestRescaledAngle);
          if(allJetsMatched){
            outside2_allMatched->fill(bestRescaledAngle);
            outsideTotal_allMatched->fill(bestRescaledAngle);
          }
        }
      }

      for (int i = 0; i<particles.size(); i++){
        sorted->fill(particleSort[i]);
      }
      writeEvent(correctJetPair1, correctJetPair2, fv_Wp, fv_Wm, fv_q1, fv_q2, fv_q3, fv_q4, massWp, massWm, selectedJets, ordered_indices);

    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      normalize(_h_dijetMass1);
      normalize(_h_dijetMass2);
      normalize(_h_dijetMass1_allMatched);
      normalize(_h_dijetMass2_allMatched);
      normalize(_h_ycutAll);
      normalize(_h_selectedJets);
      normalize(_h_mult_before);
      normalize(_h_mult_after);
      normalize(isMatched);
      normalize(orderingChoice);
      normalize(orderingChoice_allMatched);
      scale(phi_rescaled, 1.0 / sumOfWeights());

      scale(inside1, 1.0 / sumOfWeights());
      scale(outside1, 1.0 / sumOfWeights());
      scale(inside1_allMatched, 1.0 / sumOfWeights());
      scale(outside1_allMatched, 1.0 / sumOfWeights());

      scale(inside2, 1.0 / sumOfWeights());
      scale(outside2, 1.0 / sumOfWeights());
      scale(inside2_allMatched, 1.0 / sumOfWeights());
      scale(outside2_allMatched, 1.0 / sumOfWeights());

      scale(insideTotal, 1.0 / sumOfWeights());
      scale(outsideTotal, 1.0 / sumOfWeights());
      scale(insideTotal_allMatched, 1.0 / sumOfWeights());
      scale(outsideTotal_allMatched, 1.0 / sumOfWeights());

      //Find actual ratio for figure 6
      divide(*insideTotal, *outsideTotal, region_ratio);
      divide(*inside1, *inside2, inside_ratio);
      divide(*outside1, *outside2, outside_ratio);
    
      divide(*insideTotal_allMatched, *outsideTotal_allMatched, region_ratio_allMatched);
      divide(*inside1_allMatched, *inside2_allMatched, inside_ratio_allMatched);
      divide(*outside1_allMatched, *outside2_allMatched, outside_ratio_allMatched);
    
    }

    bool m_debug = false;

    // If you want to match LEP, you should use this value for the max track eta
    //const double m_maxEtaTracks = 1.738; // Using cos(theta) < 0.94, which this should correspond to

    // For these studies, currently studying the same eta cutoff for all paricles and for tracks
    // This choice is somewhat arbitrary, and we could probably go higher in eta with future detectors
    const double m_maxEtaTracks = 2.1;
    const double m_maxEta = 2.1;

    const double m_jetRadius = 0.7;

    Histo1DPtr _h_dijetMass1;
    Histo1DPtr _h_dijetMass2;

    Histo1DPtr _h_dijetMass1_allMatched;
    Histo1DPtr _h_dijetMass2_allMatched;

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
    Histo1DPtr outsideTotal;
    //Figure 5, but requiring match between each jet and a truth quark
    Histo1DPtr insideTotal_allMatched;
    Histo1DPtr outsideTotal_allMatched;

    Scatter2DPtr region_ratio;
    Scatter2DPtr inside_ratio;
    Scatter2DPtr outside_ratio;

    Scatter2DPtr region_ratio_allMatched;
    Scatter2DPtr inside_ratio_allMatched;
    Scatter2DPtr outside_ratio_allMatched;

    //Diagnostic plots
    Histo1DPtr _h_mult_nocut;

    Histo1DPtr inside1;
    Histo1DPtr inside2;
    Histo1DPtr outside1;
    Histo1DPtr outside2;

    Histo1DPtr inside1_allMatched;
    Histo1DPtr inside2_allMatched;
    Histo1DPtr outside1_allMatched;
    Histo1DPtr outside2_allMatched;

    Histo1DPtr regionA;
    Histo1DPtr regionB;
    Histo1DPtr regionC;
    Histo1DPtr regionD;
    std::vector<Histo1DPtr> regions;

    Histo1DPtr quark1;
    Histo1DPtr quark2;
    Histo1DPtr quark3;
    Histo1DPtr quark4;
    std::vector<Particle> truth_quarks;

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
    std::vector<Histo1DPtr> quarkdR;

    Histo1DPtr regionQA;
    Histo1DPtr regionQB;
    Histo1DPtr regionQC;
    Histo1DPtr regionQD;
    std::vector<Histo1DPtr> regionsQ;

    Histo1DPtr orderingChoice;
    Histo1DPtr orderingChoice_allMatched;
    Histo1DPtr isMatched;

    Histo1DPtr sorted;
    std::vector<Histo1DPtr> histos;
    std::vector<Histo1DPtr> histosRatio;

    std::ofstream* out;
  };


  RIVET_DECLARE_PLUGIN(MC_COLORRECONNECTION);

}


