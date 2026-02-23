#include "Riostream.h"
#include <string.h>
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;


// Reads input from the .dat files into root files
// The file list is taken from inputFiles.txt
// and the file format is inFileName outFileName

int main(int argc,char *argv[]){
  // Need to create a directory for the output files
  std::string outDir = "rootFiles";
  if (mkdir(outDir.c_str(), 0777) == -1){
    std::cout << "Did not make directory " << outDir << std::endl;
    cerr << "Directory creation failed with error :  " << strerror(errno) << endl;
  }
  else{
    cout << "Directory created";
  }
  

  std::string filename = "inputFiles.txt";
  std::fstream infile(filename.c_str());
  std::cout << "Reading input and output file names from " << filename << std::endl;
  std::cout << "Writing output files to " << outDir << " directory" << std::endl;

  std::string treename = "tree";
  std::string mcFile, outName;
  bool isFirst = true;
  while (infile >> mcFile >> outName){
   if( mcFile.at(0) == '#' ){
     continue;
   }
   std::cout << "Reading " << mcFile << " and printing to rootFiles/" << outName << std::endl;

   ifstream in;
   in.open(mcFile);
   
   auto f = TFile::Open(Form("%s/%s", outDir.c_str(), outName.c_str()),"RECREATE");
   TTree *newtree= new  TTree ("tree","");

   bool isJet1 = false;
   bool isJet2 = false;
   bool isJet3 = false;
   bool isJet4 = false;
   bool isConstit1 = false;
   bool isConstit2 = false;
   bool isConstit3 = false;
   bool isConstit4 = false;

   // Open the ROOT file.
   std::vector<double> jetconstituents_4_theta_1, jetconstituents_4_phi_1, jetconstituents_4_energy_1;
   std::vector<double> jetconstituents_4_theta_2, jetconstituents_4_phi_2, jetconstituents_4_energy_2;
   std::vector<double> jetconstituents_4_theta_3, jetconstituents_4_phi_3, jetconstituents_4_energy_3;
   std::vector<double> jetconstituents_4_theta_4, jetconstituents_4_phi_4, jetconstituents_4_energy_4;

   std::vector<double> Wp_theta, Wp_phi, Wp_energy;
   std::vector<double> Wm_theta, Wm_phi, Wm_energy;
   std::vector<double> truth_Wp_theta, truth_Wp_phi, truth_Wp_energy;
   std::vector<double> truth_Wm_theta, truth_Wm_phi, truth_Wm_energy;
   std::vector<double> truth_Wpq1_theta, truth_Wpq1_phi, truth_Wpq1_energy;
   std::vector<double> truth_Wpq2_theta, truth_Wpq2_phi, truth_Wpq2_energy;
   std::vector<double> truth_Wmq1_theta, truth_Wmq1_phi, truth_Wmq1_energy;
   std::vector<double> truth_Wmq2_theta, truth_Wmq2_phi, truth_Wmq2_energy;

   newtree->Branch("Wppart_theta",&Wp_theta);
   newtree->Branch("Wppart_phi",&Wp_phi);
   newtree->Branch("Wppart_energy",&Wp_energy);

   newtree->Branch("truth_Wp_theta",&truth_Wp_theta);
   newtree->Branch("truth_Wp_phi",&truth_Wp_phi);
   newtree->Branch("truth_Wp_energy",&truth_Wp_energy);

   newtree->Branch("truth_Wpq1_theta",&truth_Wpq1_theta);
   newtree->Branch("truth_Wpq1_phi",&truth_Wpq1_phi);
   newtree->Branch("truth_Wpq1_energy",&truth_Wpq1_energy);

   newtree->Branch("truth_Wpq2_theta",&truth_Wpq2_theta);
   newtree->Branch("truth_Wpq2_phi",&truth_Wpq2_phi);
   newtree->Branch("truth_Wpq2_energy",&truth_Wpq2_energy);


   newtree->Branch("Wmpart_theta",&Wm_theta);
   newtree->Branch("Wmpart_phi",&Wm_phi);
   newtree->Branch("Wmpart_energy",&Wm_energy);

   newtree->Branch("truth_Wm_theta",&truth_Wm_theta);
   newtree->Branch("truth_Wm_phi",&truth_Wm_phi);
   newtree->Branch("truth_Wm_energy",&truth_Wm_energy);

   newtree->Branch("truth_Wmq1_theta",&truth_Wmq1_theta);
   newtree->Branch("truth_Wmq1_phi",&truth_Wmq1_phi);
   newtree->Branch("truth_Wmq1_energy",&truth_Wmq1_energy);

   newtree->Branch("truth_Wmq2_theta",&truth_Wmq2_theta);
   newtree->Branch("truth_Wmq2_phi",&truth_Wmq2_phi);
   newtree->Branch("truth_Wmq2_energy",&truth_Wmq2_energy);


   newtree->Branch("jetconstituents_4_theta_1",&jetconstituents_4_theta_1);
   newtree->Branch("jetconstituents_4_phi_1",&jetconstituents_4_phi_1);
   newtree->Branch("jetconstituents_4_energy_1",&jetconstituents_4_energy_1);

   newtree->Branch("jetconstituents_4_theta_2",&jetconstituents_4_theta_2);
   newtree->Branch("jetconstituents_4_phi_2",&jetconstituents_4_phi_2);
   newtree->Branch("jetconstituents_4_energy_2",&jetconstituents_4_energy_2);

   newtree->Branch("jetconstituents_4_theta_3",&jetconstituents_4_theta_3);
   newtree->Branch("jetconstituents_4_phi_3",&jetconstituents_4_phi_3);
   newtree->Branch("jetconstituents_4_energy_3",&jetconstituents_4_energy_3);

   newtree->Branch("jetconstituents_4_theta_4",&jetconstituents_4_theta_4);
   newtree->Branch("jetconstituents_4_phi_4",&jetconstituents_4_phi_4);
   newtree->Branch("jetconstituents_4_energy_4",&jetconstituents_4_energy_4);


   std::string name; 
   int cjet = 0;
   while (std::getline(in, name)) {
      // Most of these values are not stored, but create variable for them in case this changes
      double e, phi, theta;

      // Each event is written to the tree seperately, so the tree should be filled.
      // After this, reset all the information for the new event.
      if(name.find("Event") != std::string::npos){ 
        if(!isFirst){
        newtree->Fill (); // fill the new tree with the previous event's information
        }
        isFirst = false;

        isJet1 = false;
        isJet2 = false;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = false;

        Wp_theta.clear(); Wp_phi.clear(); Wp_energy.clear();;
        Wm_theta.clear(); Wm_phi.clear(); Wm_energy.clear();;
        truth_Wp_theta.clear(); truth_Wp_phi.clear(); truth_Wp_energy.clear();;
        truth_Wm_theta.clear(); truth_Wm_phi.clear(); truth_Wm_energy.clear();;
        truth_Wpq1_theta.clear(); truth_Wpq1_phi.clear(); truth_Wpq1_energy.clear();;
        truth_Wpq2_theta.clear(); truth_Wpq2_phi.clear(); truth_Wpq2_energy.clear();;
        truth_Wmq1_theta.clear(); truth_Wmq1_phi.clear(); truth_Wmq1_energy.clear();;
        truth_Wmq2_theta.clear(); truth_Wmq2_phi.clear(); truth_Wmq2_energy.clear();;


        jetconstituents_4_theta_1.clear();
        jetconstituents_4_phi_1.clear();
        jetconstituents_4_energy_1.clear();

        jetconstituents_4_theta_2.clear();
        jetconstituents_4_phi_2.clear();
        jetconstituents_4_energy_2.clear();

        jetconstituents_4_theta_3.clear();
        jetconstituents_4_phi_3.clear();
        jetconstituents_4_energy_3.clear();

        jetconstituents_4_theta_4.clear();
        jetconstituents_4_phi_4.clear();
        jetconstituents_4_energy_4.clear();

        std::getline(in, name);
        std::istringstream line1(name);
        line1 >> theta >> phi >> e;
        Wp_theta.push_back(theta); Wp_phi.push_back(phi); Wp_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line2(name);
        line2 >> theta >> phi >> e;
        truth_Wp_theta.push_back(theta); truth_Wp_phi.push_back(phi); truth_Wp_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line3(name);
        line3 >> theta >> phi >> e;
        truth_Wpq1_theta.push_back(theta); truth_Wpq1_phi.push_back(phi); truth_Wpq1_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line4(name);
        line4 >> theta >> phi >> e;
        truth_Wpq2_theta.push_back(theta); truth_Wpq2_phi.push_back(phi); truth_Wpq2_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line5(name);
        line5 >> theta >> phi >> e;
        Wm_theta.push_back(theta); Wm_phi.push_back(phi); Wm_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line6(name);
        line6 >> theta >> phi >> e;
        truth_Wm_theta.push_back(theta); truth_Wm_phi.push_back(phi); truth_Wm_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line7(name);
        line7 >> theta >> phi >> e;
        truth_Wmq1_theta.push_back(theta); truth_Wmq1_phi.push_back(phi); truth_Wmq1_energy.push_back(e);

        std::getline(in, name);
        std::istringstream line8(name);
        line8 >> theta >> phi >> e;
        truth_Wmq2_theta.push_back(theta); truth_Wmq2_phi.push_back(phi); truth_Wmq2_energy.push_back(e);

        continue;
      }

      // Switch to reading in jet information 
      if(name.find("Jet1") != std::string::npos){ 
        isJet1 = true;
        isJet2 = false;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = false;
        continue;
      }
      if(name.find("Jet2") != std::string::npos){
        isJet1 = false;
        isJet2 = true;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = false;
        continue;
      }
      if(name.find("Jet3") != std::string::npos){
        isJet1 = false;
        isJet2 = false;
        isJet3 = true;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = false;
        continue;
      }
      if(name.find("Jet4") != std::string::npos){
        isJet1 = false;
        isJet2 = false;
        isJet3 = false;
        isJet4 = true;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = false;
        continue;
      }
      // Switch to reading in constituent information
      if(name.find("Constit1") != std::string::npos){ 
        isJet1 = false;
        isJet2 = false;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = true;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = false;
        continue;
      }
      if(name.find("Constit2") != std::string::npos){
        isJet1 = false;
        isJet2 = false;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = true;
        isConstit3 = false;
        isConstit4 = false;
        continue;
      }
      if(name.find("Constit3") != std::string::npos){
        isJet1 = false;
        isJet2 = false;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = true;
        isConstit4 = false;
        continue;
      }
      if(name.find("Constit4") != std::string::npos){
        isJet1 = false;
        isJet2 = false;
        isJet3 = false;
        isJet4 = false;
        isConstit1 = false;
        isConstit2 = false;
        isConstit3 = false;
        isConstit4 = true;
        continue;
      }

      

      std::istringstream line(name);

      // Read in the constituent information
      if(isConstit1){
        line >> theta >> phi >> e;
        jetconstituents_4_phi_1.push_back(phi);
        jetconstituents_4_theta_1.push_back(theta);
        jetconstituents_4_energy_1.push_back(e);
        continue;
      }
      if(isConstit2){
        line >> theta >> phi >> e;
        jetconstituents_4_phi_2.push_back(phi);
        jetconstituents_4_theta_2.push_back(theta);
        jetconstituents_4_energy_2.push_back(e);
        continue;
      }
      if(isConstit3){
        line >> theta >> phi >> e;
        jetconstituents_4_phi_3.push_back(phi);
        jetconstituents_4_theta_3.push_back(theta);
        jetconstituents_4_energy_3.push_back(e);
        continue;
      }
      if(isConstit4){
        line >> theta >> phi >> e;
        jetconstituents_4_phi_4.push_back(phi);
        jetconstituents_4_theta_4.push_back(theta);
        jetconstituents_4_energy_4.push_back(e);
        continue;
     }
   }
 
   in.close();

   f->cd();
   newtree->Write();   //save new tree;
   f->Write();

  }
}

