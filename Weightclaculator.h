#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>

using namespace std;
vector <float> W_HTBin(std::string FileLoc){
  const int WSize=9;
  std::string W_ROOTFiles[WSize]={"WJetsToLNu_Inc.root","WJetsToLNu_Inc.root", "WJetsToLNu_HT-100to200.root","WJetsToLNu_HT-200to400.root","WJetsToLNu_HT-400to600.root", "WJetsToLNu_HT-600to800.root","WJetsToLNu_HT-800to1200.root","WJetsToLNu_HT-1200to2500.root","WJetsToLNu_HT-2500toInf.root"};
    
  //    std::string W_ROOTFiles[WSize]={"WJetsToLNu_Inc.root","WJetsToLNu_HT-70to100.root", "WJetsToLNu_HT-100to200.root","WJetsToLNu_HT-200to400.root","WJetsToLNu_HT-400to600.root", "WJetsToLNu_HT-600to800.root","WJetsToLNu_HT-800to1200.root","WJetsToLNu_HT-1200to2500.root","WJetsToLNu_HT-2500toInf.root"};

    
  //    const int WSize=1;
  //    std::string W_ROOTFiles[WSize]={"WJetsToLNu_Inc.root"};
    
  vector<float> W_events;
  W_events.clear();
    
  for (int i=0; i <WSize;i++){
        
    TFile * File_W = new TFile((FileLoc+W_ROOTFiles[i]).c_str());
    TH1F * Histo_W = (TH1F*) File_W->Get("hcount");
    W_events.push_back(Histo_W->GetBinContent(2));
    cout<<"Number of proccessed evenets for "<<W_ROOTFiles[i]<<" = "<<Histo_W->GetBinContent(2)<<"\n";
  }
    
  return W_events ;
    
}


vector <float> WMuNu_MassBin(std::string FileLoc){
    
  const int WSize=5;
  std::string W_ROOTFiles[WSize]={"WToMuNu_M-100_.root","WToMuNu_M-200_.root","WToMuNu_M-500_.root","WToMuNu_M-1000_.root","WToMuNu_M-2000_.root"};
    
  vector<float> W_events;
  W_events.clear();
    
  for (int i=0; i <WSize;i++){
        
    TFile * File_W = new TFile((FileLoc+W_ROOTFiles[i]).c_str());
    TH1F * Histo_W = (TH1F*) File_W->Get("hcount");
    W_events.push_back(Histo_W->GetBinContent(2));
    cout<<"Number of proccessed evenets for "<<W_ROOTFiles[i]<<" = "<<Histo_W->GetBinContent(2)<<"\n";
  }
    
  return W_events ;
    
}




vector <float> WTauNu_MassBin(std::string FileLoc){
    
  const int WSize=4;
  //    std::string W_ROOTFiles[WSize]={"WToTauNu_M-100_.root","WToTauNu_M-200_.root", "WToTauNu_M-500_.root","WToTauNu_M-1000_.root","WToTauNu_M-2000_.root"};
  std::string W_ROOTFiles[WSize]={"WToTauNu_M-100_.root","WToTauNu_M-200_.root","WToTauNu_M-1000_.root","WToTauNu_M-2000_.root"};
    
  vector<float> WTauNu_events;
  WTauNu_events.clear();
    
  for (int i=0; i <WSize;i++){
        
    TFile * File_W = new TFile((FileLoc+W_ROOTFiles[i]).c_str());
    TH1F * Histo_W = (TH1F*) File_W->Get("hcount");
    WTauNu_events.push_back(Histo_W->GetBinContent(2));
    cout<<"Number of proccessed evenets for "<<W_ROOTFiles[i]<<" = "<<Histo_W->GetBinContent(2)<<"\n";
  }
    
  return WTauNu_events ;
    
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////                    //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

float weightCalc(TH1F *Histo,std::string outputName, float genHT,vector<float> W_HTbin, float WMass, vector<float> WMuNu_MassBin, vector<float> WTauNu_Massbin) {
    
  stringstream ss(outputName);
    
  string token;
  string M;
  while (getline(ss,token, '/'))  M=token;
    
  std::string FirstPart = "";
  std::string LastPart = ".root";
  std::string newOut = M.substr(FirstPart.size());
  newOut = newOut.substr(0, newOut.size() - LastPart.size());
    
    
  //    float LOtoNLO_DY = 1.230888662;
  float LOtoNLO_DY = 1; // Now we boson have pt dependent SF
  //    float LOtoNLO_W = 1.213783784;
  float LOtoNLO_W = 1;  // Now we boson have pt dependent SF
  //    float luminosity=    12900;
  //    float luminosity=    35867;
  float luminosity=    41530; //May 15th
    
    
  size_t isSingleMu = outputName.find("SingleMu");
  size_t isSingleEle = outputName.find("SingleEle");
  size_t isData = outputName.find("Data");
    
    
  if (isSingleMu != string::npos || isSingleEle!= string::npos || isData !=string::npos)   return 1;
    
    
  //    size_t isWjet = outputName.find("WJets");
  size_t isWjet = outputName.find("JetsToLNu");
  size_t isWToMuNu = outputName.find("WToMuNu");
  size_t isWToTauNu = outputName.find("WToTauNu");
}

