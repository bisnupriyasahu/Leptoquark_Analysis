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



float WScaleFactor=1.0; // computed on Apr 20th
//float TTScaleFactor=0.91;
//float WScaleFactor=1.22;
//float TTScaleFactor=0.856;



float TT_FulLep_BR= 0.1061;
float TT_SemiLep_BR= 0.4392;
float TT_Had_BR= 0.4544;


float XSection(std::string OutName) {
    
    
  ////////////////////////////////////////////////////////////////
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    
  ////////////////////////////////////////////////////////////////
    
  //    https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=398123591
    
    
    
    
    
    
    
  if (OutName.find("WJetsToLNu_Inc") != string::npos) return 50690;   // As we have large cut at Skim, this one is not needed
  else if (OutName.find("WJetsToLNu_HT-70to100") != string::npos) return 1372;
  else if (OutName.find("WJetsToLNu_HT-100to200") != string::npos) return 1343;
  else if (OutName.find("WJetsToLNu_HT-200to400") != string::npos) return 359.6;
  else if (OutName.find("WJetsToLNu_HT-400to600") != string::npos) return 48.85;
  else if (OutName.find("WJetsToLNu_HT-600to800") != string::npos) return 12.05;
  else if (OutName.find("WJetsToLNu_HT-800to1200") != string::npos) return 5.501;
  else if (OutName.find("WJetsToLNu_HT-1200to2500") != string::npos) return 1.329;
  else if (OutName.find("WJetsToLNu_HT-2500toInf") != string::npos) return 0.03216;
    
  else if (OutName.find("W1JetsToLNu") != string::npos) return 9644.5;
  else if (OutName.find("W2JetsToLNu") != string::npos) return 3144.5;
  else if (OutName.find("W3JetsToLNu") != string::npos) return 954.8;
  else if (OutName.find("W4JetsToLNu") != string::npos) return 485.6;
    
    
    
  //http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2016/204
  else if (OutName.find("WToMuNu_M-100_") != string::npos) return 163.15;
  else if (OutName.find("WToMuNu_M-200_") != string::npos) return 6.236;
  else if (OutName.find("WToMuNu_M-500_") != string::npos) return 0.2138;
  else if (OutName.find("WToMuNu_M-1000_") != string::npos) return 0.01281;
  else if (OutName.find("WToMuNu_M-2000_") != string::npos) return 5.56e-04;
  //    else if (OutName.find("WToMuNu_M-3000_") != string::npos) return 2.904e-05;
    
    
    
  else if (OutName.find("WToTauNu_M-100_") != string::npos) return 163.15;
  else if (OutName.find("WToTauNu_M-200_") != string::npos) return 6.236;
  else if (OutName.find("WToTauNu_M-500_") != string::npos) return 0.2138;
  else if (OutName.find("WToTauNu_M-1000_") != string::npos) return 0.01281;
  else if (OutName.find("WToTauNu_M-2000_") != string::npos) return 5.56e-04;
    
    
    
  //    else if (OutName.find("DYJetsToLL_Inc") != string::npos) return 4895 * 1.012; // As we have large cut at Skim, this one is not needed
  else if (OutName.find("DYJetsToLL_M-50_HT-70to100") != string::npos) return 175.3; // from EXO-16-049
  else if (OutName.find("DYJetsToLL_M-50_HT-100to200") != string::npos) return 148;
  else if (OutName.find("DYJetsToLL_M-50_HT-200to400") != string::npos) return 40.94;
  else if (OutName.find("DYJetsToLL_M-50_HT-400to600") != string::npos) return 5.497;
  else if (OutName.find("DYJetsToLL_M-50_HT-600to800") != string::npos) return 1.354;
  else if (OutName.find("DYJetsToLL_M-50_HT-800to1200") != string::npos) return 0.625;
  else if (OutName.find("DYJetsToLL_M-50_HT-1200to2500") != string::npos) return 0.151;
  else if (OutName.find("DYJetsToLL_M-50_HT-2500toInf") != string::npos) return 0.003647;
    
    
  //Di-boson   Pythia is not useful in this analysis
  //    else if (OutName.find("WW_pythia") != string::npos) return 115.0;
  //    else if (OutName.find("WZ_pythia") != string::npos) return 47.13;
  //    else if (OutName.find("ZZ_pythia") != string::npos) return 16.523;
    
    

    
    
  else if (OutName.find("WWTo2L2Nu_powheg") != string::npos) return  12.178 ;
  else if (OutName.find("WWTo4Q_powheg") != string::npos) return  51.723 ;
  else if (OutName.find("WWTo1LNuQQ_powheg") != string::npos) return  49.997 ;
  //    else if (OutName.find("WWTo1L1Nu2Q_amcatnloFXFX_madspin") != string::npos) return 49.997  ;
    
                           
  else if (OutName.find("WZTo1L3Nu_amcatnloFXFX") != string::npos) return 3.033e+00;  //Not available  Just added after CWR
  else if (OutName.find("WZTo2L2Q_amcNLO") != string::npos) return  5.595 ;
  //    else if (OutName.find("WZTo2Q2Nu_amcatnloFXFX") != string::npos) return 10.000; // NOTE THIS IS JUST MY ESTIMATION small effect Not available
  else if (OutName.find("WZTo3LNu_amcNLO") != string::npos) return  4.42965 ;
  //       else if (OutName.find("WZToLNu2Q_powheg") != string::npos) return  10.71 ;
  else if (OutName.find("WZTo1L1Nu2Q_amcNLO") != string::npos) return  10.71 ;
  // The missing one is WZto4Q

    


    

    
    
    
  else if (OutName.find("ZZTo2L2Nu_powheg") != string::npos) return  0.564 ;
  //    else if (OutName.find("ZZTo2L2Q_powheg") != string::npos) return  3.22 ;
  else if (OutName.find("ZZTo2L2Q_amcNLO") != string::npos) return  3.22 ;
  //    else if (OutName.find("ZZTo2Q2Nu_powheg") != string::npos) return  4.04 ;
  else if (OutName.find("ZZTo4L_powheg") != string::npos) return  1.212 ;
  //    else if (OutName.find("ZZTo4Q_amcatnloFXFX") != string::npos) return 6.000;   // NOTE THIS IS JUST MY ESTIMATION
  //    else if (OutName.find("ZZTo2L2Q_13TeV_amcatnloFXFX_madspin") != string::npos) return  3.22 ;
  // The missing one is ZZto4Nu
    
    
  //    else if (OutName.find("WZTo2L2Q_13TeV_amcatnloFXFX_madspin") != string::npos) return  5.595 ;
  //    else if (OutName.find("WZTo3LNu_amcNLO") != string::npos) return   4.42965;
    
    
    
    
  //Need to check these
    
    
  //    else if (OutName.find("") != string::npos) return   ;
    
    

    
  //SingleTop
  else if (OutName.find("ST_t-channel_antitop") != string::npos) return 80.95;
  else if (OutName.find("ST_t-channel_top") != string::npos) return 136.02;
  else if (OutName.find("ST_tW_antitop_5f") != string::npos) return 35.6;
  else if (OutName.find("ST_tW_top_5f") != string::npos) return 35.6;
  else if (OutName.find("ST_s_channel") != string::npos) return 3.36 ;
    
    
  else if (OutName.find("TTJets_DiLeptonic") != string::npos) return (831.76*TT_FulLep_BR);
  else if (OutName.find("TTJets_Hadronic") != string::npos) return (831.76*TT_Had_BR);
  else if (OutName.find("TTJets_SemiLeptonic") != string::npos) return (831.76*TT_SemiLep_BR);
    
  //    else if (OutName.find("TTJets") != string::npos) return (831.76);
    
    
    
    
    
    
    
  //    //    https://twiki.cern.ch/twiki/bin/view/CMS/Exo2015LQ1AndLQ2Analyses
  //    else if (OutName.find("skimed_lq200") != string::npos) return 60.6;
  //    else if (OutName.find("skimed_lq250") != string::npos) return 20.3;
  //    else if (OutName.find("skimed_lq300") != string::npos) return     8.05E+00;
  //    else if (OutName.find("skimed_lq350") != string::npos) return     3.58E+00;
  //    else if (OutName.find("skimed_lq400") != string::npos) return     1.74E+00;
  //    else if (OutName.find("skimed_lq450") != string::npos) return     9.05E-01;
  //    else if (OutName.find("skimed_lq500") != string::npos) return     4.96E-01;
  //    else if (OutName.find("skimed_lq550") != string::npos) return     2.84E-01;
  //    else if (OutName.find("skimed_lq600") != string::npos) return     1.69E-01;
  //    else if (OutName.find("skimed_lq650") != string::npos) return     1.03E-01;
  //    else if (OutName.find("skimed_lq700") != string::npos) return     6.48E-02;
  //    else if (OutName.find("skimed_lq750") != string::npos) return     4.16E-02;
  //    else if (OutName.find("skimed_lq800") != string::npos) return     2.73E-02;
  //    else if (OutName.find("skimed_lq850") != string::npos) return     1.82E-02;
  //    else if (OutName.find("skimed_lq900") != string::npos) return     1.23E-02;
  //    else if (OutName.find("skimed_lq950") != string::npos) return     8.45E-03;
  //    else if (OutName.find("skimed_lq1000") != string::npos) return     5.86E-03;
  //    else if (OutName.find("skimed_lq1050") != string::npos) return     4.11E-03;
  //    else if (OutName.find("skimed_lq1100") != string::npos) return     2.91E-03;
  //    else if (OutName.find("skimed_lq1150") != string::npos) return     2.08E-03;
  //    else if (OutName.find("skimed_lq1200") != string::npos) return     1.50E-03;
  //    else if (OutName.find("skimed_lq1250") != string::npos) return     1.09E-03;
  //    else if (OutName.find("skimed_lq1300") != string::npos) return     7.95E-04;
  //    else if (OutName.find("skimed_lq1350") != string::npos) return     5.85E-04;
  //    else if (OutName.find("skimed_lq1400") != string::npos) return     4.33E-04;
  //    else if (OutName.find("skimed_lq1450") != string::npos) return     3.21E-04;
  //    else if (OutName.find("skimed_lq1500") != string::npos) return     2.40E-04;
  //
  //
  //
    
  else if (OutName.find("Codex") != string::npos ) return      1.0;
    
  else if (OutName.find("EWK_DYToLL") != string::npos ) return      3.987;
    
    
  else if (OutName.find("QCD_Pt-20toInf_MuEnrichedPt15") != string::npos) return     720648000  * 0.00042 ;
  else if (OutName.find("QCD") != string::npos) return     720648000  * 0.00042 ;
    
  //    https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#W_jets
    
  else if (OutName.find("WJetsToLNu_FXFX") != string::npos) return  61526.7   ;
  else if (OutName.find("WJetsToLNu_Pt-50to100_FXFX") != string::npos) return  8053   ;
  else if (OutName.find("WJetsToLNu_Pt-100to250_FXFX") != string::npos) return  676.3   ;
  else if (OutName.find("WJetsToLNu_Pt-250to400_FXFX") != string::npos) return  23.94   ;
  else if (OutName.find("WJetsToLNu_Pt-400to600_FXFX") != string::npos) return  3.031   ;
  else if (OutName.find("WJetsToLNu_Pt-600toInf_FXFX") != string::npos) return  0.4524   ;
    
    
    
  else if (OutName.find("DYJetsToLL_M-50_FXFX") != string::npos) return          5765.4 ;
  else if (OutName.find("DYJetsToLL_Pt-100to250_FXFX") != string::npos) return   83.12 ;
  else if (OutName.find("DYJetsToLL_Pt-250to400_FXFX") != string::npos) return   3.047 ;
  else if (OutName.find("DYJetsToLL_Pt-400to650_FXFX") != string::npos) return   0.3921 ;
  else if (OutName.find("DYJetsToLL_Pt-650toInf_FXFX") != string::npos) return   0.03636 ;
    
  else if (OutName.find("LQ") != string::npos ) return      1.0;
    
    
    
  else if (OutName.find("ZJetsToNuNu_HT-100to200") != string::npos) return    280.35  * 1.23       ;
  else if (OutName.find("ZJetsToNuNu_HT-200to400") != string::npos) return     77.67  * 1.23     ;
  else if (OutName.find("ZJetsToNuNu_HT-400to600") != string::npos) return      10.73 * 1.23     ;
  else if (OutName.find("ZJetsToNuNu_HT-600to800") != string::npos) return      2.559  * 1.23     ;
  else if (OutName.find("ZJetsToNuNu_HT-800to1200") != string::npos) return     1.1796   * 1.23    ;
  else if (OutName.find("ZJetsToNuNu_HT-1200to2500") != string::npos) return     0.28833   * 1.23    ;
  else if (OutName.find("ZJetsToNuNu_HT-2500toInf") != string::npos) return      0.006945    * 1.23     ;
    
    
    
    
    
  else {
    cout<<"\n\n*********\nNot Listed in XSection menu !!!! Watch cout    "<<OutName<< "\n\n*********\n";
    return 1;
  }
}

vector <float> W_HTBin(std::string FileLoc){

  const int WSize=9;
  std::string W_ROOTFiles[WSize]={"WJetsToLNu_Inc.root","WJetsToLNu_Inc.root", "WJetsToLNu_HT-100to200.root","WJetsToLNu_HT-200to400.root","WJetsToLNu_HT-400to600.root", "WJetsToLNu_HT-600to800.root","WJetsToLNu_HT-800to1200.root","WJetsToLNu_HT-1200to2500.root","WJetsToLNu_HT-2500toInf.root"};
  std::cout<<"cominginsie whtbin1"<<std::endl;   
    
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
  //  int W_events = 9901992;  
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
    
  //  vector<float> W_events;
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
//float weightCalc(TH1F *Histo,std::string outputName, float genHT,int W_HTbin, float WMass, int WMuNu_MassBin, int WTauNu_Massbin) {
    
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




  //##################################################################
  //   Stitching for  W sample
  //##################################################################
    
  if (isWjet != string::npos && WMass <= 100){  //FIXME
    if (genHT <= 70)  return  luminosity * XSection("WJetsToLNu_Inc") / W_HTbin[0];
    else if (genHT > 70 && genHT <= 100)  return   luminosity * XSection("WJetsToLNu_Inc") / W_HTbin[1];
    else if (genHT > 100 && genHT <= 200)  return   luminosity * XSection("WJetsToLNu_HT-100to200") / W_HTbin[2];
    else if (genHT > 200 && genHT <= 400)  return   luminosity * XSection("WJetsToLNu_HT-200to400") / W_HTbin[3];
    else if (genHT > 400 && genHT <= 600)  return   luminosity * XSection("WJetsToLNu_HT-400to600") / W_HTbin[4];
    else if (genHT > 600 && genHT <= 800)  return   luminosity * XSection("WJetsToLNu_HT-600to800") / W_HTbin[5];
    else if (genHT > 800 && genHT <= 1200)  return   luminosity * XSection("WJetsToLNu_HT-800to1200") / W_HTbin[6];
    else if (genHT > 1200 && genHT <= 2500)  return   luminosity * XSection("WJetsToLNu_HT-1200to2500") / W_HTbin[7];
    else if (genHT > 2500 )  return   luminosity * XSection("WJetsToLNu_HT-2500toInf") / W_HTbin[8];
    else   {cout<<"**********   wooow  ********* There is a problem here\n";return 0;}
  }    
  //    if (isWjet != string::npos && WMass <= 100){  FIXME
  //        if (genHT <= 70)  return  luminosity * XSection("WJetsToLNu_Inc") / W_HTbin[0];
  //        else if (genHT > 70 && genHT <= 100)  return   luminosity * XSection("WJetsToLNu_HT-70to100") / W_HTbin[1];
  //        else if (genHT > 100 && genHT <= 200)  return   luminosity * XSection("WJetsToLNu_HT-100to200") / W_HTbin[2];
  //        else if (genHT > 200 && genHT <= 400)  return   luminosity * XSection("WJetsToLNu_HT-200to400") / W_HTbin[3];
  //        else if (genHT > 400 && genHT <= 600)  return   luminosity * XSection("WJetsToLNu_HT-400to600") / W_HTbin[4];
  //        else if (genHT > 600 && genHT <= 800)  return   luminosity * XSection("WJetsToLNu_HT-600to800") / W_HTbin[5];
  //        else if (genHT > 800 && genHT <= 1200)  return   luminosity * XSection("WJetsToLNu_HT-800to1200") / W_HTbin[6];
  //        else if (genHT > 1200 && genHT <= 2500)  return   luminosity * XSection("WJetsToLNu_HT-1200to2500") / W_HTbin[7];
  //        else if (genHT > 2500 )  return   luminosity * XSection("WJetsToLNu_HT-2500toInf") / W_HTbin[8];
  //        else   {cout<<"**********   wooow  ********* There is a problem here\n";return 0;}
  //    }
  //    
    
  //    if (isWjet != string::npos && WMass <= 100){
  //        
  //        if (genNumJet_ == 0) return  luminosity * XSection("WJetsToLNu_Inc") / W_HTbin[0];
  //        else if (genNumJet_ == 1) return  luminosity * XSection("W1JetsToLNu") / W_HTbin[1];
  //        else if (genNumJet_ == 2) return  luminosity * XSection("W2JetsToLNu") / W_HTbin[2];
  //        else if (genNumJet_ == 3) return  luminosity * XSection("W3JetsToLNu") / W_HTbin[3];
  //        else if (genNumJet_ == 4) return  luminosity * XSection("W4JetsToLNu") / W_HTbin[4];
  //                else   {cout<<"**********   wooow  ********* There is a problem here\n";return 0;}
  //}

    
    
  else if (isWToMuNu != string::npos){
    if  (WMass > 100 && WMass <= 200)  return   luminosity /  (WMuNu_MassBin[0] /XSection("WToMuNu_M-100_"));
    else if  (WMass > 200 && WMass <= 500)  return   luminosity /  (WMuNu_MassBin[0] /XSection("WToMuNu_M-100_") +  WMuNu_MassBin[1]/ XSection("WToMuNu_M-200_")) ;
    else if  (WMass > 500 && WMass <= 1000)  return   luminosity /  (WMuNu_MassBin[0] /XSection("WToMuNu_M-100_") +  WMuNu_MassBin[1]/ XSection("WToMuNu_M-200_") +  WMuNu_MassBin[2]/ XSection("WToMuNu_M-500_")) ;
    else if  (WMass > 1000 && WMass <= 2000)  return   luminosity /  (WMuNu_MassBin[0] /XSection("WToMuNu_M-100_") +  WMuNu_MassBin[1]/ XSection("WToMuNu_M-200_") +  WMuNu_MassBin[2]/ XSection("WToMuNu_M-500_") + WMuNu_MassBin[3]/ XSection("WToMuNu_M-1000_")) ;
    else if  (WMass > 2000 )  return   luminosity / (WMuNu_MassBin[0] /XSection("WToMuNu_M-100_") +  WMuNu_MassBin[1]/ XSection("WToMuNu_M-200_") +  WMuNu_MassBin[2]/ XSection("WToMuNu_M-500_") + WMuNu_MassBin[3]/ XSection("WToMuNu_M-1000_") + WMuNu_MassBin[4]/ XSection("WToMuNu_M-2000_")) ;
    else   {cout<<"**********   wooow1  ********* There is a problem here\n";return 0;}
  }
    
    
    
    
    
    
  else if (isWToTauNu != string::npos){
    if  (WMass > 100 && WMass <= 200)  return   luminosity /  (WTauNu_Massbin[0] /XSection("WToMuNu_M-100_"));
    else if  (WMass > 200 && WMass <= 1000)  return   luminosity /  (WTauNu_Massbin[0] /XSection("WToMuNu_M-100_") +  WTauNu_Massbin[1]/ XSection("WToMuNu_M-200_")) ;
    else if  (WMass > 1000 && WMass <= 2000)  return   luminosity /  (WTauNu_Massbin[0] /XSection("WToMuNu_M-100_") +  WTauNu_Massbin[1]/ XSection("WToMuNu_M-200_") +  WTauNu_Massbin[2]/ XSection("WToMuNu_M-1000_")) ;
    else if  (WMass > 2000 )  return   luminosity / (WTauNu_Massbin[0] /XSection("WToMuNu_M-100_") +  WTauNu_Massbin[1]/ XSection("WToMuNu_M-200_") +   WTauNu_Massbin[2]/ XSection("WToMuNu_M-1000_") + WTauNu_Massbin[3]/ XSection("WToMuNu_M-2000_")) ;
    else   {cout<<"**********   wooow  ********* There is a problem here\n";return 0;}
  }
    
    
    
  else
    return luminosity * XSection(newOut)*1.0 / Histo->GetBinContent(2);
}

