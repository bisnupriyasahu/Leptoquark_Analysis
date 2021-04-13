#define MuonClass_cxx
#include "MuonClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <vector>
#include "TStopwatch.h"
#include <cstring>
#include <list>
#include <TStyle.h>

#include <string>
#include <ostream>
using namespace std;

int main(int argc, const char* argv[])
{

  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
      return 1;
    }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
      return 1;
    }
  MuonClass t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);

  return 0;


}

void MuonClass::Loop(Long64_t maxEvents, int reportEvery)
{
  //  const char* file2;
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;

  auto numOf_c_quark=0;
  auto numOf_s_quark=0;

  myMap1 = new std::map<std::string, TH1F*>();
  myMap2 = new map<string, TH2F*>();

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "Total entries: " << nentries << std::endl;
  Long64_t nentriesToCheck = nentries;

  
  //########################################                                                                                                                                                          
  // Muon Id, Iso, Trigger and Tracker Eff files                                                                                                                                                      
  //########################################                                                                                                                                                          
  TH2F** HistoMuId=FuncHistMuId();
  TH2F** HistoMuIso=FuncHistMuIso();
  TH2F** HistoMuTrg=FuncHistMuTrigger();
  TGraphAsymmErrors * HistoMuTrack=FuncHistMuTrack();
  //########################################                                                                                                                                                          
  // Electron MVA IdIso files                                                                                                                                                                         
  //########################################                                                                                                                                                          
  /* TH2F * HistoEleMVAIdIso90= FuncHistEleMVAId("Tot");
  TH2F * HistoEleMVAIdIso90_EffMC= FuncHistEleMVAId("MC");
  TH2F * HistoEleMVAIdIso90_EffData= FuncHistEleMVAId("Data");
  */
  //########################################                                                                                                                                                          
  // W and DY K-factor files  (Bin-based K-factor)                                                                                                                                                    
  //########################################                                                                                                                                                          
  std::string ROOTLocHT= "/nfs_scratch/bsahu/LQ_Analysis/CMSSW_10_2_18/src/MuonAnalyzer/Muon_Analyzer/test/";
  //   vector<float> W_HTBinROOTFiles = W_HTBin(ROOTLocHT);                                                                                                                                           
  //  std::cout<<"COMING AFTER WJETS"<<std::endl;                                                                                                                                                     
  // vector<float> WMuNu_MassBinROOTFiles = WMuNu_MassBin(ROOTLocHT);                                                                                                                                 
  // vector<float> WTauNu_MassBinROOTFiles = WTauNu_MassBin(ROOTLocHT);                                                                                                                               
  TFile * MassDepKFactor=TFile::Open("k_fakNNLO_use.root");
  TH1F* HistMassDepKFactor= (TH1F*) MassDepKFactor->Get("k_fak_mean");
  
  
  //########################################                                                                                                                                                          
  // Btagging scale factor and uncertainties                                                                                                                                                          
  //########################################                                                                                                                                                          
  TH2F ** Btagg_TT=FuncHistBTagSF();
  


  //########################################
  // Pileup files
  //########################################
  TH1F *  HistoPUData =HistPUData();
  // Need a fix for PU distribution
  //    TH1F * HistoPUMC=HistPUMC();
  //        TH1F *  HistoPUMC =HistPUData();
  size_t isDataXXX;// = InputROOT.find("Data");
  bool check_data=0;
  if (isDataXXX!= string::npos)  check_data=1;
  //TH1F * HistoPUMC=HistPUMC(check_data,f_Double);




  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;

  Long64_t nbytes = 0, nb = 0;
  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      auto thereisSquark= false;
      for  (auto ij=0 ; ij < mcPID->size(); ij++){
        // if (mcStatus->at(ij)==23 || fabs(mcPID->at(ij))==13)     cout << mcPID->at(ij) <<"  "<< mcStatus->at(ij) <<" " << mcPt->at(ij) <<"  "<< mcMomPID->at(ij)<<  " " << mcGMomPID->at(ij)<<"\n";    
        if (fabs(mcPID->at(ij))==3  && mcStatus->at(ij)==23 ) {numOf_s_quark++; thereisSquark= true;}
        //                            cout <<" pdgid is "<< mcPID->at(ij)<<"\n";                                                                                                                          
         if (fabs(mcPID->at(ij))==4  && mcStatus->at(ij)==23 ) numOf_c_quark++;
        //                            cout <<" pdgid is "<< mcPID->at(ij)<<"\n";                                                                                                                          
      }
      if (thereisSquark) continue;
      //      auto Squark= false;                                                                                                                                                                          
      float jetES[3]={-1,0,1};
      std::string ResolJet_Cat[3] = {"JetERDown", "", "JetERUp"};
      std::string ScaleJet_Cat[3] = {"JetESDown", "", "JetESUp"};
      std::string ScaleMETUE_Cat[5] = {"METUESDown", "", "METUESUp","METJESDown","METJESUp"};




      //      std::cout<<"coming here 1b"<<std::endl;                                                                                                                                                     
      //if (Squark) continue;                                                                                                                                                                             

      //###############################################################################################                                                                                                   
      //  MET Filters (only on Data)                                                                                                                                                                      
      //###############################################################################################                                                                                                   
      if (isData && (metFilters!=0)) continue;   //FIXME                                                                                                                                                  
      //std::cout<<"coming here 1c"<<std::endl;                                                                                                                                                           


      //###############################################################################################                                                                                                   
      //      Trigger Requirement                                                                                                                                                                         
      //###############################################################################################                                                                                                   
      bool PassTrigger = ((HLTEleMuX >> 21 & 1)==1); //   else if (name.find("HLT_Mu50_v") != string::npos) bitEleMuX = 21;                                                                               
      if (! PassTrigger) continue;
      //      std::cout<<"coming here 1d"<<std::endl; 

      //###############################################################################################                                                                                                   
      //  This part is to avoid of the duplicate of mu-j pair from one events                                                                                                                             
      //###############################################################################################                                                                                                   
      std::vector<string> HistNamesFilled;
      HistNamesFilled.clear();


      //###############################################################################################                                                                                                   
      //  GenInfo                                                                                                                                                                                         
      //###############################################################################################                                                                                                   
      vector<float>  genInfo=GeneratorInfo();                    

                                                                                                                                     
      //            //######################## Top Pt Reweighting
      float TopPtReweighting = 1;
      size_t isTTJets; //= InputROOT.find("TTJets");
      if (isTTJets!= string::npos) TopPtReweighting = genInfo[0];

            
      //            //######################## W K-factor
      float WBosonPt=0;
      float WBosonMass=0;
      float WBosonKFactor=1;
            
      WBosonPt=genInfo[1];
      WBosonMass=genInfo[3];
            
      size_t isWJetsToLNu_Inc; //= InputROOT.find("WJetsToLNu_Inc");
      size_t isWJets; //= InputROOT.find("JetsToLNu");
      size_t isWToMuNu; //= (InputROOT.find("WToMuNu") );
      size_t isWToTauNu; // = (InputROOT.find("WToTauNu") );
            
      if (WBosonMass > 100 && (isWToMuNu!= string::npos || isWToTauNu!=string::npos)) WBosonKFactor=HistMassDepKFactor->GetBinContent(int(WBosonMass)/10 +1); //Mass binned K-factor
      if (WBosonMass <= 100 && isWJets!= string::npos  )WBosonKFactor= FuncBosonKFactor("W1Cen") + FuncBosonKFactor("W2Cen") * WBosonPt; //HT binned & inclusive K-factor
            
      //................................................................................................................
      //................................................................................................................
      if (isWJets!= string::npos && WBosonMass > 100) continue;
      if (isWJetsToLNu_Inc!= string::npos && (genHT > 100.0)) continue;
      //................................................................................................................
      //................................................................................................................
     
      //            //######################## Z K-factor
      float ZBosonPt=0;
      float ZBosonKFactor=1;
      size_t isDYJets; // = InputROOT.find("DYJets");
      ZBosonPt=genInfo[2];
      if (isDYJets!= string::npos) ZBosonKFactor= FuncBosonKFactor("Z1Cen") + FuncBosonKFactor("Z2Cen") * ZBosonPt;


      //###############################################################################################
      //  Lumi, GEN & PileUp Weight
      //###############################################################################################
            
      float LumiWeight = 1;
      float GetGenWeight=1;
      float PUWeight = 1;

      if (!isData){
                
	//######################## Lumi Weight
	//if (HistoTot) LumiWeight = weightCalc(HistoTot, InputROOT,genHT, W_HTBinROOTFiles, WBosonMass, WMuNu_MassBinROOTFiles,WTauNu_MassBinROOTFiles);
	//######################## Gen Weight
	GetGenWeight = genWeight;
	//std::cout<<"coming till lumiweight "<<std::endl;
	//######################## PileUp Weight
	//                int puNUmmc=int(puTrue->at(0)*10);

	int puNUmmc = int(puTrue->at(0)*5);
	//std::cout<<"coming till punummc "<<std::endl;
	int puNUmdata = int(puTrue->at(0)*5);
	//std::cout<<"coming till pummdata "<<std::endl;
	float PUMC_; //HistoPUMC->GetBinContent(puNUmmc+1);
	//std::cout<<"coming till pumc_ "<<PUMC_<<std::endl;
	float PUData_ = HistoPUData->GetBinContent(puNUmdata+1);
	PUMC_ = PUData_; //HistoPUMC->GetBinContent(puNUmmc+1);
	
	//std::cout<<"PUData_"<<PUData_<<"PUMC_"<<PUMC_<<std::endl;
	//if (PUMC_ == 0)
	//std::cout<<"coming till pudata inside if "<<std::endl;
	//cout<<"PUMC_ is zero!!! & num pileup= "<< puTrue->at(0)<<"\n";
	//else
	//PUWeight = PUData_/PUMC_;
	PUWeight = 1;
	//std::cout<<"coming till pudata inside if "<<std::endl;
	//if  (puNUmmc == 1 || puNUmmc == 0) continue;
	
      }

    
      //std::cout<<"coming after loop "<<std::endl;

      //############################################################################################                                                                                                      
      //   Final Total Weight                                                                                                                                                                             
      //############################################################################################         
      float TotalWeight_withTopPtRW = LumiWeight * GetGenWeight * PUWeight * TopPtReweighting * WBosonKFactor * ZBosonKFactor ;
      float TotalWeight_NoTopPtRW = LumiWeight * GetGenWeight * PUWeight * WBosonKFactor * ZBosonKFactor ;
      //###########       numTau   ###########################################################
      int numTau= getNumTau();
            
      //###########       Ele Veto   ###########################################################
      int numElectron= getNumElectron();
      float ElectronCor= 1;
      //      std::cout<<"coming ele veto "<<std::endl;
            
      //###########       BTag SF   ###########################################################
      float FinalBTagSF=FuncFinalBTagSF(isData,Btagg_TT,BJetPtCut,CSVCut);
            
      //###########       numBJet   ###########################################################
      int numBJet = getnumBJets(BJetPtCut,CSVCut);
      //          std::cout<<"coming numbjet "<<std::endl;
      //###########       numJet   ###########################################################
      int numJet = getnumJets(SimpleJetPtCut);
            
      //###########       numZboson   ###########################################################
      int numZboson = getNumZBoson();


      //###############################################################################################
      //  Some Histogram Filling
      //###############################################################################################
      plotFill("_WeightLumi",LumiWeight,10000,0,1000);
      plotFill("_WeightGen",GetGenWeight,10000,0,1000);
      plotFill("_WeightPU",PUWeight,10000,0,1000);
      plotFill("_WeightTopPtReweighting",TopPtReweighting,100,0,2);
      plotFill("_WeightWBosonKFactor",WBosonKFactor,500,0,5);
      plotFill("_ZeightWBosonKFactor",ZBosonKFactor,500,0,5);
            
      plotFill("_TotalWeight_withTopPtRW",TotalWeight_withTopPtRW,1000,0,100);
      plotFill("_TotalWeight_NoTopPtRW",TotalWeight_NoTopPtRW,1000,0,100);
      plotFill("_nVtx_NoPUCorr",nVtx,60,0,60);
      plotFill("_nVtx_PUCorr",nVtx,60,0,60,PUWeight);
      plotFill("_WBosonPt",WBosonPt,150,0,1500,PUWeight);
      plotFill("_FinalBTagSF", FinalBTagSF,200,0,2);
      //      std::cout<<"coming numzboson "<<std::endl;
      for (int qq=0; qq < 60;qq++){
	if ((HLTEleMuX >> qq & 1) == 1)
	  plotFill("_HLT",qq,60,0,60);
      }





       tree->Fill();

      if (jentry%reportEvery == 0)
	{
	  //float timeLeft = ((sw.RealTime() / 60.0)/(jentry+1)) * (nentriesToCheck-jentry-1);                                                                                                     
	  std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
	}
    }
  if((nentriesToCheck-1)%reportEvery != 0)
    std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop(); // stop stopwatch                                                                                                                                                                         
  std::cout<<"Total number of events: "<<nTotal<<std::endl;
  //fileName = new TFile(file2, "RECREATE")
  std::cout<<"coming here 1a "<<std::endl;

    fileName = new TFile(File2, "RECREATE");
    fileName->cd();
  std::cout<<"coming here 1b "<<std::endl;
  map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
  map<string, TH1F*>::const_iterator jMap1 = myMap1->end();

  for (; iMap1 != jMap1; ++iMap1)
    nplot1(iMap1->first)->Write();

  map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
  map<string, TH2F*>::const_iterator jMap2 = myMap2->end();

  for (; iMap2 != jMap2; ++iMap2)
    nplot2(iMap2->first)->Write();
  fileName->Close();
  cout<< "numOf_c_quark " << numOf_c_quark << "  numOf_s_quark " << numOf_s_quark<<"\n";
}
/* 
void MuonClass::BookHistos(const char* file2)
		     {
  myMap1 = new std::map<std::string, TH1F*>();
  myMap2 = new map<string, TH2F*>();
  fileName = new TFile(file2, "RECREATE");
  fileName->cd();
  tree = new TTree("ADD","ADD");
  //tree->Branch("event_","std::vector<unsigned int>",&event_);
  //tree->Branch("event_info","std::vector<double>",&event_info);
  //    fileName->cd();
  

  std::cout<<"coming here 1b "<<std::endl;        
  map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
  map<string, TH1F*>::const_iterator jMap1 = myMap1->end();

  for (; iMap1 != jMap1; ++iMap1)
    nplot1(iMap1->first)->Write();

  map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
  map<string, TH2F*>::const_iterator jMap2 = myMap2->end();
  
  for (; iMap2 != jMap2; ++iMap2)
    nplot2(iMap2->first)->Write();
  //  fileName->Close(); 		         			 
			 

		     }

*/
