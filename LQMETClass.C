#define LQMETClass_cxx
//#include "LQMETClass.h"
#include "function.h"
#include <string>
#include <ostream>
#include <vector>

int main(int argc, const char* argv[])
{

  auto numOf_c_quark=0;
  auto numOf_s_quark=0;  
  TH1F *_MuPt;
  _MuPt = new TH1F("_MuPt","",100,0,200);  
  std::string out = *(argv + 1);
  
  cout << "\n\n\n OUTPUT NAME IS:    " << out << endl;     //PRINTING THE OUTPUT name
  TFile *fout = TFile::Open(out.c_str(), "RECREATE");
  
  myMap1 = new std::map<std::string, TH1F*>();
  myMap2 = new map<string, TH2F*>();
  
  std::vector<string> input;
  for (int f = 3; f < argc; f++) {
    input.push_back(*(argv + f));
    
    //    cout <<"INPUT NAME IS:   " << input[f - 3] << "\n";
  }
  
  
  
    
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
  TH2F * HistoEleMVAIdIso90= FuncHistEleMVAId("Tot");
  TH2F * HistoEleMVAIdIso90_EffMC= FuncHistEleMVAId("MC");
  TH2F * HistoEleMVAIdIso90_EffData= FuncHistEleMVAId("Data");
 
 

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



  //########################################################################################################################################################
  //########################################################################################################################################################
  //########################################################################################################################################################
  //                                                  Loop over inout ROOT files
  //########################################################################################################################################################
  //########################################################################################################################################################
  //########################################################################################################################################################
  for (int k = 0; k < input.size(); k++) {
    
    TFile *f_Double = TFile::Open(input[k].c_str());
    
    //    cout << "\n  Now is running on ------->   " << std::string(f_Double->GetName()) << "\n";
    
    std::string InputROOT= std::string(f_Double->GetName());
    TFile * myFile = TFile::Open(f_Double->GetName());
    TH1F * HistoTot = (TH1F*) myFile->Get("hcount");
    //std::cout<<"HistoTot entries  "<<HistoTot->GetEntries()<<std::endl;
    TTree *  Run_Tree;
    //    std::cout<<"coming tiil here"<<std::endl;        
    Run_Tree= Xttree(f_Double);
    
    
    //########################################
    // Pileup files
    //########################################
  
    TH1F *  HistoPUData =HistPUData();
    // Need a fix for PU distribution
    
    //    TH1F * HistoPUMC=HistPUMC();
    //        TH1F *  HistoPUMC =HistPUData();
    size_t isDataXXX = InputROOT.find("Data");
    bool check_data=0;
    if (isDataXXX!= string::npos)  check_data=1;
    //    TH1F * HistoPUMC=HistPUMC(check_data,f_Double);
    
    //#####################################################################
    //#####################################################################
    //                           Loop over Events in each ROOT files
    //#####################################################################
    //#####################################################################
    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntriesFast();
    //   cout<<"nentries_wtn===="<<nentries_wtn<<"\n";
    
       for (Int_t i = 0; i < nentries_wtn; i++) {
      //                    for (Int_t i = 0; i < 10000; i++) {
      auto thereisSquark= false;
      Run_Tree->GetEntry(i);
      if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
      fflush(stdout);
      for  (auto ij=0 ; ij < mcPID->size(); ij++){
	// if (mcStatus->at(ij)==23 || fabs(mcPID->at(ij))==13)     cout << mcPID->at(ij) <<"  "<< mcStatus->at(ij) <<" " << mcPt->at(ij) <<"  "<< mcMomPID->at(ij)<<  " " << mcGMomPID->at(ij)<<"\n";
	if (fabs(mcPID->at(ij))==3  && mcStatus->at(ij)==23 ) {numOf_s_quark++; thereisSquark= true;}
	//                            cout <<" pdgid is "<< mcPID->at(ij)<<"\n";
	if (fabs(mcPID->at(ij))==4  && mcStatus->at(ij)==23 ) numOf_c_quark++;
	//                            cout <<" pdgid is "<< mcPID->at(ij)<<"\n";
      }
       if (thereisSquark) continue;
      
       
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
       size_t isTTJets = InputROOT.find("TTJets");
      if (isTTJets!= string::npos) TopPtReweighting = genInfo[0];

            
      //            //######################## W K-factor
      float WBosonPt=0;
      float WBosonMass=0;
      float WBosonKFactor=1;
            
      WBosonPt=genInfo[1];
      WBosonMass=genInfo[3];
            
      size_t isWJetsToLNu_Inc = InputROOT.find("WJetsToLNu_Inc");
      size_t isWJets = InputROOT.find("JetsToLNu");
      size_t isWToMuNu = (InputROOT.find("WToMuNu") );
      size_t isWToTauNu = (InputROOT.find("WToTauNu") );
            
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
      size_t isDYJets = InputROOT.find("DYJets");
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
	//	if (HistoTot) LumiWeight = weightCalc(HistoTot, InputROOT,genHT, W_HTBinROOTFiles, WBosonMass, WMuNu_MassBinROOTFiles,WTauNu_MassBinROOTFiles);
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
	//	if (PUMC_ == 0)
	//	std::cout<<"coming till pudata inside if "<<std::endl;
	//cout<<"PUMC_ is zero!!! & num pileup= "<< puTrue->at(0)<<"\n";
	//else
	//	PUWeight = PUData_/PUMC_;
	PUWeight = 1;
	//std::cout<<"coming till pudata inside if "<<std::endl;
	//	if  (puNUmmc == 1 || puNUmmc == 0) continue;
	
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
      int numBJet = numBJets(BJetPtCut,CSVCut);
      //          std::cout<<"coming numbjet "<<std::endl;
      //###########       numJet   ###########################################################
      int numJet=numJets(SimpleJetPtCut);
            
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
      //      std::cout<<"coming before loop  of neve "<<std::endl;

     

      //############################################################################################                                                                                                      
      //###########       Loop over MuJet events   #################################################                                                                                                      
      //###############################################################################################                                                                                                  

      TLorentzVector Mu4Momentum, Jet4Momentum,KJet4Momentum,LQ4Momentum,Mu24Momentum;
      for  (int imu=0 ; imu < nMu; imu++){
                
	float IsoMu=muPFChIso->at(imu)/muPt->at(imu);
	if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
	  IsoMu= ( muPFChIso->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
                
	bool MuPtCut = muPt->at(imu) > LeptonPtCut_ && fabs(muEta->at(imu)) < 2.4 ;
	bool MuId=( (muIDbit->at(imu) >> 2 & 1)  && fabs(muD0->at(imu)) < 0.045 && fabs(muDz->at(imu)) < 0.2); //Tight Muon Id
                

	if (! MuPtCut || !MuId ) continue;
        //        	std::cout<<"coming inside imuloop 1"<<std::endl;
  
	//###########       numMu   ###########################################################                                                                                                           
	int numMu= getNumMu();
        _MuPt->Fill(muPt->at(imu));

	plotFill("_NumMU", numMu,10,0,10);


	
	float MuonCor=getCorrFactorMuon94X(isData,  muPt->at(imu), muEta->at(imu) , HistoMuId,HistoMuIso,HistoMuTrg,HistoMuTrack);
                
	Mu4Momentum.SetPtEtaPhiM(muPt->at(imu),muEta->at(imu),muPhi->at(imu),MuMass);

	//	std::cout<<"coming inside imuloop 2"<<std::endl;
	//###########   compute ST and RecoHT    ######################################################                                                                                               

	float recoHT=0;

	for (int ijet= 0 ; ijet < nJet ; ijet++){
	  // if (jetPFLooseId->at(ijet) > 0.5 && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5)                                          
	  if (((*jetID)[ijet]>>0&1 == 1) && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5)
	    recoHT += jetPt->at(ijet);
	  //std::cout<<"recoHT   = "<<recoHT<<std::endl;
	}

	float ST=recoHT+muPt->at(imu);
	//std::cout<<"ST  = "<<ST<<std::endl;
	

	//###########    loop over  Jet    #########################################                                                                                                                  
	for (int ijet= 0 ; ijet < nJet ; ijet++){
	  //std::cout<<"coming inside ijet "<<std::endl;
	  //std::cout<<"jetpt at ijetis here "<<jetPt->at(ijet)<<"jeteta is "<<jetEta->at(ijet)<<"jetPhi at ijet"<<jetPhi->at(ijet)<<"jet energy"<<jetE->at(ijet)<<std::endl;
	  Jet4Momentum.SetPtEtaPhiE(jetPt->at(ijet),jetEta->at(ijet),jetPhi->at(ijet),jetE->at(ijet));

	  
	  bool goodJet = (((*jetID)[ijet]>>0&1 == 1) && jetPt->at(ijet) > JetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5);
	  if (! goodJet) continue;
                    
	  LQ4Momentum=Jet4Momentum + Mu4Momentum;

	  //###############################################################################################
	  //  Isolation Categorization
	  //###############################################################################################
	  bool LepPassIsolation= IsoMu < LeptonIsoCut;
                    
	  const int size_isoCat = 2;
	  bool Isolation = LepPassIsolation;
	  bool AntiIsolation =  !LepPassIsolation;
                    
	  bool Iso_category[size_isoCat] = {Isolation, AntiIsolation};
	  std::string iso_Cat[size_isoCat] = {"_Iso", "_AntiIso"};
	  //###############################################################################################
	  //  MT Categorization
	  //###############################################################################################
	  float tmass_MuMet= TMass_F(muPt->at(imu), muPt->at(imu)*cos(muPhi->at(imu)),muPt->at(imu)*sin(muPhi->at(imu)) , pfMET, pfMETPhi);
	  float tmass_JetMet= TMass_F(jetPt->at(ijet), jetPt->at(ijet)*cos(jetPhi->at(ijet)),jetPt->at(ijet)*sin(jetPhi->at(ijet)) , pfMET, pfMETPhi);
	  float tmass_LQMet= TMass_F(LQ4Momentum.Pt(), LQ4Momentum.Px(),LQ4Momentum.Py(), pfMET, pfMETPhi);
                    
	  //        	std::cout<<"coming inside imuloop end "<<std::endl;
	  const int size_mTCat = 5;	  
  
	  bool NoMT = 1;
	  bool HighMT = (tmass_MuMet > 100);
	  bool MT50To150=(tmass_MuMet > 50 && tmass_MuMet <= 150);
	  bool MTMore300=tmass_MuMet > 300 ;
	  bool MTMore500=tmass_MuMet > 500 ;
                    
	  bool MT_category[size_mTCat] = {NoMT,HighMT,MT50To150,MTMore300,MTMore500};
	  std::string MT_Cat[size_mTCat] = {"_NoMT","_HighMT","_MT50To150","_MT300","_MT500"};
                    
	  //###############################################################################################
	  //  dPhi Jet_MET Categorization
	  //###############################################################################################
                    
	  const int size_jetMetPhi = 1;
	  bool HighDPhi = (deltaPhi(Jet4Momentum.Phi(),pfMETPhi) >= 0.5 && deltaPhi(Mu4Momentum.Phi(),pfMETPhi) >= 0.5  );
	  bool jetMetPhi_category[size_jetMetPhi] = {HighDPhi};
	  std::string jetMetPhi_Cat[size_jetMetPhi] = {"_HighDPhi"};
                    

	  //###############################################################################################
	  //  TTbar & DY control region Categorization
	  //###############################################################################################
	  const int size_CR = 3;
                    
	  bool signalRegion = numTau+numZboson + numElectron  + numBJet < 1;
                    
	  bool TTcontrolRegion_DiLep = (numTau + numZboson < 1 && numElectron > 0 );
	  if (TTcontrolRegion_DiLep) {
	    FinalBTagSF=1;
	    ElectronCor= getElectronCor(HistoEleMVAIdIso90);
	  }
                    
	  bool TTcontrolRegion_SingleLep = (numTau+numZboson + numElectron  < 1  && numBJet >= 1);
                    
	  bool region_category[size_CR] = {signalRegion,TTcontrolRegion_DiLep,TTcontrolRegion_SingleLep};
	  std::string region_Cat[size_CR] = {"", "_ttbarCRDiLep","_ttbarCRSingleLep"};
	  //	  /###############################################################################################
	  //  Top Pt Reweighting Cat: The SF is meant to correct only the shape of the pt(top) distribution- not the amount of generated events ( you have to consider that the average weight is not 1 ! ) So we define two category for ttbar events
	  
	  int size_topPtRW =2;
	  
	  float TotalWeight[2] = {TotalWeight_withTopPtRW,TotalWeight_NoTopPtRW};
	  std::string topPtRW[2] = {"", "_NoTopRW"};
                    
	  if (isTTJets == string::npos) size_topPtRW = 1; // If the sample in not ttbar, don't care about new category
                    
	  //###############################################################################################	
 
	  std::string CHL="MuJet";
	        
	  plotFill("NumTau", numTau,10,0,10);
	  plotFill("NumElectron", numElectron,10,0,10);
	  plotFill("NumZboson", numZboson,10,0,10);
	  // plotFill("numElectron", numElectron,10,0,10);
	  plotFill("Weight_Mu", MuonCor,200,0,2);
	  plotFill("Weight_Ele", ElectronCor,200,0,2);
	  plotFill("TotalWeight_Mu",TotalWeight[0]*MuonCor,10000,0,1000);
	  plotFill("TotalNonLumiWeight_Mu",TotalWeight[0]*MuonCor/LumiWeight,10000,0,1000);	   

	   for (int iso = 0; iso < size_isoCat; iso++) {
	    if (Iso_category[iso]) {
	      for (int imt = 0; imt < size_mTCat; imt++) {
		if (MT_category[imt]) {
		  for (int jpt = 0; jpt < size_jetMetPhi; jpt++) {
		    if (jetMetPhi_category[jpt]) {
		      for (int iCR = 0; iCR < size_CR; iCR++) {
			if (region_category[iCR]) {
			  for (int itopRW = 0; itopRW < size_topPtRW; itopRW++) {
                                                        
			    ////////////////////////////////////   Naming the Histogram
			    float FullWeight = TotalWeight[itopRW] * MuonCor * ElectronCor * FinalBTagSF ;
			    std::string FullStringName = topPtRW[itopRW] + MT_Cat[imt] + jetMetPhi_Cat[jpt] +  region_Cat[iCR] + iso_Cat[iso]  ;
                                                        
			    //##################
			    //This check is used to make sure that each event is just filled once for any of the categories ==> No doube-counting of events  (this is specially important for ttbar events where we have many jets and leptons)
			    if (!( std::find(HistNamesFilled.begin(), HistNamesFilled.end(), FullStringName) != HistNamesFilled.end())){
			      HistNamesFilled.push_back(FullStringName);
                                                            
			       if (1){
				//##################
			                                        
			                                 
                                                                
				plotFill(CHL+"_LQMass"+FullStringName,LQ4Momentum.M(),200,0,2000,FullWeight);
				plotFill(CHL+"_LepPt"+FullStringName,muPt->at(imu),2000,0,2000,FullWeight);
				
				plotFill(CHL+"_tmass_MuMet"+FullStringName,tmass_MuMet,200,0,2000,FullWeight);
				plotFill(CHL+"_MET"+FullStringName,pfMET,200,0,2000,FullWeight);
				plotFill(CHL+"_METPhi"+FullStringName,pfMETPhi,400,-4,4,FullWeight);
				plotFill(CHL+"_JetPt"+FullStringName,jetPt->at(ijet) ,2000,0,2000,FullWeight);
				plotFill(CHL+"_JetEta"+FullStringName,jetEta->at(ijet),120,-3,3,FullWeight);
				
				plotFill(CHL+"_LepEta"+FullStringName,muEta->at(imu),100,-2.5,2.5,FullWeight);
				plotFill(CHL+"_tmass_LQMet"+FullStringName,tmass_LQMet,200,0,2000,FullWeight);
				
				plotFill(CHL+"_dPhi_Jet_Met"+FullStringName,deltaPhi(Jet4Momentum.Phi(),pfMETPhi),160,0,3.2,FullWeight);
				plotFill(CHL+"_dPhi_Mu_Met"+FullStringName,deltaPhi(Mu4Momentum.Phi(),pfMETPhi),160,0,3.2,FullWeight);
				plotFill(CHL+"_dPhi_Mu_Jet"+FullStringName,deltaPhi(Mu4Momentum.Phi(),Jet4Momentum.Phi()),160,0,3.2,FullWeight);
				plotFill(CHL+"_NumJet"+FullStringName,numJet,10,0,10,FullWeight);
				plotFill(CHL+"_NumBJet"+FullStringName,numBJet,10,0,10,FullWeight);
				plotFill(CHL+"_WBosonMass"+FullStringName,WBosonMass,200,0,2000,FullWeight);
				plotFill(CHL+"_genHT"+FullStringName,genHT,200,0,2000,FullWeight);
				plotFill(CHL+"_recoHT"+FullStringName,recoHT,300,0,3000,FullWeight);
				plotFill(CHL+"_ST"+FullStringName,recoHT+muPt->at(imu),300,0,3000,FullWeight);
				plotFill(CHL+"_dR_Mu_Jet"+FullStringName,Jet4Momentum.DeltaR(Mu4Momentum),500,0,5,FullWeight);
				plotFill(CHL+"_LepPhi"+FullStringName,muPhi->at(imu),400,-4,4,FullWeight);
				plotFill("_FinalBTagSF"+FullStringName, FinalBTagSF,200,0,2);
				
				plotFill(CHL+"_nVtx"+FullStringName,nVtx,100,0,100,FullWeight);
				if (PUWeight!=0)   plotFill(CHL+"_nVtx_NoPU"+FullStringName,nVtx,100,0,100,FullWeight/ PUWeight);
				plotFill("_FinalBTagSF"+FullStringName, FinalBTagSF,200,0,2);
                                                                
				
				
			       }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	 }
	  
	  //###############################################################################################
	  //  Doing EleTau Analysis
	  //###############################################################################################
                    
	  

	  	}//for jet loop in jet4momentum                                                                                                                                                             
	
      }//for loop imu                                                                                                                                                                                     
      
       }//nentry                                                                                                                                                                                             
    
    //if((nentriesToCheck-1)%reportEvery != 0)
    //std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
    //sw.Stop();
    //  fileName = new TFile(file2, "RECREATE");
     
  }
   
fout->cd();
 _MuPt->Write(); 
  map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
  map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
    
  for (; iMap1 != jMap1; ++iMap1)
    nplot1(iMap1->first)->Write();
    
  map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
  map<string, TH2F*>::const_iterator jMap2 = myMap2->end();
    
  for (; iMap2 != jMap2; ++iMap2)
    nplot2(iMap2->first)->Write();
    
  fout->Close();
  cout<< "numOf_c_quark " << numOf_c_quark << "  numOf_s_quark " << numOf_s_quark<<"\n";
  }//mainfunction
