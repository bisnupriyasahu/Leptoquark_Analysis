#define LQMETClass_cxx
#include "LQMETClass.h"
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
#include <map>

using namespace std;
std::vector<string> input;



int main(int argc, const char* argv[])
{

  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
      return 1;
    }

  std::cout<<"comming here 2feb"<<std::endl;
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
      return 1;
    }
  myMap1 = new std::map<std::string, TH1F*>();
  myMap2 = new map<string, TH2F*>();
  

    



std::cout<<"hello to input to addinputroot"<<std::endl;
    //for (int f = 4; f < argc; f++) {

    //cout <<"INPUT NAME IS:   " << input[argv[2]] << "\n";
  
  
  LQMETClass t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void LQMETClass::Loop(Long64_t maxEvents, int reportEvery)
{


  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  //  TLorentzVector Mu4Momentum, Jet4Momentum,LQ4Momentum,Mu24Momentum;
  //  TLorentzVector Mu4Momentum,Jet4MomentumNonSmear, Jet4Momentum,KJet4Momentum,NewJet4Collection,LQ;



  //########################################
  // Muon Id, Iso, Trigger and Tracker Eff files
  //########################################
  TH2F** HistoMuId=FuncHistMuId();
   TH2F** HistoMuIso=FuncHistMuIso();
  TH1F** HistoMuTrg=FuncHistMuTrigger();
  TGraphAsymmErrors * HistoMuTrack=FuncHistMuTrack();





  /*
  Mu4Momentum.Clear();
  Jet4Momentum.Clear();
  LQ4Momentum.Clear();
  Mu24Momentum.Clear();
*/  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "Total entries: " << nentries << std::endl;
  Long64_t nentriesToCheck = nentries;

  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;

  Long64_t nbytes = 0, nb = 0;
  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();


  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {//      auto Squark= false;                                         
      float jetES[3]={-1,0,1};
      std::string ResolJet_Cat[3] = {"JetERDown", "", "JetERUp"};
      std::string ScaleJet_Cat[3] = {"JetESDown", "", "JetESUp"};
      std::string ScaleMETUE_Cat[5] = {"METUESDown", "", "METUESUp","METJESDown","METJESUp"};
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //if (ientry % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", ientry, nentries);                                                                                                  
      //fflush(stdout);                                                                                                                                                                                   
      /* for  (auto k=0 ; k < mcPID->size(); k++){                                                                                                                                                        
         if (fabs(mcPID->at(k))==3  && mcStatus->at(k)==23  && fabs(mcMomPID->at(k))== 9000009 ) {numOf_s_quark++; Squark= true;}                                                                         
         if (fabs(mcPID->at(k))==4  && mcStatus->at(k)==23  && fabs(mcMomPID->at(k))== 9000008 ) numOf_c_quark++;                                                                                         
     } //for k                                                                                                                                                                                            


      */
      std::cout<<"coming here 1b"<<std::endl;
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
      std::cout<<"coming here 1d"<<std::endl;


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
      // size_t isTTJets = InputROOT.find("TTJets");
      //if (isTTJets!= string::npos) TopPtReweighting = genInfo[0];

      //            //######################## W K-factor
      float WBosonPt=0;
      float WBosonMass=0;
      float WBosonKFactor=1;
      float WBosonKFactor_ewkUp=1;
      float WBosonKFactor_ewkDown=1;
            
            
      WBosonPt=genInfo[1];
      WBosonMass=genInfo[3];
      /*
      size_t isWJetsToLNu_Inc = InputROOT.find("WJetsToLNu_Inc");
      size_t isWJets = InputROOT.find("WJets");
      size_t isWToMuNu = (InputROOT.find("WToMuNu") );
      size_t isWToTauNu = (InputROOT.find("WToTauNu") );
      */



      //###############################################################################################
      //  Lumi, GEN & PileUp Weight
      //###############################################################################################
            
      float LumiWeight = 1;
      if (!isData){
                
	//######################## Lumi Weight
	//if (HistoTot) LumiWeight = weightCalc(HistoTot, InputROOT,genHT, W_HTBinROOTFiles, WBosonMass, WMuNu_MassBinROOTFiles,WTauNu_MassBinROOTFiles);


      }

      //############################################################################################                                                                                                      
      //   Final Total Weight                                                                                                                                                                             
      //############################################################################################         




      //###########       numTau   ###########################################################
      //int numTau= getNumTau();



      //###########       numBJet   ###########################################################                                                                                                      
      int numBJet=getNumBJets(BJetPtCut,CSVCut);
      std::cout<<"coming here 2a"<<std::endl;


      //###########       numJet   ###########################################################                                                                                                           
      int numJet=getNumJets(SimpleJetPtCut);
      std::cout<<"coming here 2b"<<std::endl;


      //###########       numZboson   ###########################################################                                                                                                         
      int numZboson = getNumZBoson();
      std::cout<<"coming here 2d"<<std::endl;


      //#######################################################################################
      //  Some Histogram Filling
      //#######################################################################################
      //      plotFill("_WeightLumi",LumiWeight,1000,0,10);








      //############################################################################################                                                                                                      
      //###########       Loop over MuJet events   #################################################                                                                                                      
      //###############################################################################################                                                                                                  

      for  (int imu=0 ; imu < nMu; imu++){
                
	float IsoMu=muPFChIso->at(imu)/muPt->at(imu);
	if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
	  IsoMu= ( muPFChIso->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
	//                IsoMu= ( muPFChIso->at(imu)/muPt->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
                
	bool MuPtCut = muPt->at(imu) > LeptonPtCut_ && fabs(muEta->at(imu)) < 2.4 ;
	bool MuIdIso=( (muIDbit->at(imu) >> 2 & 1)  && fabs(muD0->at(imu)) < 0.045 && fabs(muDz->at(imu)) < 0.2); //Tight Muon Id
                
	if (! MuPtCut || !MuIdIso ) continue;
                
                
	float LepCor=getCorrFactorMuon94X(isData,  muPt->at(imu), muEta->at(imu) , HistoMuId,HistoMuIso,HistoMuTrg,HistoMuTrack);
 
	TLorentzVector Mu4Momentum,Jet4MomentumNonSmear, Jet4Momentum,KJet4Momentum,NewJet4Collection,LQ;
	Mu4Momentum.SetPtEtaPhiM(muPt->at(imu),muEta->at(imu),muPhi->at(imu),MuMass);
                
	//mu_Mu4Momentum->Fill(Mu4Momentum);
	//	Muon_IsoMu->Fill(IsoMu);                                                                                                                                                           
	// Muon_MuPtCut->Fill(MuPtCut.at(imu));                                                                                                                                                      
	// Muon_MuId->Fill(MuId.at(imu));                                                                                                                                                             
	//Muon_MuPt->Fill(muPt->at(imu));                                                                                                                                                             
	//Muon_MuEta->Fill(muEta->at(imu));                                                                                                                                                           
	//Muon_MuPhi->Fill(muPhi->at(imu));                                                                                                                                                           


	/*	//###########   compute ST and RecoHT    ######################################################                                                                                               

	float recoHT=0;

	for (int ijet= 0 ; ijet < nJet ; ijet++){
	  // if (jetPFLooseId->at(ijet) > 0.5 && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5)                                          
	  if ((*jetID)[ijet]>>0&1 == 1 && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5)
	    recoHT += jetPt->at(ijet);
	
	  std::cout<<"recoHT   = "<<recoHT<<std::endl;
	}

	float ST=recoHT+muPt->at(imu);
	std::cout<<"ST  = "<<ST<<std::endl;
	

	//###########    loop over  Jet    #########################################                                                                                                                  
	*/
	for (int ijet= 0 ; ijet < nJet ; ijet++){
	 
    
	  float JetSmearResolution[3]={1,1,1};
	  if (!isData){
	    // NOT available for all samples need to apply for signal though
	    //                        JetSmearResolution[0]=jetP4SmearDo->at(ijet);
	    //                        JetSmearResolution[1]=jetP4Smear->at(ijet);
	    //                        JetSmearResolution[2]=jetP4SmearUp->at(ijet);
	    //
	    JetSmearResolution[0]=1.0;
	    JetSmearResolution[1]=1.0;
	    JetSmearResolution[2]=1.0;
                        
	  }
	  float UESMET[5]={pfMET_T1UESDo,pfMET,pfMET_T1UESUp,pfMET_T1JESDo,pfMET_T1JESUp};
	  float UESMETPhi[5]={pfMETPhi_T1UESDo,pfMETPhi,pfMETPhi_T1UESUp,pfMETPhi_T1JESDo,pfMETPhi_T1JESUp};


	  for (int jetRes=0;jetRes<3;jetRes++){
	    for (int metUE=0; metUE < 5; metUE++){
	      for (int jetScl=0;jetScl<3;jetScl++){
                                    
                                    
		// This is to check that we only make the plots only either JES or MET is applied (not both of them simultaneously)!
		if (jetRes!=1 && metUE != 1  ) continue;
		if (jetRes!=1 && jetScl != 1 ) continue;
		if (metUE != 1 && jetScl != 1 ) continue;
                                    
                                    
		Jet4MomentumNonSmear.SetPtEtaPhiE(jetPt->at(ijet),jetEta->at(ijet),jetPhi->at(ijet),jetE->at(ijet));
		Jet4Momentum=Jet4MomentumNonSmear*JetSmearResolution[jetRes];
                                    
                                    
		NewJet4Collection.SetPtEtaPhiE(Jet4Momentum.Pt()*(1+ jetES[jetScl]*jetJECUnc->at(ijet)) ,jetEta->at(ijet),jetPhi->at(ijet),Jet4Momentum.E()*(1+ jetES[jetScl]*jetJECUnc->at(ijet)));
                                    
                                    
                                    
		bool goodJet = ((*jetID)[ijet]>>0&1 == 1 && NewJet4Collection.Pt() > JetPtCut && fabs(NewJet4Collection.Eta() ) < 2.4 && NewJet4Collection.DeltaR(Mu4Momentum) > 0.5);
		if (! goodJet) continue;
		float jetMET=UESMET[metUE];
		float jetMETPhi=UESMETPhi[metUE];
                                    
		float jetMET_x = UESMET[metUE] * TMath::Cos(UESMETPhi[metUE]) - (Jet4Momentum.Px()- NewJet4Collection.Px()) ;
		float jetMET_y = UESMET[metUE] * TMath::Sin(UESMETPhi[metUE]) - (Jet4Momentum.Py()- NewJet4Collection.Py()) ;
		jetMET = sqrt (pow(jetMET_x,2)+ pow(jetMET_y,2));
		jetMETPhi = atan(jetMET_y / jetMET_x);
		if (UESMETPhi[metUE] > (TMath::Pi() / 2)) jetMETPhi += TMath::Pi();
		if (UESMETPhi[metUE] < (-TMath::Pi() / 2)) jetMETPhi -= TMath::Pi();
                                    
                                    
		std::cout<<"her it comes LQ"<<std::endl;                                    
		LQ=NewJet4Collection + Mu4Momentum;



		//###############################################################################################
		// Apply all veto cuts
		//###############################################################################################
		//                                    if ((numTau+numElectron +numZboson + numBJet) > 0) continue;
		//###############################################################################################
		//  dPhi Jet_MET Categorization (To get a relaxed cut for QCD template)
		//###############################################################################################
		const int size_jetMetPhi = 2;
                                    
		bool HighDPhi = deltaPhi(NewJet4Collection.Phi(),jetMETPhi) > 0.5 && deltaPhi(Mu4Momentum.Phi(),jetMETPhi) > 0.5;
		bool noDPhi = 1;
                                    
                                    
		bool jetMetPhi_category[size_jetMetPhi] = {HighDPhi,noDPhi};
		std::string jetMetPhi_Cat[size_jetMetPhi] = {"", "_NoDPhi"};
                                    


		//###############################################################################################
		//  Isolation Categorization
		//###############################################################################################
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
		float tmass_LQMet= TMass_F(LQ.Pt(), LQ.Px(),LQ.Py(), pfMET, pfMETPhi);
		float tmass_MuMet= TMass_F(muPt->at(imu), muPt->at(imu)*cos(muPhi->at(imu)),muPt->at(imu)*sin(muPhi->at(imu)) , jetMET, jetMETPhi);
		const int size_mTCat = 3;
		bool MT100 = tmass_MuMet > 100;
		bool MT50To150=(tmass_MuMet > 50 && tmass_MuMet <= 150);
		bool MT400 = tmass_MuMet > 400;
		bool MT_category[size_mTCat] = {MT100,MT50To150,MT400};
		std::string MT_Cat[size_mTCat] = {"_MT100","_MT50To150","_MT400"};


		std::string CHL="MuJet";
                                    
		plotFill("Weight_Mu", LepCor,200,0,2);

	      }
	    }
	  }



	}//for jet loop in jet4momentum                                                                                                                                                             

      }//for loop imu                                                                                                                                                                                     

    }//jentry                                                                                                                                                                                             
 
  if((nentriesToCheck-1)%reportEvery != 0)
    std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop();
  //  fileName = new TFile(file2, "RECREATE");
  fileName->cd();
  map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
  map<string, TH1F*>::const_iterator jMap1 = myMap1->end();

  for (; iMap1 != jMap1; ++iMap1)
    nplot1(iMap1->first)->Write();
  map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
  map<string, TH2F*>::const_iterator jMap2 = myMap2->end();

  for (; iMap2 != jMap2; ++iMap2)
    nplot2(iMap2->first)->Write();

  fileName->Close();



}//class 


void LQMETClass::Histograms(const char* file2)
{
  std::cout<<"coming here 10"<<std::endl;
  //   fileName = new TFile(file2, "RECREATE");
  
//    fileName->cd();
  /*
  map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
  map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
    
  for (; iMap1 != jMap1; ++iMap1)
    nplot1(iMap1->first)->Write();
  map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
  map<string, TH2F*>::const_iterator jMap2 = myMap2->end();
    
  for (; iMap2 != jMap2; ++iMap2)
    nplot2(iMap2->first)->Write();

*/
  //  fileName->Close();

  //////////histograms filling/////////////                                                                                                                                                              
  //  Muon_IsoMu=new TH1F("Muon_IsoMu","muiso",10,0,10);
  /*  Muon_MuPt=new TH1F("Muon_MuPt","",100,0,1000);
  Muon_MuEta=new TH1F("Muon_MuEta","",100,-5,5);
  Muon_MuPhi=new TH1F("Muon_MuPhi","",50,-5,5);
    
 mu_Mu4Momentum=new TH1F("mu_Mu4Momentum",100,0,1000);
  *///  MUJ_ST=new TH1F("MUJ_ST","",100,0,1000);
  // MUJ_recoHT=new TH1F("MUJ_recoHT","",100,0,1000);
  /*Jet_JetEta=new TH1F("Jet_JetEta","",100,-5,5);
  Jet_JetPhi=new TH1F("Jet_JetPhi","",50,-5,5);
  */
}
