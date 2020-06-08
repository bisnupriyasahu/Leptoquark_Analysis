
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  This code is to pre-select the events. This ise used as a base for background estimation, ...
////////////////////////////////////////////////////////////////////






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

  std::cout<<"comming here 2feb"<<std::endl;
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


  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  TLorentzVector Mu4Momentum, Jet4Momentum,LQ4Momentum,Mu24Momentum;
  Mu4Momentum.Clear();
  Jet4Momentum.Clear();
  LQ4Momentum.Clear();
  Mu24Momentum.Clear();
  Long64_t nentries = fChain->GetEntriesFast();
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
    {
      //      auto Squark= false;
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
	std::cout<<"coming here 1c"<<std::endl; 
      //###########       Trigger Requirement ###########################################################
      bool PassTrigger = ((HLTEleMuX >> 21 & 1)==1); //   else if (name.find("HLT_Mu50_v") != string::npos) bitEleMuX = 21;
      if (! PassTrigger) continue;


	std::cout<<"coming here 1d"<<std::endl; 
	
      //############################################################################################
      //   Final Total Weight
      //############################################################################################

	    //###########       numBJet   ###########################################################
      int numBJet=getNumBJets(BJetPtCut,CSVCut);
      	std::cout<<"coming here 2a"<<std::endl; 
	//###########       numJet   ###########################################################
      int numJet=getNumJets(SimpleJetPtCut);
            	std::cout<<"coming here 2b"<<std::endl; 
      //###########       numZboson   ###########################################################
		int numZboson = getNumZBoson();
		 
		
	std::cout<<"coming here 2d"<<std::endl; 
      //############################################################################################
      //###########       Loop over MuJet events   #################################################
      //############################################################################################
     


      //      Muon_nMu->Fill(nMu);
      
      for  (int imu=0 ; imu < nMu; imu++){
	     float IsoMu=muPFChIso->at(imu)/muPt->at(imu);
	     if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
	       IsoMu= ( muPFChIso->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
	     
	     bool MuPtCut = muPt->at(imu) > LeptonPtCut_ && fabs(muEta->at(imu)) < 2.4 ;
	     bool MuId=( (muIDbit->at(imu) >> 2 & 1)  && fabs(muD0->at(imu)) < 0.045 && fabs(muDz->at(imu)) < 0.2); //Tight Muon Id
	     std::cout<<"coming here 3"<<std::endl; 
	     
	     
	     if (! MuPtCut || !MuId ) continue;
	     
	     //float MuonCor=getCorrFactorMuon94X(isData,  muPt->at(imu), muEta->at(imu) , HistoMuId,HistoMuIso,HistoMuTrg,HistoMuTrack);
	     
	     Mu4Momentum.SetPtEtaPhiM(muPt->at(imu),muEta->at(imu),muPhi->at(imu),MuMass);
	     
	     //	Muon_IsoMu->Fill(IsoMu.at(imu));
	     //  Muon_MuPtCut->Fill(MuPtCut.at(imu));
	     //	Muon_MuId->Fill(MuId.at(imu));
	     //Muon_MuPt->Fill(muPt->at(imu));
	     //Muon_MuEta->Fill(muEta->at(imu));
	     //Muon_MuPhi->Fill(muPhi->at(imu));
	     
	     
	     //###########   compute ST and RecoHT    ######################################################
	     
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
	     	     
	       for (int ijet= 0 ; ijet < nJet ; ijet++){
	       Jet4Momentum.SetPtEtaPhiE(jetPt->at(ijet),jetEta->at(ijet),jetPhi->at(ijet),jetE->at(ijet));
	       std::cout<<"her it comes jet4mom"<<std::endl;
	       //	  bool goodJet = ( jetPt->at(ijet) > JetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5);
	       // bool goodJet = (jetPFLooseId->at(ijet) > 0.5 && jetPt->at(ijet) > JetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5)
	       bool goodJet = ((*jetID)[ijet]>>0&1 == 1 && jetPt->at(ijet) > JetPtCut && fabs(jetEta->at(ijet)) < 2.4 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5);
	       if (! goodJet) continue;
	       LQ4Momentum=Jet4Momentum + Mu4Momentum;
	       std::cout<<"her it comes LQ"<<std::endl;
	       
	       
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
                    

	       }//for jet loop in jet4momentum
	     
      }//for loop imu
      
    }//jentry

      if((nentriesToCheck-1)%reportEvery != 0)
	std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
      sw.Stop();
}//class






void MuonClass::Histograms(const char* file2)
{
  std::cout<<"coming here 10"<<std::endl; 
  fileName = new TFile(file2, "RECREATE");
  fileName->cd();   
    //////////histograms filling/////////////   
    Muon_nMu=new TH1F("Muon_nMu","",10,0,10);
    Muon_MuPt=new TH1F("Muon_MuPt","",100,0,1000);
    Muon_MuEta=new TH1F("Muon_MuEta","",100,-5,5);
    Muon_MuPhi=new TH1F("Muon_MuPhi","",50,-5,5);
    MUJ_ST=new TH1F("MUJ_ST","",100,0,1000);
    MUJ_recoHT=new TH1F("MUJ_recoHT","",100,0,1000);
    Jet_JetEta=new TH1F("Jet_JetEta","",100,-5,5);
    Jet_JetPhi=new TH1F("Jet_JetPhi","",50,-5,5);
}
