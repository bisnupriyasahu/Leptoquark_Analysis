#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "math.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <iostream>
#include <map>
#include "TLorentzVector.h"
#include "LQMETClass.h"
#include "corrector.h"
#include "plotfill_histo.h"
#include "Weightclaculator.h"
using namespace std;




float deltaPhi(float a, float b) {
  float result = a - b;
  while (result > M_PI) result -= 2 * M_PI;
  while (result <= -M_PI) result += 2 * M_PI;
  return fabs(result);
}

float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
  return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}


float dR_(float ieta, float iphi, float jeta, float jphi){
    
  float deta=ieta-jeta;
  float dphi=deltaPhi(iphi,jphi);
  return sqrt(pow(deta,2)+pow(dphi,2));
}


TTree *  Xttree( TFile * f_Double){
    
  //            TTree *Run_Tree = (TTree*) f_Double->Get("ggNtuplizer/EventTree");
  TTree *Run_Tree = (TTree*) f_Double->Get("phoJetNtuplizer/eventTree");
  //  std::cout<<"coming tiil here"<<std::endl;  
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(6);
  //std::cout<<"coming tiil here1"<<std::endl;  
  //########################################   General Info
  //std::cout<<"coming tiil here2"<<std::endl;  
  Run_Tree->SetBranchAddress("isData", &isData);
  //std::cout<<"coming tiil here2a"<<std::endl;  
  Run_Tree->SetBranchAddress("run", &run);

  Run_Tree->SetBranchAddress("lumis", &lumis);
  Run_Tree->SetBranchAddress("event", &event);

  Run_Tree->SetBranchAddress("genWeight",&genWeight);
  Run_Tree->SetBranchAddress("HLTEleMuX", &HLTEleMuX);
  Run_Tree->SetBranchAddress("puTrue", &puTrue);
  Run_Tree->SetBranchAddress("nVtx",&nVtx);
  //std::cout<<"coming tiil here2a"<<std::endl;  
  //########################################   MC Info
  Run_Tree->SetBranchAddress("nMC", &nMC);
  Run_Tree->SetBranchAddress("mcPID", &mcPID);
  Run_Tree->SetBranchAddress("mcStatus", &mcStatus);
  Run_Tree->SetBranchAddress("mcPt", &mcPt );
  Run_Tree->SetBranchAddress("mcEta", &mcEta );
  Run_Tree->SetBranchAddress("mcPhi", &mcPhi );
  Run_Tree->SetBranchAddress("mcE", &mcE );
  Run_Tree->SetBranchAddress("mcMass", &mcMass );
  //Run_Tree->SetBranchAddress("mcMomPID", &mcMomPID );
  //Run_Tree->SetBranchAddress("mcGMomPID", &mcGMomPID );
  Run_Tree->SetBranchAddress("mcStatusFlag",&mcStatusFlag);
  //std::cout<<"coming tiil here3"<<std::endl;  
  //########################################   Tau Info
  Run_Tree->SetBranchAddress("nTau", &nTau);
  Run_Tree->SetBranchAddress("tau_Pt"  ,&tau_Pt);
  Run_Tree->SetBranchAddress("tau_Eta"  ,&tau_Eta);
  Run_Tree->SetBranchAddress("tau_Phi"  ,&tau_Phi);
  Run_Tree->SetBranchAddress("tau_Mass"  ,&tau_Mass);
  Run_Tree->SetBranchAddress("tau_Charge"  ,&tau_Charge);
  //  Run_Tree->SetBranchAddress("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding);
  //Run_Tree->SetBranchAddress("tauByTightMuonRejection3", &tauByTightMuonRejection3);
  //Run_Tree->SetBranchAddress("tauByLooseMuonRejection3", &tauByLooseMuonRejection3);
  /* Run_Tree->SetBranchAddress("tauByMVA6MediumElectronRejection"  ,&tauByMVA6MediumElectronRejection);
  Run_Tree->SetBranchAddress("tauByLooseCombinedIsolationDeltaBetaCorr3Hits",&tauByLooseCombinedIsolationDeltaBetaCorr3Hits);
  Run_Tree->SetBranchAddress("tauByMediumCombinedIsolationDeltaBetaCorr3Hits",&tauByMediumCombinedIsolationDeltaBetaCorr3Hits);
  Run_Tree->SetBranchAddress("tauByMVA6LooseElectronRejection", &tauByMVA6LooseElectronRejection);
  */Run_Tree->SetBranchAddress("tau_Dxy",&tau_Dxy);
  Run_Tree->SetBranchAddress("tau_DecayMode",&tau_DecayMode);
  /*  Run_Tree->SetBranchAddress("tauByLooseIsolationMVArun2v1DBoldDMwLT",&tauByLooseIsolationMVArun2v1DBoldDMwLT);
  Run_Tree->SetBranchAddress("tauByVLooseIsolationMVArun2v1DBoldDMwLT",&tauByVLooseIsolationMVArun2v1DBoldDMwLT);
  */Run_Tree->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017" ,&tau_byIsolationMVArun2017v2DBoldDMwLTraw2017);
  //std::cout<<"coming tiil here4"<<std::endl;  
  //########################################   Mu Info
  Run_Tree->SetBranchAddress("nMu", &nMu);
  Run_Tree->SetBranchAddress("muPt"  ,&muPt);
  Run_Tree->SetBranchAddress("muEta"  ,&muEta);
  Run_Tree->SetBranchAddress("muPhi"  ,&muPhi);
  Run_Tree->SetBranchAddress("muIsoTrk", &muIsoTrk);
  Run_Tree->SetBranchAddress("muCharge",&muCharge);
  Run_Tree->SetBranchAddress("muIDbit",&muIDbit);//NEW
  Run_Tree->SetBranchAddress("muPFChIso", &muPFChIso);
  Run_Tree->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
  Run_Tree->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  Run_Tree->SetBranchAddress("muPFPUIso", &muPFPUIso);
  Run_Tree->SetBranchAddress("muD0",&muD0);
  Run_Tree->SetBranchAddress("muDz",&muDz);
  //std::cout<<"coming tiil here5"<<std::endl; 
  //########################################   Ele Info
  Run_Tree->SetBranchAddress("nEle", &nEle);
  Run_Tree->SetBranchAddress("eleCharge", &eleCharge);
  Run_Tree->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent);
  Run_Tree->SetBranchAddress("eleE", &eleE);
  Run_Tree->SetBranchAddress("eleSCE", &eleSCE);
  //Run_Tree->SetBranchAddress("eleESEn", &eleESEn);
  //Run_Tree->SetBranchAddress("eleESEnP1", &eleESEnP1);
  //Run_Tree->SetBranchAddress("eleESEnP2", &eleESEnP2);
  Run_Tree->SetBranchAddress("eleD0", &eleD0);
  Run_Tree->SetBranchAddress("eleDz", &eleDz);
  Run_Tree->SetBranchAddress("eleSIP", &eleSIP);
  Run_Tree->SetBranchAddress("elePt", &elePt);
  Run_Tree->SetBranchAddress("eleEta", &eleEta);
  Run_Tree->SetBranchAddress("elePhi", &elePhi);
  Run_Tree->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5);
  Run_Tree->SetBranchAddress("eleCalibEt", &eleCalibEt);
  Run_Tree->SetBranchAddress("eleCalibE", &eleCalibE);
  Run_Tree->SetBranchAddress("eleSCEta", &eleSCEta);
  Run_Tree->SetBranchAddress("eleSCPhi", &eleSCPhi);
  Run_Tree->SetBranchAddress("eleSCRawE", &eleSCRawE);
  Run_Tree->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth);
  Run_Tree->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth);
  Run_Tree->SetBranchAddress("eleHoverE", &eleHoverE);
  Run_Tree->SetBranchAddress("eleEoverP", &eleEoverP);
  //Run_Tree->SetBranchAddress("eleEoverPout", &eleEoverPout);
  Run_Tree->SetBranchAddress("eleEoverPInv", &eleEoverPInv);
  Run_Tree->SetBranchAddress("eleBrem", &eleBrem);
  Run_Tree->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx);
  Run_Tree->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx);
  //Run_Tree->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo);
  Run_Tree->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5);
  Run_Tree->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5);
  Run_Tree->SetBranchAddress("eleConvVeto", &eleConvVeto);
  Run_Tree->SetBranchAddress("eleMissHits", &eleMissHits);
  //Run_Tree->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR);
  Run_Tree->SetBranchAddress("elePFChIso", &elePFChIso);
  Run_Tree->SetBranchAddress("elePFPhoIso", &elePFPhoIso);
  Run_Tree->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
  //Run_Tree->SetBranchAddress("elePFPUIso", &elePFPUIso);
  //Run_Tree->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso);
  //Run_Tree->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso);
  //Run_Tree->SetBranchAddress("elePFMiniIso", &elePFMiniIso);
  Run_Tree->SetBranchAddress("eleMVAIsoID", &eleMVAIsoID);
  Run_Tree->SetBranchAddress("eleMVAnoIsoID", &eleMVAnoIsoID);
  Run_Tree->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx);
  /*  Run_Tree->SetBranchAddress("eleE1x5", &eleE1x5);
  Run_Tree->SetBranchAddress("eleE2x5", &eleE2x5);
  Run_Tree->SetBranchAddress("eleE5x5", &eleE5x5);
  Run_Tree->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5);
  Run_Tree->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5);
  Run_Tree->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5);
  */Run_Tree->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5);
  /*Run_Tree->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed);
  Run_Tree->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt);
  Run_Tree->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt);
  Run_Tree->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt);
  Run_Tree->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt);
  Run_Tree->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt);
  Run_Tree->SetBranchAddress("elecaloEnergy", &elecaloEnergy);
  Run_Tree->SetBranchAddress("eleTrkdxy", &eleTrkdxy);
  Run_Tree->SetBranchAddress("eleKFHits", &eleKFHits);
  Run_Tree->SetBranchAddress("eleKFChi2", &eleKFChi2);
  */ Run_Tree->SetBranchAddress("eleIDbit", &eleIDbit);
    
  //########################################   Jet Info
  Run_Tree->SetBranchAddress("nJet",&nJet);
  Run_Tree->SetBranchAddress("jetPt",&jetPt);
  Run_Tree->SetBranchAddress("jetEta",&jetEta);
  Run_Tree->SetBranchAddress("jetPhi",&jetPhi);
  Run_Tree->SetBranchAddress("jetE",&jetE);
  Run_Tree->SetBranchAddress("jetCSV2BJetTags",&jetCSV2BJetTags);
  //  Run_Tree->SetBranchAddress("jetPFLooseId",&jetPFLooseId);
  Run_Tree->SetBranchAddress("jetID",&jetID);
  Run_Tree->SetBranchAddress("jetPUID",&jetPUID);
  Run_Tree->SetBranchAddress("jetRawPt",&jetRawPt);
  Run_Tree->SetBranchAddress("jetJECUnc",&jetJECUnc);
  Run_Tree->SetBranchAddress("jetRawE",&jetRawE);
  Run_Tree->SetBranchAddress("jetHadFlvr",&jetHadFlvr);
  Run_Tree->SetBranchAddress("jetP4Smear",&jetP4Smear);
  Run_Tree->SetBranchAddress("jetP4SmearUp",&jetP4SmearUp);
  Run_Tree->SetBranchAddress("jetP4SmearDo",&jetP4SmearDo);
    
    
    
  //########################################   MET Info
  Run_Tree->SetBranchAddress("pfMET",&pfMET);
  Run_Tree->SetBranchAddress("pfMET_T1UESUp",&pfMET_T1UESUp);
  Run_Tree->SetBranchAddress("pfMET_T1UESDo",&pfMET_T1UESDo);
  Run_Tree->SetBranchAddress("pfMET_T1JESUp",&pfMET_T1JESUp);
  Run_Tree->SetBranchAddress("pfMET_T1JESDo",&pfMET_T1JESDo);
    
  Run_Tree->SetBranchAddress("pfMETPhi",&pfMETPhi);
  Run_Tree->SetBranchAddress("pfMETPhi_T1UESUp",&pfMETPhi_T1UESUp);
  Run_Tree->SetBranchAddress("pfMETPhi_T1UESDo",&pfMETPhi_T1UESDo);
  Run_Tree->SetBranchAddress("pfMETPhi_T1JESUp",&pfMETPhi_T1JESUp);
  Run_Tree->SetBranchAddress("pfMETPhi_T1JESDo",&pfMETPhi_T1JESDo);
    
  Run_Tree->SetBranchAddress("metFilters",&metFilters);
  Run_Tree->SetBranchAddress("genHT",&genHT);
    
  Run_Tree->SetBranchAddress("pdfSystWeight",&pdfSystWeight);
  //  Run_Tree->SetBranchAddress("pdfSystWeight",&pdfSystWeightId);
  Run_Tree->SetBranchAddress("pdfWeight",&pdfWeight);
    
    
  return Run_Tree;
}


//########################################
// Pileup files
//########################################



TH1F *  HistPUData(){
  //std::cout<<"endPUDAta-----------------------------------------------------------------------"<<std::endl;
  TFile * PUData= TFile::Open("rootfile/pileup_hists/Data_nPU_new.root");
  TH1F * HistoPUData= (TH1F *) PUData->Get("pileup");
  HistoPUData->Rebin(2);
  HistoPUData->Scale(1.0/HistoPUData->Integral());
  cout << "HistoPUData integral= "<<HistoPUData->Integral()<<"\n";
  return HistoPUData;
}


TH1F *  HistPUMC(bool isData,TFile *f_Double){
  if (isData) return 0;
  else{
    //    TFile * PUMC= TFile::Open("../interface/pileup-hists/mcMoriondPU.root");
    //    TH1F * HistoPUMC= (TH1F *) PUMC->Get("pileup");
    //cout << "Coming in histpumc 1"<<"\n";
    TFile * PUMC= TFile::Open(f_Double->GetName());
    std::cout<<"PUMC= inputfile name"<<f_Double->GetName()<<std::endl;
    TH1F * HistoPUMC= (TH1F *) PUMC->Get("puTrue");
    //std::cout << "PUMC histogram"<<HistoPUMC->Integral()<<std::endl;
    //cout << "Coming in histpumc 2 "<<"\n";
    //    std::cout << "HistoPUMC"<<HistoPUMC->Integral()<<std::endl;
    HistoPUMC->Scale(1.0/HistoPUMC->Integral());
    //    std::cout << "HistoPUMC"<<HistoPUMC->Integral()<<std::endl;
    //cout << "Coming in histpumc 3 "<<"\n";
    return HistoPUMC;
  }
  return 0;
}



//########################################
// Muon Id, Iso, Trigger and Tracker Eff files
//########################################

TH2F**  FuncHistMuId(){
  //  std::cout<<"MU ID -----------------------------------------------------------------------"<<std::endl;                                                                                                   
  TFile * MuCorrId_BCDEF= TFile::Open(("rootfile/pileup_hists/RunABCD_SF_ID.root"));
  //  TH2F * HistoMuId_BCDEF= (TH2F *) MuCorrId_BCDEF->Get("NUM_TightID_DEN_genTracks_pt_abseta");
  TH2F * HistoMuId_BCDEF= (TH2F *) MuCorrId_BCDEF->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");
  static TH2F* HistoMuId[1]={HistoMuId_BCDEF};
  return  HistoMuId;
}

TH2F** FuncHistMuIso(){

  //  std::cout<<"MU Iso -----------------------------------------------------------------------"<<std::endl;                                                                                                   
  TFile * MuCorrIso_BCDEF= TFile::Open(("rootfile/pileup_hists/RunABCD_SF_ISO.root"));
  //TH2F * HistoMuIso_BCDEF= (TH2F *) MuCorrIso_BCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");
  //  TH2F * HistoMuIso_BCDEF= (TH2F *) MuCorrIso_BCDEF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_alleta_bin1");
  TH2F * HistoMuIso_BCDEF= (TH2F *) MuCorrIso_BCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");
  static  TH2F* HistoMuIso[1]={HistoMuIso_BCDEF};
  return HistoMuIso;}


TH2F** FuncHistMuTrigger(){
  //  std::cout<<"Trigger -----------------------------------------------------------------------"<<std::endl;                                                                                                  
  TFile * MuCorrTrg_BCDEF= TFile::Open(("rootfile/pileup_hists/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root"));
  //TH1F * HistoMuTrg_BCDEF= (TH1F *) MuCorrTrg_BCDEF->Get("Mu50_EtaBins/eta_ratio");
  TH2F * HistoMuTrg_BCDEF= (TH2F *) MuCorrTrg_BCDEF->Get("Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/pt_abseta_ratio");
  static TH2F* HistoMuTrg[2]={HistoMuTrg_BCDEF};
  return HistoMuTrg;
}


TGraphAsymmErrors* FuncHistMuTrack(){
  std::cout<<"track -----------------------------------------------------------------------"<<std::endl;                                                                                                   
  TFile * MuCorrTrack= TFile::Open(("rootfile/pileup_hists/Tracking_EfficienciesAndSF_BCDEFGH.root"));
  TGraphAsymmErrors * HistoMuTrack= (TGraphAsymmErrors *) MuCorrTrack->Get("ratio_eff_eta3_dr030e030_corr");
  return HistoMuTrack;
}


TH2F * FuncHistEleMVAId(std::string type){
  //      TFile * EleCorrMVAIdIso90= TFile::Open(("rootfile/pileup_hists/egammaEffi.txt_EGM2D.root")); //This is for 2016
  TFile * EleCorrMVAIdIso90= TFile::Open(("rootfile/pileup_hists/2018_ElectronMVA90.root")); //This is for 2018
  // TFile * EleCorrMVAIdIso90= TFile::Open(("../interface/pileup-hists/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp90iso.root"));
  TH2F * HistoEleMVAIdIso90= (TH2F *) EleCorrMVAIdIso90->Get("EGamma_SF2D");
  TH2F * HistoEleMVAIdIso90_EffMC= (TH2F *) EleCorrMVAIdIso90->Get("EGamma_EffMC2D");
  TH2F * HistoEleMVAIdIso90_EffData= (TH2F *) EleCorrMVAIdIso90->Get("EGamma_EffData2D");
  
  if (type.find("Tot") != string::npos)
    return HistoEleMVAIdIso90;
  else if (type.find("MC") != string::npos)
    return HistoEleMVAIdIso90_EffMC;
  else if (type.find("Data") != string::npos)
    return HistoEleMVAIdIso90_EffData;
  else
    return 0;
    
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float MuMass= 0.10565837;
float eleMass= 0.000511;


int getNumTau(){
    
    
  int numTau=0;
  for  (int itau=0 ; itau < nTau; itau++){
        
    if (tau_Pt->at(itau) < 20  || fabs(tau_Eta->at(itau)) > 2.3 ) continue;
        
    //bool TauIdIso =  taupfTausDiscriminationByDecayModeFinding->at(itau) > 0.5 && tauByLooseMuonRejection3->at(itau) > 0 && tauByMVA6LooseElectronRejection->at(itau) > 0 && tauByLooseIsolationMVArun2v1DBoldDMwLT->at(itau) > 0;
    bool TauIdIso =  0;//taupfTausDiscriminationByDecayModeFinding->at(itau) > 0.5 && tauByLooseMuonRejection3->at(itau) > 0 && tauByMVA6LooseElectronRejection->at(itau) > 0 && tauByLooseIsolationMVArun2v1DBoldDMwLT->at(itau) > 0;
    if (!TauIdIso) continue;
    numTau++;
  }
  return numTau;
}

int getNumElectron(){
    
    
  //            https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_recipes_for_2016
  int numElectron=0;
  float ElectronCor=1;
  float ElectronEffVeto=1;
  for  (int jele=0 ; jele < nEle; jele++){
        
    if ( elePt->at(jele) < 15 || fabs(eleEta->at(jele)) > 2.5) continue;
        
    bool eleMVAIdExtra= false;
    if (fabs (eleSCEta->at(jele)) <= 0.8 && eleMVAIsoID->at(jele) >   -0.83  ) eleMVAIdExtra= true;
    else if (fabs (eleSCEta->at(jele)) >  0.8 &&fabs (eleSCEta->at(jele)) <=  1.5 && eleMVAIsoID->at(jele) >   -0.77  ) eleMVAIdExtra= true;
    else if ( fabs (eleSCEta->at(jele)) >=  1.5 && eleMVAIsoID->at(jele) >  -0.69  ) eleMVAIdExtra= true;
    else eleMVAIdExtra= false;
        
        
    if (eleMVAIdExtra)
      numElectron++;
  }
  return numElectron;
    
}

int getNumZBoson(){

  int numZboson=0;
  for (int xmu=0; xmu< nMu; xmu++){
    for (int ymu=xmu+1; ymu< nMu; ymu++){
            
      TLorentzVector Mu4Momentum_0,Mu4Momentum_1,Z4Momentum;
      Mu4Momentum_0.SetPtEtaPhiM(muPt->at(xmu),muEta->at(xmu),muPhi->at(xmu),MuMass);
      Mu4Momentum_1.SetPtEtaPhiM(muPt->at(ymu),muEta->at(ymu),muPhi->at(ymu),MuMass);
      Z4Momentum=Mu4Momentum_1+Mu4Momentum_0;
            
      float IsoMu1=muPFChIso->at(xmu)/muPt->at(xmu);
      if ( (muPFNeuIso->at(xmu) + muPFPhoIso->at(xmu) - 0.5* muPFPUIso->at(xmu) )  > 0.0)
	IsoMu1= ( muPFChIso->at(xmu)/muPt->at(xmu) + muPFNeuIso->at(xmu) + muPFPhoIso->at(xmu) - 0.5* muPFPUIso->at(xmu))/muPt->at(xmu);
            
      float IsoMu2=muPFChIso->at(ymu)/muPt->at(ymu);
      if ( (muPFNeuIso->at(ymu) + muPFPhoIso->at(ymu) - 0.5* muPFPUIso->at(ymu) )  > 0.0)
	IsoMu2= ( muPFChIso->at(ymu)/muPt->at(ymu) + muPFNeuIso->at(ymu) + muPFPhoIso->at(ymu) - 0.5* muPFPUIso->at(ymu))/muPt->at(ymu);
            
      if ( muPt->at(xmu) > 60 && muPt->at(ymu) > 15 &&  (muIDbit->at(xmu) >> 1 & 1) & (muIDbit->at(ymu) >> 1 & 1) & IsoMu1 < 0.25  && IsoMu2 < 0.25 && Z4Momentum.M() > 80 && Z4Momentum.M()< 100  && (muCharge->at(xmu) * muCharge->at(ymu) < 0))
	numZboson++;
    }
  }
    
  return numZboson;
}



//###########       bJet multiplicity   ###########################################################
int numBJets( float BJetPtCut, float CSVCut){

  int numBJet=0;
  for (int ijet= 0 ; ijet < nJet ; ijet++){
    // std::cout<<"comin inside numbjet for loop"<<std::endl;
    if((jetID->at(ijet)>>0&1) == 1 && jetPt->at(ijet) > BJetPtCut && fabs(jetEta->at(ijet)) < 2.4 && jetCSV2BJetTags->at(ijet) >  CSVCut)
      numBJet++;
    //    std::cout<<"comin inside numbjet for loop"<<std::endl;
}
  return numBJet;
}


//###########       Jet multiplicity   ###########################################################
int numJets( float SimpleJetPtCut){
  int numJet=0;
  for (int ijet= 0 ; ijet < nJet ; ijet++){
    if ((jetID->at(ijet)>>0&1) == 1 && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4)
      numJet++;
  }
  return numJet;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



float compTopPtWeight(float topPt) {
  //Updated values for 13 TeV
  const float a =  0.0615  ;
  const float b =  -0.0005 ;
  //    const float a =  0.148 ;
  //    const float b =  -0.00129;
  //    const float a = 0.156;
  //    const float b = -0.00137;
  return TMath::Exp(a + b * topPt);
}

float compTopPtWeight(float top1Pt, float top2Pt) {
  //std::cout << "<compTopPtWeight>:" << std::endl;
  float topPtWeight2 = compTopPtWeight(top1Pt) * compTopPtWeight(top2Pt);
  //std::cout << " top1Pt = " << top1Pt << ", top2Pt = " << top2Pt << ": topPtWeight2 = " << topPtWeight2 << std::endl;
  return ( topPtWeight2 > 0.) ? TMath::Sqrt(topPtWeight2) : 0.;
}


std::vector<float> GeneratorInfo(){                                                                                                                                                         
 

  vector<float>    infoGen;                                                                                                                                                                               
  infoGen.clear();                                                                                                                                                                                        
  float GenTopPt=0;                                                                                                                                                                                       
  float GenAntiTopPt=0;                                                                                                                                                                                   
  float TopPtReweighting = 1;                                                                                                                                                                             
  float WBosonPt=0;                                                                                                                                                                                       
  float WBosonKFactor=1;                                                                                                                                                                                  
  float ZBosonPt=0;                                                                                                                                                                                     
  float ZBosonKFactor=1;                                                                                                                                                                                 
  int modPDGId=-10;                                                                                                                                                                                     
  int AntimodPDGId=-10;                                                                                                                                                                                   
  float WBosonMass=0;                                                                                                                                                                                     
  TLorentzVector GenMu4Momentum,GenAntiMu4Momentum, WGEN4Momentum, MUGEN4Momentum, NUGEN4Momentum;                                                                                                      
  for (int igen=0;igen < nMC; igen++){                                                                                                                                                                   
    if (mcPID->at(igen) == 6 && mcStatus->at(igen) ==62) GenTopPt=mcPt->at(igen) ;                                                                                                                        
    if (mcPID->at(igen) == -6 && mcStatus->at(igen) ==62) GenAntiTopPt=mcPt->at(igen);                                                                                                                    
    //W Pt
    if (fabs(mcPID->at(igen)) ==24   && mcStatus->at(igen) ==22)  {WBosonPt= mcPt->at(igen); WBosonMass=mcMass->at(igen);}
    if ( fabs(mcPID->at(igen)) ==13 && mcStatus->at(igen) ==1 )  {MUGEN4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen));}
    if ( fabs(mcPID->at(igen)) ==14  && mcStatus->at(igen) ==1)  {NUGEN4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen));}
        
    //W Pt
    if (fabs(mcPID->at(igen)) ==24   && mcStatus->at(igen) ==22)  {WBosonPt= mcPt->at(igen); WBosonMass=mcMass->at(igen);}
    if ( fabs(mcPID->at(igen)) ==13 && mcStatus->at(igen) ==1 )  {MUGEN4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen));}
    if ( fabs(mcPID->at(igen)) ==14  && mcStatus->at(igen) ==1)  {NUGEN4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen));}
   
    //Z Pt
    if (fabs(mcPID->at(igen)) ==23)  ZBosonPt= mcPt->at(igen); //FIXME somethime we do not have Z in the DY events
    //    if ( mcPID->at(igen) ==13  )  {GenMu4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen)); modPDGId=mcPID->at(igen);}
    //if ( mcPID->at(igen) ==-13  )  {GenAntiMu4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen)); AntimodPDGId=mcPID->at(igen);}

    if ( mcPID->at(igen) ==13  )  {GenMu4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen));}
    if ( mcPID->at(igen) ==-13  )  {GenAntiMu4Momentum.SetPtEtaPhiM(mcPt->at(igen),mcEta->at(igen),mcPhi->at(igen),mcMass->at(igen));}

  }
  TopPtReweighting=compTopPtWeight(GenTopPt, GenAntiTopPt);
  if (ZBosonPt ==0)
    ZBosonPt=(GenMu4Momentum+GenAntiMu4Momentum).Pt();  //This is a temp solution to the above problem
 
  if (WBosonPt==0)
    WBosonPt = (MUGEN4Momentum+NUGEN4Momentum).Pt();
            
  //######################## Top Pt Reweighting
  infoGen.push_back(TopPtReweighting);
    
  //######################## W Pt
  infoGen.push_back(WBosonPt);
    
  //######################## Z Pt
  infoGen.push_back(ZBosonPt);
    
  //######################## W Mass
  infoGen.push_back(WBosonMass);
    
  return infoGen;


}                                            

//########################################
// Btagging scale factor This part is the outcome of the CodexAnalyzer_BTagEff.cc
//########################################

TH2F **  FuncHistBTagSF(){
  TFile * TTEff= TFile::Open(("rootfile/BTagSF/TTJets.root"));
  TH2F * TTSF0_btagged= (TH2F *) TTEff->Get("BSF_FLV0_Btagged");
  TH2F * TTSF0_total= (TH2F *) TTEff->Get("BSF_FLV0_Total");
  TH2F * TTSF5_btagged= (TH2F *) TTEff->Get("BSF_FLV5_Btagged");
  TH2F * TTSF5_total= (TH2F*) TTEff->Get("BSF_FLV5_Total");
    
  static TH2F * Btagg_TT[4]={TTSF0_btagged,TTSF0_total,TTSF5_btagged,TTSF5_total};
    
  return Btagg_TT;
}

//        TFile * DataEff= TFile::Open(("OutFiles_BTagSF/Data.root"));
//        TH2F * DataSF0_btagged= (TH2F *) DataEff->Get("BSF_FLV0_Btagged");
//        TH2F * DataSF0_total= (TH2F *) DataEff->Get("BSF_FLV0_Total");
//        TH2F * DataSF5_btagged= (TH2F *) DataEff->Get("BSF_FLV5_Btagged");
//        TH2F * DataSF5_total= (TH2F *) DataEff->Get("BSF_FLV5_Total");




//########################################
// W and DY K-factor files  (FIT-based K-factor)
//########################################

TFile * kfactorW=TFile::Open("rootfile/kfactor_W.root");
TH1F* HistkfactorW= (TH1F*) kfactorW->Get("KFcator");
float kf_W_1=HistkfactorW->GetBinContent(1);
float kf_W_2=HistkfactorW->GetBinContent(2);


TFile * kfactorWUp=TFile::Open("rootfile/kfactor_monoJet_WUp.root");
TH1F* HistkfactorWUp= (TH1F*) kfactorWUp->Get("KFcator");
float kf_W_1Up=HistkfactorWUp->GetBinContent(1);
float kf_W_2Up=HistkfactorWUp->GetBinContent(2);

TFile * kfactorWDown=TFile::Open("rootfile/kfactor_monoJet_WDown.root");
TH1F* HistkfactorWDown= (TH1F*) kfactorWDown->Get("KFcator");
float kf_W_1Down=HistkfactorWDown->GetBinContent(1);
float kf_W_2Down=HistkfactorWDown->GetBinContent(2);



TFile * kfactorZ=TFile::Open("rootfile/kfactor_Z.root");
TH1F* HistkfactorZ= (TH1F*) kfactorZ->Get("KFcator");
float kf_Z_1=HistkfactorZ->GetBinContent(1);
float kf_Z_2=HistkfactorZ->GetBinContent(2);


TFile * kfactorZUp=TFile::Open("rootfile/kfactor_monoJet_ZUp.root");
TH1F* HistkfactorZUp= (TH1F*) kfactorZUp->Get("KFcator");
float kf_Z_1Up=HistkfactorZUp->GetBinContent(1);
float kf_Z_2Up=HistkfactorZUp->GetBinContent(2);

TFile * kfactorZDown=TFile::Open("rootfile/kfactor_monoJet_ZDown.root");
TH1F* HistkfactorZDown= (TH1F*) kfactorZDown->Get("KFcator");
float kf_Z_1Down=HistkfactorZDown->GetBinContent(1);
float kf_Z_2Down=HistkfactorZDown->GetBinContent(2);


float FuncBosonKFactor(std::string X){
    
  if (X.find("W1Cen") != string::npos)
    return kf_W_1;
  else if (X.find("W2Cen") != string::npos)
    return kf_W_2;
  else if (X.find("W1Up") != string::npos)
    return kf_W_1Up;
  else if (X.find("W2Up") != string::npos)
    return kf_W_2Up;
  else if (X.find("W1Down") != string::npos)
    return kf_W_1Down;
  else if (X.find("W2Down") != string::npos)
    return kf_W_2Down;
    
  else if (X.find("Z1Cen") != string::npos)
    return kf_Z_1;
  else if (X.find("Z2Cen") != string::npos)
    return kf_Z_2;
  else if (X.find("Z1Up") != string::npos)
    return kf_Z_1Up;
  else if (X.find("Z2Up") != string::npos)
    return kf_Z_2Up;
  else if (X.find("Z1Down") != string::npos)
    return kf_Z_1Down;
  else if (X.find("Z2Down") != string::npos)
    return kf_Z_2Down;
    
    
  else
    return 0;
}





//###########       getBtagEfficiency using TTbar samples from CodexAnalyzer_BTagEff.cc   ###########################################################
float getBtagEfficiency(bool isData, bool passCSV, float pt, float eta, TH2F ** Btagg_TT){
    
    
  if ( isData) return 1;
    
    
    
  int ptBIN;
  if ( pt < 50 ) ptBIN=1;
  if (pt >= 50 && pt < 70 ) ptBIN=2;
  if (pt >= 70 && pt < 100 ) ptBIN=3;
  if (pt >= 100 && pt < 140) ptBIN=4;
  if (pt >= 140 && pt < 200) ptBIN=5;
  if (pt >= 200 && pt < 300) ptBIN=6;
  if (pt >= 300 && pt < 600) ptBIN=7;
  if (pt >= 600 ) ptBIN=8;
    
  int etaBIN;
  if (eta >= 0 && eta < 0.8 ) etaBIN=1;
  if (eta >= 0.8 && eta < 1.5 ) etaBIN=2;
  if (eta >= 1.5 ) etaBIN=3;
    
    
    
  TH2F * TTSF0_btagged=Btagg_TT[0];
  TH2F * TTSF0_total=Btagg_TT[1];
  TH2F * TTSF5_btagged=Btagg_TT[2];
  TH2F * TTSF5_total=Btagg_TT[3];
    
  //    cout << "Btag efficiency is = "<< pt << " ptBIN " <<ptBIN << "   "<<eta << " etaBIN " << etaBIN <<"  ratio=  " <<TTSF0_btagged->GetBinContent(ptBIN,etaBIN) << "    "<<TTSF0_total->GetBinContent(ptBIN,etaBIN) <<"\n";
    
    
  if (passCSV)
    return  TTSF5_btagged->GetBinContent(ptBIN,etaBIN)*1.0/TTSF5_total->GetBinContent(ptBIN,etaBIN);
  else
    return  TTSF0_btagged->GetBinContent(ptBIN,etaBIN)*1.0/TTSF0_total->GetBinContent(ptBIN,etaBIN);
    
    
}


float FuncFinalBTagSF(bool isData, TH2F ** Btagg_TT, float BJetPtCut, float CSVCut){
 
  float EffJet =1;
  float SF=1;
  float P_Data_P_mc=1;
  float FinalBTagSF=1;
    
  if (isData) return 1;
  for (int ijet= 0 ; ijet < nJet ; ijet++){
        
    float HadronFlavor= isData ? 1 : jetHadFlvr->at(ijet);
        
    if (/*jetPFLooseId->at(ijet) > 0.5 &&*/ jetPt->at(ijet) > BJetPtCut && fabs(jetEta->at(ijet)) < 2.4 ){
            
            
      if ( jetCSV2BJetTags->at(ijet) >  CSVCut ){
	EffJet= getBtagEfficiency( isData, 1,  jetPt->at(ijet), fabs(jetEta->at(ijet)), Btagg_TT);
	SF= GetBJetSF(isData, jetPt->at(ijet), jetPt->at(ijet), HadronFlavor);
	P_Data_P_mc=SF*EffJet/EffJet;
                
                
      }
      else{
	EffJet= getBtagEfficiency( isData, 0,  jetPt->at(ijet), fabs(jetEta->at(ijet)), Btagg_TT);
	SF=GetBJetSF(isData,jetPt->at(ijet), jetPt->at(ijet), HadronFlavor);
	P_Data_P_mc=(1-SF*EffJet)/(1-EffJet);
                
      }
      FinalBTagSF *=P_Data_P_mc;
    }
        
    //        FinalBTagSF *=P_Data_P_mc; //  Seemd to ve a BUGGGGGGG  May16
  }
  return FinalBTagSF;


}



float getElectronCor(TH2F * HistoEleMVAIdIso90){
    
  float ElectronCor=1;
  for  (int jele=0 ; jele < nEle; jele++){
        
    if ( elePt->at(jele) < 15 || fabs(eleEta->at(jele)) > 2.5) continue;
        
    bool eleMVAIdExtra= false;
    if (fabs (eleSCEta->at(jele)) <= 0.8 && eleMVAIsoID->at(jele) >   -0.83  ) eleMVAIdExtra= true;
    else if (fabs (eleSCEta->at(jele)) >  0.8 &&fabs (eleSCEta->at(jele)) <=  1.5 && eleMVAIsoID->at(jele) >   -0.77  ) eleMVAIdExtra= true;
    else if ( fabs (eleSCEta->at(jele)) >=  1.5 && eleMVAIsoID->at(jele) >  -0.69  ) eleMVAIdExtra= true;
    else eleMVAIdExtra= false;
        
        
        
    //        if (!(eleMVAIdExtra )) {
    //            ElectronEffVeto= ElectronEffVeto * getEffVetoMVA90WPElectron80X(isData,  elePt->at(jele),eleSCEta->at(jele),    HistoEleMVAIdIso90 , HistoEleMVAIdIso90_EffMC,HistoEleMVAIdIso90_EffData);
    //            continue;
    //        }
        
    if (eleMVAIdExtra)
      ElectronCor=getCorrFactorMVA90WPElectron94X(isData,  elePt->at(jele),eleSCEta->at(jele),  HistoEleMVAIdIso90 );
        
    break;
  }
  return ElectronCor;
}
