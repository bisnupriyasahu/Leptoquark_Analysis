#ifndef MuonClass_h
#define MuonClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
// Header file for the classes stored in the TTree if any.                                                                                    
#include "TLorentzVector.h" 
#include <iostream>
#include <fstream>
#include <sstream>
#include "vector"
#include <vector>
#include <TH2.h>
#include <TH1.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TChain.h>
#include "math.h"
#include "TNtuple.h"
#include <stdio.h>
#include <utility>

using namespace std;

class MuonClass {
 public :
  //TTree          *fChain;                                                                        
  TTree *tree1; 
  Int_t           fCurrent;                                                                                 
  TTree          *fChain;



  TFile *fileName;
  TTree *tree;

  TH1F *MUJ_recoHT;
  TH1F *MUJ_ST;
  TH1F *Muon_nMu;
  TH1F *Muon_MuPt;
  TH1F *Muon_MuEta;
  TH1F *Muon_MuPhi;
  TH1F *Jet_JetEta;
  TH1F *Jet_JetPhi;
  int getNumJets(float SimpleJetPtCut); 
  int getNumBJets(float BJetPtCut,float CSVCut);
  int getNumZBoson();
  float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi);
  //  auto numOf_c_quark=0;
  //auto numOf_s_quark=0;



  // Fixed size dimensions of array or collections stored in the TTree if any.                                                                 

  // Declaration of leaf types  



  Int_t           nMu;
  vector<float>   *muPt;
  vector<float>   *muEn;
  vector<float>   *muEta;
  vector<float>   *muPhi;
  vector<int>     *muCharge;
  vector<int>     *muType;
  vector<unsigned short> *muIDbit;
  vector<float>   *muD0;
  vector<float>   *muDz;
  vector<float>   *muSIP;
  vector<float>   *muChi2NDF;
  vector<float>   *muInnerD0;
  vector<float>   *muInnerDz;
  vector<int>     *muTrkLayers;
  vector<int>     *muPixelLayers;
  vector<int>     *muPixelHits;
  vector<int>     *muMuonHits;
  vector<int>     *muStations;
  vector<int>     *muMatches;
  vector<int>     *muTrkQuality;
  vector<float>   *muIsoTrk;
  vector<float>   *muPFChIso;
  vector<float>   *muPFPhoIso;
  vector<float>   *muPFNeuIso;
  vector<float>   *muPFPUIso;
  vector<float>   *muPFMiniIso;
  vector<unsigned int> *muFiredTrgs;
  vector<unsigned int> *muFiredL1Trgs;
  vector<float>   *muInnervalidFraction;
  vector<float>   *musegmentCompatibility;
  vector<float>   *muchi2LocalPosition;
  vector<float>   *mutrkKink;
  vector<float>   *muBestTrkPtError;
  vector<float> *muBestTrkPt;
  Bool_t isData;
  ULong64_t HLTEleMuX;

  //////////////////////MET////////////////////////////////////////////

  Int_t           metFilters;
  Float_t         pfMET;
  Float_t         pfMETPhi;
  Float_t         pfMET_T1JERUp;
  Float_t         pfMET_T1JERDo;
  Float_t         pfMET_T1JESUp;
  Float_t         pfMET_T1JESDo;
  Float_t         pfMET_T1UESUp;
  Float_t         pfMET_T1UESDo;
  Float_t         pfMETPhi_T1JESUp;
  Float_t         pfMETPhi_T1JESDo;
  Float_t         pfMETPhi_T1UESUp;
  Float_t         pfMETPhi_T1UESDo;





  //////////////////jets//////////////////////////////////////////



  Int_t           nJet;
  vector<float>   *jetPt;
  vector<float>   *jetE;
  vector<float>   *jetEta;
  vector<float>   *jetPhi;
  vector<float>   *jetRawPt;
  vector<float>   *jetRawEn;
  vector<float>   *jetMt;
  vector<float>   *jetArea;
  vector<float>   *jetLeadTrackPt;
  vector<float>   *jetLeadTrackEta;
  vector<float>   *jetLeadTrackPhi;
  vector<int>     *jetLepTrackPID;
  vector<float>   *jetLepTrackPt;
  vector<float>   *jetLepTrackEta;
  vector<float>   *jetLepTrackPhi;
  vector<float>   *jetCSV2BJetTags;
  vector<float>   *jetDeepCSVTags_b;
  vector<float>   *jetDeepCSVTags_bb;
  vector<float>   *jetDeepCSVTags_c;
  vector<float>   *jetDeepCSVTags_udsg;
  vector<bool>    *jetPFLooseId;
  vector<int>     *jetID;
  vector<float>   *jetPUID;
  vector<int>     *jetPUFullID;
  vector<float>   *jetJECUnc;
  vector<ULong64_t> *jetFiredTrgs;
  vector<float>   *jetCHF;
  vector<float>   *jetNHF;
  vector<float>   *jetCEF;
  vector<float>   *jetNEF;
  vector<int>     *jetNCH;
  vector<int>     *jetNNP;
  vector<float>   *jetMUF;
  vector<float>   *jetVtxPt;
  vector<float>   *jetVtxMass;
  vector<float>   *jetVtxNtrks;
  vector<float>   *jetVtx3DVal;
  vector<float>   *jetVtx3DSig;



  float LeptonPtCut_=60;
  float LeptonIsoCut=0.15;
  float JetPtCut=100;
  float MuMass= 0.10565837;
  float SimpleJetPtCut=30;
  float BJetPtCut=30;
  float CSVCut=   0.9693   ; 
  //  float LeptonPtCut_=60;


  // List of branches    
  TBranch        *b_nMu;   //! 
  TBranch        *b_muPt;
  TBranch        *b_muEn;
  TBranch        *b_muEta;
  TBranch        *b_muPhi;
  TBranch        *b_muCharge;
  TBranch        *b_muType;
  TBranch        *b_muIDbit;
  TBranch        *b_muD0;
  TBranch        *b_muDz;
  TBranch        *b_muSIP;
  TBranch        *b_muChi2NDF;
  TBranch        *b_muInnerD0;
  TBranch        *b_muInnerDz;
  TBranch        *b_muTrkLayers;
  TBranch        *b_muPixelLayers;
  TBranch        *b_muPixelHits;
  TBranch        *b_muMuonHits;
  TBranch        *b_muStations;
  TBranch        *b_muMatches;
  TBranch        *b_muTrkQuality;
  TBranch        *b_muIsoTrk;
  TBranch        *b_muPFChIso;
  TBranch        *b_muPFPhoIso;
  TBranch        *b_muPFNeuIso;
  TBranch        *b_muPFPUIso;
  TBranch        *b_muPFMiniIso;
  TBranch        *b_muFiredTrgs;
  TBranch        *b_muFiredL1Trgs;
  TBranch        *b_muInnervalidFraction;
  TBranch        *b_musegmentCompatibility;
  TBranch        *b_muchi2LocalPosition;
  TBranch        *b_mutrkKink;
  TBranch        *b_muBestTrkPtError;
  TBranch        *b_muBestTrkPt;
  TBranch        *b_HLTEleMuX;
  TBranch        *b_isData;
  TBranch        *b_metFilters; 
  TBranch        *b_pfMET;   //!                                                                                                     
  TBranch        *b_pfMETPhi;   //!                                                                                                  
  TBranch        *b_pfMET_T1JERUp;   //!                                                                                             
  TBranch        *b_pfMET_T1JERDo;   //!                                                                                             
  TBranch        *b_pfMET_T1JESUp;   //!                                                                                             
  TBranch        *b_pfMET_T1JESDo;   //!                                                                                             
  TBranch        *b_pfMET_T1UESUp;   //!                                                                                             
  TBranch        *b_pfMET_T1UESDo;   //!                                                                                             
  TBranch        *b_pfMETPhi_T1JESUp;   //!                                                                                          
  TBranch        *b_pfMETPhi_T1JESDo;   //!                                                                                          
  TBranch        *b_pfMETPhi_T1UESUp;   //!                                                                                          
  TBranch        *b_pfMETPhi_T1UESDo;   //! 
  TBranch        *b_nJet;   //!                                                                                                               
  TBranch        *b_jetPt;   //!                                                                                                              
  TBranch        *b_jetE;   //!                                                                                                              
  TBranch        *b_jetEta;   //!                                                                                                             
  TBranch        *b_jetPhi;   //!                                                                                                             
  TBranch        *b_jetRawPt;   //!                                                                                                           
  TBranch        *b_jetRawEn;   //!                                                                                                           
  TBranch        *b_jetMt;   //!                                                                                                              
  TBranch        *b_jetArea;   //!                                                                                                            
  TBranch        *b_jetLeadTrackPt;   //!                                                                                                     
  TBranch        *b_jetLeadTrackEta;   //!                                                                                                    
  TBranch        *b_jetLeadTrackPhi;   //!                                                                                                    
  TBranch        *b_jetLepTrackPID;   //!                                                                                                     
  TBranch        *b_jetLepTrackPt;   //!                                                                                                      
  TBranch        *b_jetLepTrackEta;   //!                                                                                                     
  TBranch        *b_jetLepTrackPhi;   //!                                                                                                     
  TBranch        *b_jetCSV2BJetTags;   //!                                                                                                    
  TBranch        *b_jetDeepCSVTags_b;   //!                                                                                                   
  TBranch        *b_jetDeepCSVTags_bb;   //!                                                                                                  
  TBranch        *b_jetDeepCSVTags_c;   //!                                                                                                   
  TBranch        *b_jetDeepCSVTags_udsg;   //!                                                                                                
  TBranch        *b_jetPFLooseId;   //!                                                                                                       
  TBranch        *b_jetID;   //!                                                                                                              
  TBranch        *b_jetPUID;   //!                                                                                                            
  TBranch        *b_jetPUFullID;   //!                                                                                                        
  TBranch        *b_jetJECUnc;   //!                                                                                                          
  TBranch        *b_jetFiredTrgs;   //!                                                                                                       
  TBranch        *b_jetCHF;   //!               
  TBranch        *b_jetNHF;   //!                                                                                                             
  TBranch        *b_jetCEF;   //!                                                                                                             
  TBranch        *b_jetNEF;   //!                                                                                                             
  TBranch        *b_jetNCH;   //!                                                                                                             
  TBranch        *b_jetNNP;   //!                                                                                                             
  TBranch        *b_jetMUF;   //!                                                                                                             
  TBranch        *b_jetVtxPt;   //!                                                                                                           
  TBranch        *b_jetVtxMass;   //!                                                                                                         
  TBranch        *b_jetVtxNtrks;   //!                                                                                                        
  TBranch        *b_jetVtx3DVal;   //!                                                                                                        
  TBranch        *b_jetVtx3DSig;   //!  



  MuonClass(const char* file1, const char* file2);
  virtual ~MuonClass();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TChain *tree);
  virtual void     Loop(Long64_t maxEvents, int reportEvery);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void Histograms(const char* file2);
  /////////////////////////////////////number functions declaration///////////////////////////////////////
 
  



  /////////////////////////////////////number functions declaration///////////////////////////////////////                                                                                               
 




};

#endif

 
#ifdef MuonClass_cxx
MuonClass::MuonClass(const char* file1, const char* file2)
{

  TChain *chain = new TChain("phoJetNtuplizer/eventTree");
  TString path = file1;


  TSystemDirectory sourceDir("hi",path);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter next(fileList);
  TSystemFile* filename;
  int fileNumber = 0;
  int maxFiles = -1;
  std::cout<<"path:"<<path<<std::endl;
  //std::cout<<"COMMING IN BEFORE PATH"<<std::endl;
  while ((filename = (TSystemFile*)next()) && fileNumber >  maxFiles)
    {
      // std::cout<<"comin 1"<<std::endl;
	if(fileNumber > 1)                                                                                                                    
	  {                                                                                                                              
	    TString dataset = "";
	    TString  FullPathInputFile = (path+filename->GetName());
	    TString name = filename->GetName();
	   //   std::cout<<"comin 2"<<std::endl;
	    if(name.Contains(dataset))
	      {
		std::cout<<"FullPathInputFile:"<<FullPathInputFile<<std::endl;
		chain->Add(FullPathInputFile);
	      }//name dataset
	  }//file no
	    fileNumber++;
    }//whileloop
	//  std::cout<<"comin 3"<<std::endl;
	
  
  std::cout<<"All files added."<<std::endl;
  std::cout<<"Initializing chain."<<std::endl;
  Init(chain);
  Histograms(file2);
}//Muonclass
      



MuonClass::~MuonClass()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  fileName->cd();
  fileName->Write();
  fileName->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////         Number functonf for Jets,BJets,ZBoson            //////////////////////////////////////////////////////////////////////



//###########       bJet multiplicity   ###########################################################

int MuonClass::getNumBJets(float BJetPtCut,float CSVCut)
{
  int numBJet=0;
  for (int ijet= 0 ; ijet < nJet ; ijet++){
    //   if (jetPFLooseId->at(ijet) > 0.5 && jetPt->at(ijet) > BJetPtCut && fabs(jetEta->at(ijet)) < 2.4  && jetCSV2BJetTags->at(ijet) >  CSVCut)
    if((*jetID)[ijet]>>0&1 == 1 && jetPt->at(ijet) > BJetPtCut && fabs(jetEta->at(ijet)) < 2.4 && jetCSV2BJetTags->at(ijet) >  CSVCut)
      numBJet++;
    //(jetID->at(ijet)>>1&1) ==1;
    //(jetID->at(ijet)>>0&&1 ==1) && (jetID->at(ijet)>>0&&1 > 0.5
    std::cout<<"numBJet  ijet "<<ijet<<std::endl;
    std::cout<<"numBJet   "<<numBJet<<std::endl;
    //    std::cout<<"jetID->at(ijet) >>0&&1 ==1  "<<((jetID->at(ijet)>>0&1) == 1)<<std::endl;
  }
  return numBJet;
}

//###########       Jet multiplicity   ###########################################################
int MuonClass::getNumJets(float SimpleJetPtCut){
  int numJet=0;
  
  //std::cout<<"jetPt->at(ijet) "<<jetPt->at(ijet)<<std::endl;
  std::cout<<"SimpleJetPtCut      "<<SimpleJetPtCut<<std::endl;
  std::cout<<"nJet   "<<nJet<<std::endl;
  for (int ijet=0;ijet<nJet;ijet++){
    std::cout<<"iJet   "<<ijet<<std::endl;   
    // if (jetPFLooseId->at(ijet) > 0.5 && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4 )
    // if (jetPFLooseId->at(ijet) > 0.5 )
    if ((*jetID)[ijet]>>0&1 == 1 && jetPt->at(ijet) > SimpleJetPtCut && fabs(jetEta->at(ijet)) < 2.4)
      numJet++;
    std::cout<<"numJet   "<<numJet<<std::endl;
    // std::cout<<"jetPFLooseId->at(ijet)   "<<jetPFLooseId->at(ijet)<<std::endl;
  }
  return numJet;
}

//###########       ZBoson multiplicity   ###########################################################

Int_t MuonClass::getNumZBoson(){    
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//##################################              TMass_F                ########################################################################
float MuonClass::TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
  return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t MuonClass::GetEntry(Long64_t entry)
{
    // Read contents of entry.                                                                                                                     
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
  }
  Long64_t MuonClass::LoadTree(Long64_t entry)
  {
    // Set the environment to read one entry                                                                                                       
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
    }
    return centry;
  }

  void MuonClass::Init(TChain *tree)
  {
    // Set object pointer
    muPt = 0;
    muEn = 0;
    muEta = 0;
    muPhi = 0;
    muCharge = 0;
    muType = 0;
    muIDbit = 0;
    muD0 = 0;
    muDz = 0;
    muSIP = 0;
    muChi2NDF = 0;
    muInnerD0 = 0;
    muInnerDz = 0;
    muTrkLayers = 0;
    muPixelLayers = 0;
    muPixelHits = 0;
    muMuonHits = 0;
    muStations = 0;
    muMatches = 0;
    muTrkQuality = 0;
    muIsoTrk = 0;
    muPFChIso = 0;
    muPFPhoIso = 0;
    muPFNeuIso = 0;
    muPFPUIso = 0;
    muFiredTrgs = 0;
    muFiredL1Trgs = 0;
    muInnervalidFraction = 0;
    musegmentCompatibility = 0;
    muchi2LocalPosition = 0;
    mutrkKink = 0;
    muBestTrkPtError = 0;
    muBestTrkPt = 0;
    jetPt = 0;
    jetE = 0;
    jetEta = 0;
    jetPhi = 0;
    jetRawPt = 0;
    jetRawEn = 0;
    jetMt = 0;
    jetArea = 0;
    jetLeadTrackPt = 0;
    jetLeadTrackEta = 0;
    jetLeadTrackPhi = 0;
    jetLepTrackPID = 0;
    jetLepTrackPt = 0;
    jetLepTrackEta = 0;
    jetLepTrackPhi = 0;
    jetCSV2BJetTags = 0;
    jetDeepCSVTags_b = 0;
    jetDeepCSVTags_bb = 0;
    jetDeepCSVTags_c = 0;
    jetDeepCSVTags_udsg = 0;
    jetPFLooseId = 0;
    jetID = 0;
    jetPUID = 0;
    jetPUFullID = 0;
    jetJECUnc = 0;
    jetFiredTrgs = 0;
    jetCHF = 0;
    jetNHF = 0;
    jetCEF = 0;
    jetNEF = 0;
    jetNCH = 0;
    jetNNP = 0;
    jetMUF = 0;
    jetVtxPt = 0;
    jetVtxMass = 0;
    jetVtxNtrks = 0;
    jetVtx3DVal = 0;
    jetVtx3DSig = 0;






    // Set branch addresses and branch pointers                                                                                                 
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);


    std::cout<<"comming"<<std::endl;

    fChain->SetBranchAddress("isData", &isData, &b_isData);
    fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
    fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
    fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
    fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
    fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
    fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
    fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
    fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
    fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
    fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
    fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
    fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
    fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
    fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);
    fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
    fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
    fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
    fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
    fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
    fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
    fChain->SetBranchAddress("muType", &muType, &b_muType);
    fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
    fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
    fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
    fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
    fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
    fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
    fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
    fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
    fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
    fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
    fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
    fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
    fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
    fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
    fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
    fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
    fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
    fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
    fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
    fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
    fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
    fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
    fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
    fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
    fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
    fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
    fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
    fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
    fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
    fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
    fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
    fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
    fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
    fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
    fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
    fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
    fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
    fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
    fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
    fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
    fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
    fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
    fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
    fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
    fChain->SetBranchAddress("jetDeepCSVTags_b", &jetDeepCSVTags_b, &b_jetDeepCSVTags_b);
    fChain->SetBranchAddress("jetDeepCSVTags_bb", &jetDeepCSVTags_bb, &b_jetDeepCSVTags_bb);
    fChain->SetBranchAddress("jetDeepCSVTags_c", &jetDeepCSVTags_c, &b_jetDeepCSVTags_c);
    fChain->SetBranchAddress("jetDeepCSVTags_udsg", &jetDeepCSVTags_udsg, &b_jetDeepCSVTags_udsg);
    fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
    fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
    fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
    fChain->SetBranchAddress("jetPUFullID", &jetPUFullID, &b_jetPUFullID);
    fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
    fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
    fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
    fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
    fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
    fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
    fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
    fChain->SetBranchAddress("jetNNP", &jetNNP, &b_jetNNP);
    fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
    fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
    fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
    fChain->SetBranchAddress("jetVtxNtrks", &jetVtxNtrks, &b_jetVtxNtrks);
    fChain->SetBranchAddress("jetVtx3DVal", &jetVtx3DVal, &b_jetVtx3DVal);
    fChain->SetBranchAddress("jetVtx3DSig", &jetVtx3DSig, &b_jetVtx3DSig);


    Notify();
  }




    Bool_t MuonClass::Notify()
    {
      // The Notify() function is called when a new file is opened. This                                                                          
      // can be either for a new TTree in a TChain or when when a new TTree                                                                       
      // is started when using PROOF. It is normally not necessary to make changes                                                                
      // to the generated code, but the routine can be extended by the                                                                            
      // user if needed. The return value is currently not used.                                                                                  

      return kTRUE;
    }

    void MuonClass::Show(Long64_t entry)
    {
      // Print contents of entry.                                                                                                                    
      // If entry is not specified, print current entry                                                                                              
      if (!fChain) return;
      fChain->Show(entry);
    }
    Int_t MuonClass::Cut(Long64_t entry)
    {
      // This function may be called from Loop.                                                                                              
      // returns  1 if entry is accepted.                                                                                                    
      // returns -1 otherwise.                                                                                                                       
      return 1;
    }
#endif // #ifdef MuonClass_cxx                                                                                                               


