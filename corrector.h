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
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
using namespace std;



float LandauFunc(float x){
    
  //return Landau(float x, Double_t mpv = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
    
  float Land= 9.794 * TMath::Landau(x, 210, 137);
  if (x > 200)  Land= 9.794 * TMath::Landau(200, 210, 137);
    
  return Land;
    
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////      getCorrFactorMuon94X       /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int ptBIN=0;
int etaBIN=0;
int etaPOINT=-1;


/////////////////////////////////////////////////////
//  Muon Id Correction 94X
////////////////////////////////////////////////////////////


float Cor94X_ID_Mu(float pt,float eta, TH2F* HistoId){
  if (pt > 100 ) pt=100;
  return HistoId->GetBinContent(HistoId->GetXaxis()->FindBin(pt),HistoId->GetYaxis()->FindBin(fabs(eta)));
}


/////////////////////////////////////////////////////
//  Muon Iso Correction 94X
////////////////////////////////////////////////////////////
float Cor94X_Iso_Mu(float pt, float eta , TH2F * HistoIso) {
  if (pt > 100 ) pt=100;
  return HistoIso->GetBinContent(HistoIso->GetXaxis()->FindBin(pt),HistoIso->GetYaxis()->FindBin(fabs(eta)));
}


///////////////////////////////////////////////
//  Muon Trigger Correction 94X  eta-binned only
////////////////////////////////////////////////////////////

float Cor94X_Trigger_Mu_onlyEta(float eta, TH1F* HistoTrg ){
  return HistoTrg->GetBinContent(HistoTrg->GetXaxis()->FindBin(eta));
}


///////////////////////////////////////////////                                                                                                                                                           
//  Muon Trigger Correction 94X  pt-eta binned                                                                                                                                                           
////////////////////////////////////////////////////////////

float Cor94X_Trigger_Mu_EtaPt(float pt,float eta, TH2F* HistoTrg ){
  if (pt > 1000) pt=1000;
  return HistoTrg->GetBinContent(HistoTrg->GetXaxis()->FindBin(pt),HistoTrg->GetYaxis()->FindBin(fabs(eta)));
}



/////////////////////////////////////////////////////
//  Muon TRK Correction 94X
////////////////////////////////////////////////////////////
float Cor94X_TRK_Mu_Full2016(float eta, TGraphAsymmErrors * graph ) {
  // take the ratio_eff_eta3_dr030e030_corr histogram (as function of eta).
  double * ipoint=NULL;
  double x_point,y_sf;
  ipoint=graph->GetX();
  if (eta >=(ipoint[0]-graph->GetErrorXlow(0))  && eta < (ipoint[0]+graph->GetErrorXhigh(0)) )  etaPOINT=0 ;
  else if (eta >= (ipoint[1]-graph->GetErrorXlow(1))  && eta < (ipoint[1]+graph->GetErrorXhigh(1)) ) etaPOINT=1;
  else if (eta >= (ipoint[2]-graph->GetErrorXlow(2)) && eta < (ipoint[2]+graph->GetErrorXhigh(2)) ) etaPOINT=2;
  else if (eta >= (ipoint[3]-graph->GetErrorXlow(3)) && eta < (ipoint[3]+graph->GetErrorXhigh(3)) ) etaPOINT=3;
  else if (eta >= (ipoint[4]-graph->GetErrorXlow(4)) && eta < (ipoint[4]+graph->GetErrorXhigh(4)) ) etaPOINT=4;
  else if (eta >= (ipoint[5]-graph->GetErrorXlow(5)) && eta < (ipoint[5]+graph->GetErrorXhigh(5)) ) etaPOINT=5;
  else if (eta >= (ipoint[6]-graph->GetErrorXlow(6)) && eta < (ipoint[6]+graph->GetErrorXhigh(6)) ) etaPOINT=6;
  else if (eta >= (ipoint[7]-graph->GetErrorXlow(7)) && eta < (ipoint[7]+graph->GetErrorXhigh(7)) ) etaPOINT=7;
  else if (eta >= (ipoint[8]-graph->GetErrorXlow(8)) && eta < (ipoint[8]+graph->GetErrorXhigh(8)) ) etaPOINT=8;
  else if (eta >= (ipoint[9]-graph->GetErrorXlow(9)) && eta < (ipoint[9]+graph->GetErrorXhigh(9)) ) etaPOINT=9;
  else if (eta >= (ipoint[10]-graph->GetErrorXlow(10)) && eta < (ipoint[10]+graph->GetErrorXhigh(10)) ) etaPOINT=10;
  else if (eta >= (ipoint[11]-graph->GetErrorXlow(11)) && eta < (ipoint[11]+graph->GetErrorXhigh(11)) ) etaPOINT=11;
  else if (eta >= (ipoint[12]-graph->GetErrorXlow(12)) && eta < (ipoint[12]+graph->GetErrorXhigh(12)) ) etaPOINT=12;
  else if (eta >= (ipoint[13]-graph->GetErrorXlow(13)) && eta < (ipoint[13]+graph->GetErrorXhigh(13)) ) etaPOINT=13;
  else if (eta >= (ipoint[14]-graph->GetErrorXlow(14)) && eta < (ipoint[14]+graph->GetErrorXhigh(14)) ) etaPOINT=14;
  else return 1;
    
  graph->GetPoint(etaPOINT,x_point,y_sf);
  return y_sf;
    
}


//////////////////////////////////////////////////////////////
//  Electron Id/Iso Correction 74X  from HTT group
////////////////////////////////////////////////////////////

//https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
float Cor94X_IDIso_Ele(float pt, float eta,  TH2F * HistoEleSF){
    
  if (pt >= 10 && pt < 20 ) ptBIN=1;
  if (pt >= 20 && pt < 35 ) ptBIN=2;
  if (pt >= 35 && pt < 50) ptBIN=3;
  if (pt >= 50 && pt < 90) ptBIN=4;
  if (pt >= 90 && pt < 150) ptBIN=5;
  if (pt >= 150 ) ptBIN=6;
    
  if (eta >= -2.5 && eta < -2 ) etaBIN=1;
  if (eta >= -2 && eta < -1.566 ) etaBIN=2;
  if (eta >= -1.566 && eta < -1.444) etaBIN=3;
  if (eta >= -1.444 && eta < -0.800) etaBIN=4;
  if (eta >= -0.800 && eta < 0 ) etaBIN=5;
  if (eta >= 0 && eta < 0.800 ) etaBIN=6;
  if (eta >= 0.800 && eta < 1.444 ) etaBIN=7;
  if (eta >= 1.444 && eta < 1.566 ) etaBIN=8;
  if (eta >= 1.566 && eta < 2 ) etaBIN=9;
  if (eta >= 2 && eta < 2.5 ) etaBIN=10;
    

  //std::cout<<"coming inside Cor94X_IDIso_Ele() "<<std::endl;      
  //    float SF_0p5= HistoEleSF0p5->GetBinContent(etaBIN, ptBIN);
  //    float SF_5= HistoEleSF5->GetBinContent(etaBIN, ptBIN);
  //    float FinalSF= 0.05 * SF_0p5  + 0.95 * SF_5;   //approximation of 10/fb data
    
  //    cout << pt << "  " << eta << " "<<SF_0p5<<"  "<<SF_5 << "  "<< FinalSF <<"\n";
  //    cout<< "--->Electron    pt= "<<pt<<  "   eta" <<eta<< "  SF="<<HistoEleSF->GetBinContent(etaBIN, ptBIN)<<"\n";
  return HistoEleSF->GetBinContent(etaBIN, ptBIN);
  
}


float getCorrFactorMuon94X(bool isData, float pt, float eta, TH2F ** HistoId, TH2F ** HistoIso,TH2F ** HistoTrg, TGraphAsymmErrors * graph){
  if (isData)
    return 1;
  else{
        
    float Weighted_IDSF=Cor94X_ID_Mu(pt,eta,HistoId[0]);

    float Weighted_IsoSF=Cor94X_Iso_Mu(pt,eta,HistoIso[0]);

    //float Weighted_TriggerSF=Cor94X_Trigger_Mu_onlyEta(eta,HistoTrg[0]);
    float Weighted_TriggerSF=Cor94X_Trigger_Mu_EtaPt(pt,eta,HistoTrg[0]);
    float Tracking_SF=Cor94X_TRK_Mu_Full2016(eta, graph);
                             
    return (Weighted_IDSF * Weighted_IsoSF * Tracking_SF * Weighted_TriggerSF);
    // std::cout<<"coming inside  getCorrFactorMuon94X  "<<std::endl;     
 }
    
}



float GetBJetSF(bool isData, float x, float jetEta, float jetHadFlvr){
    
  //    cout<< "--->Btag    pt= "<<x<<  "   jetHadFlvr" <<jetHadFlvr<< "  SF="<<1.05636+0.000920353*x+-7.85916e-07*x*x+1.92221e-11*x*x*x <<"   "<< 0.931535+(1.40704e-05*x) <<"\n";
    
  if (isData) return 1;
  else {
    if (jetHadFlvr ==0 ){
      //            cout <<"GetBJetSF jetHadFlvr-->0  pt is"<<x<< " SF= "<<0.943355+8.95816/(x*x)+0.000240703*x<<"\n";
      return 0.943355+8.95816/(x*x)+0.000240703*x;
    }
    else if (jetHadFlvr ==4 || jetHadFlvr ==5 ){
      //            cout <<"GetBJetSF jetHadFlvr-->4/5  pt is"<<x<< " SF= "<<0.91423*((1.+(0.00958053*x))/(1.+(0.010132*x)))<<"\n";
      return 0.91423*((1.+(0.00958053*x))/(1.+(0.010132*x)));
    }
    else
      return 1;
    //  std::cout<<"coming inside GetBJetSF "<<std::endl;    
  }
}



// Loose WP of MVA Ele
float getCorrFactorMVA90WPElectron94X(bool isData, float pt, float eta,    TH2F * HistoEleSF ){
  
  if (isData)
    return 1;
  else
    return Cor94X_IDIso_Ele(pt,eta,HistoEleSF);
}
