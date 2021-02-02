#ifndef plotfill_histo_h
#define  plotfill_histo_h


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
#include "TSystem.h"
//#include "myevent.h"
//#include "LinkDef.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
using namespace std;

//****************************************************
map<string, TH1F*>* myMap1;
map<string, TH2F*>* myMap2;
//**********************************************

TH1F* nplot1(string name) {
  if (myMap1->find(name) != myMap1->end())
    return (*myMap1)[name];
  else
    return 0;
}

TH2F* nplot2(string name) {
  if (myMap2->find(name) != myMap2->end())
    return (*myMap2)[name];
  else
    return 0;
}
//****************************************************



void plotFill(string name, float x, int nx, float nxmin, float nxmax, double weight=1) {
  if (myMap1->find(name) == myMap1->end())
    (*myMap1)[name] = new TH1F(name.c_str(), name.c_str(), nx, nxmin, nxmax);
  (*myMap1)[name]->SetDefaultSumw2();
  (*myMap1)[name]->Fill(x,weight);
}




#endif
