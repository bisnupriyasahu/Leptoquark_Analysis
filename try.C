#include <string>
#include <iostream>
#include <vector>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <vector>
#include "TFile.h"

int main(int argc, char** argv) {

  using namespace std;

  std::string out = *(argv + 1);
  std::cout << "\n\n\n OUTPUT NAME IS:    " << out << std::endl;
  TFile *fout = TFile::Open(out.c_str(), "RECREATE");

  std::vector<string> input;
  for (int f = 2; f < argc; f++) {
    input.push_back(*(argv + f));
    std::cout <<"INPUT NAME IS:   " << input[f - 2] << "\n";
  }



}
