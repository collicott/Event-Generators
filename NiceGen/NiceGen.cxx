#ifndef __CINT__

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TLine.h"
#include "TArrow.h"
#include "TPaveText.h"
#include <fstream>
using namespace std;

// Some physical constants.
#include "physics.h"
// Base Particle class
#include "BasePart.h"
// Base Generator class
#include "BaseGen.h"
// Pi0 Generator class
#include "Pi0PhotGen.h"
// Compton Generator class
#include "ComptonGen.h"

int main()
{
  gRandom->SetSeed(0);

  Int_t i=0, j=0;

  Float_t PolG = 0.7, PolT = 0.8;
  Float_t BeamLo, BeamHi;

  TString proc[2] = {"Compton","Pi0 Photoproduction"};
  TString type[3] = {"Norm","Incoh","Coher"};
  TString targ[6] = {"p","n","3He","4He","12C","16O"};
  //TString reco[6] = {"p","n","3He","4He","12C","16O"};
  TString name;

  Float_t recomass[6] = {kMP_MEV,kMN_MEV,kM_HE3_MEV,kM_HE4_MEV,kM_C12_MEV,kM_O16_MEV};
  
  Int_t procsel, typesel, targsel, recosel, NEvn;
  
  cout << "Choose process:" << endl;
  cout << "1) Compton" << endl;
  cout << "2) Pi0 Photoproduction" << endl;
  cout << "3) Random" << endl;
  cout << "--------------------" << endl;

  cin >> procsel;
  procsel--;

  cout << "Choose type:" << endl;
  cout << "1) Normal" << endl;
  cout << "2) Incoherent" << endl;
  cout << "3) Coherent" << endl;
  cout << "--------------------" << endl;

  cin >> typesel;
  typesel--;

  if(typesel==0) targsel = 0;
  else{
    cout << "Choose target:" << endl;
    cout << "1) 3He" << endl;
    cout << "2) 4He" << endl;
    cout << "3) 12C" << endl;
    cout << "4) 16O" << endl;
    cout << "--------------------" << endl;

    cin >> targsel;
    targsel++;

    PolT = 0;
  }
  
  if(typesel==1){
    cout << "Choose recoil:" << endl;
    cout << "1) p" << endl;
    cout << "2) n" << endl;
    cout << "--------------------" << endl;

    cin >> recosel;
    recosel--;
  }
  else recosel = targsel;

  cout << "Enter minimum tagged photon energy (MeV)" << endl;

  cin >> BeamLo;

  cout << "Enter maximum tagged photon energy (MeV)" << endl;

  cin >> BeamHi;

  cout << "Enter number of events:" << endl;

  cin >> NEvn;

  name = (type[typesel]+" "+proc[procsel]+" on "+targ[targsel]);

  cout << endl << (name+" selected") << endl << endl;
  cout << "--------------------------------------------------" << endl << endl;
  
  TF1 *f1 = new TF1("f1", "1/x", BeamLo, BeamHi);

  TString file1 = "out/hist.root";
  TString file2 = "out/ntpl.root";

  ComptonGen *cgen;
  Pi0PhotGen *pgen;
  
  if(proc[procsel]=="Compton"){
    cgen = new ComptonGen(name,targ[targsel],recomass[recosel],BeamLo,BeamHi);
    cgen->Init(PolG, PolT);

    while(i<NEvn){
      if(cgen->NewEvent(f1->GetRandom())) i++;
      else{
	j++;
	continue;
      }
    }
  }
  else if(proc[procsel]=="Pi0 Photoproduction"){
    pgen = new Pi0PhotGen(name,targ[targsel],recomass[recosel],BeamLo,BeamHi);
    pgen->Init(PolG, PolT);

    while(i<NEvn){
      if(pgen->NewEvent(f1->GetRandom())) i++;
      else{
	j++;
	continue;
      }
    }
  }
  else{
    cout << "Invalid process selection" << endl;
    return 0;
  }
  
  cout << i << " good events, " << j << " bad events." << endl << endl;
  cout << "--------------------------------------------------" << endl << endl;

  if(proc[procsel]=="Compton"){
    cgen->SaveHists(file1);
    cgen->SaveNtuple(file2);
    delete cgen;
  }
  else if(proc[procsel]=="Pi0 Photoproduction"){
    pgen->SaveHists(file1);
    pgen->SaveNtuple(file2);
    delete pgen;
  }

  return 0;
}

#endif
