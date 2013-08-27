class Pi0PhotGen : public BaseGen {
 protected:
  BasePart Photon, Target, Pi0, Recoil, Decay1, Decay2;
  Bool_t Incoh;
  Float_t Ang[3][181], Unp[3][181], TOb[3][181], EOb[3][181], FOb[3][181];
  TH2F *hPvP1, *hPvP1Lo, *hPvP1Hi, *hPvP2, *hPvP2Lo, *hPvP2Hi, *hMiM;
 public:
  Pi0PhotGen(TString, TString, Float_t, Float_t, Float_t);
  ~Pi0PhotGen();
  void Init(Float_t, Float_t);
  Bool_t NewEvent(Float_t);
  void Reset();
  void SaveHists(TString);
};

Pi0PhotGen::Pi0PhotGen(TString name, TString target, Float_t recomass, Float_t beamlo, Float_t beamhi) : BaseGen(name, beamlo, beamhi), Photon("Pi0Phot"), Target("Pi0Targ"), Pi0("Pi0Pi0"), Recoil("Pi0Reco"), Decay1("Pi0Dec1"), Decay2("Pi0Dec2") {

  cout << ("Constructing generator for "+ProcName) << endl;

  TargName = target;

  if(name.Contains("Incoh")) Incoh = kTRUE;

  else Incoh = kFALSE;

  Photon.Mass = 0;
  Target.Mass = recomass;
  Pi0.Mass = kMPI0_MEV;
  Recoil.Mass = recomass;
  Decay1.Mass = 0;
  Decay2.Mass = 0;

  Int_t ptag[4] = {14,7,1,1};

  InitNtuple(4,ptag);
  
  hPvP1 = new TH2F("hPvP1","Proton vs Decay Photon1 Angle",181,0,181,181,0,181);
  hPvP1->GetXaxis()->SetTitle("Photon Angle (deg)");
  hPvP1->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");

  hPvP1Lo = new TH2F("hPvP1Lo","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvP1Lo->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvP1Lo->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");
  hPvP1Lo->SetMarkerColor(2);

  hPvP1Hi = new TH2F("hPvP1Hi","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvP1Hi->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvP1Hi->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");

  hPvP2 = new TH2F("hPvP2","Proton vs Decay Photon2 Angle",181,0,181,181,0,181);
  hPvP2->GetXaxis()->SetTitle("Photon Angle (deg)");
  hPvP2->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");

  hPvP2Lo = new TH2F("hPvP2Lo","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvP2Lo->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvP2Lo->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");
  hPvP2Lo->SetMarkerColor(2);

  hPvP2Hi = new TH2F("hPvP2Hi","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvP2Hi->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvP2Hi->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");

  hMiM = new TH2F("hMissM","Missing Mass vs Proton K",201,0,201,201,0,201);
  hMiM->GetXaxis()->SetTitle("Proton Kinetic Energy (MeV)");
  hMiM->GetYaxis()->SetTitle("Missing Mass (MeV)");

};

Pi0PhotGen::~Pi0PhotGen(){

  cout << ("Deleting generator for "+ProcName) << endl;

};

void Pi0PhotGen::Init(Float_t PolG, Float_t PolT){

  cout << "Loading Pi0 Photoproduction data files" << endl;
  
  TString Text[3] = {"par/Pi0_cm_", "par/Pi0_cm_", "par/Pi0_cm_"};
  Int_t BeamE[3] = {BeamLo, BeamMd, BeamHi};

  const char *File[3];
  Int_t AngN, PhiN;

  Float_t val, max;

  for(Int_t i=0; i<3; i++){   
    Text[i] += BeamE[i];
    File[i] = Text[i];
    AngN = 0;
    max = 0;
    
    ifstream fin(File[i]);
    while(!fin.eof()){
      fin>>Ang[i][AngN];
      fin>>Unp[i][AngN];
      fin>>TOb[i][AngN];
      fin>>EOb[i][AngN];
      fin>>FOb[i][AngN];
      
      AngN++;
    }
    fin.close();
    
    for(AngN=0; AngN<181; AngN++){
      for(PhiN=0; PhiN<181; PhiN++){
	val = (Unp[i][AngN]*(1+PolT*TOb[i][AngN]*sin(-PhiN*kD2R)+PolG*PolT*FOb[i][AngN]*cos(-PhiN*kD2R)));
	CrossSec->Fill(BeamE[i], AngN, PhiN, val);
	if(val>max) max = val;
      }
    }
    CrossMax->Fill(BeamE[i], max);
  }

  FillCross();

  cout << "--------------------------------------------------" << endl << endl;
  cout << "Running" << endl << endl;

};

Bool_t Pi0PhotGen::NewEvent(Float_t BeamE){

  Reset();
  NewVertex();

  Photon.SetP4Lab(BeamE,BeamE,0,0);
  Target.SetP4Lab(Target.Mass,0,0,0);

  Collision2B(Photon, Target, Pi0, Recoil);

  if(Weight(BeamE,Pi0.Theta,Pi0.Phi)){
    return kFALSE;
  }
  
  if(Incoh){
    SpecModel(Target);
    Collision2B(Photon, Target, Pi0, Recoil);    
  }

  if(Pi0.Ener<Pi0.Mass){
    return kFALSE;
  }

  Decay2B(Pi0, Decay1, Decay2);

  Photon.HistLab();
  Target.HistLab();
  Pi0.HistLab();
  Recoil.HistLab();
  Decay1.HistLab();
  Decay2.HistLab();

  Photon.FillNtuple(var,3);
  Recoil.FillNtuple(var,8);
  Pi0.FillNtuple(var,13);
  Decay1.FillNtuple(var,18);
  Decay2.FillNtuple(var,23);
  
  // Fill ntuple
  h1->Fill(var);

  hPvP1->Fill(Decay1.Theta,Recoil.Theta);
  hPvP2->Fill(Decay2.Theta,Recoil.Theta);
  if(Recoil.KEner<40){
    hPvP1Lo->Fill(Decay1.Theta,Recoil.Theta);
    hPvP2Lo->Fill(Decay2.Theta,Recoil.Theta);
  }
  else{
    hPvP1Hi->Fill(Decay1.Theta,Recoil.Theta);
    hPvP2Hi->Fill(Decay2.Theta,Recoil.Theta);
  }

  hMiM->Fill(Recoil.KEner,(Photon.KEner-(Pi0.KEner+Recoil.KEner)));

  return kTRUE;
  
};

void Pi0PhotGen::Reset(){

  Photon.Reset();
  Target.Reset();
  Pi0.Reset();
  Recoil.Reset();
  Decay1.Reset();
  Decay2.Reset();
  
};

void Pi0PhotGen::SaveHists(TString file){

  cout << ("Saving "+ProcName+" histograms") << endl;

  TFile hfile(file, "RECREATE", "MC_Hists_File");

  Photon.WriteHists();
  Target.WriteHists();
  Pi0.WriteHists();
  Recoil.WriteHists();
  Decay1.WriteHists();
  Decay2.WriteHists();

  TCanvas *cPvP1 = new TCanvas("cPvP1", "Proton vs Photon Angle", 800, 1000);

  TLine *xlo = new TLine(20,0,20,180);
  TLine *xhi = new TLine(160,0,160,180);
  TLine *ylo = new TLine(0,20,180,20);
  TLine *yhi = new TLine(0,160,180,160);

  hPvP1Lo->Draw();
  hPvP1Hi->Draw("same");
  xlo->Draw("same");
  xhi->Draw("same");
  ylo->Draw("same");
  yhi->Draw("same");

  TCanvas *cPvP2 = new TCanvas("cPvP2", "Proton vs Photon Angle", 800, 1000);

  hPvP2Lo->Draw();
  hPvP2Hi->Draw("same");
  xlo->Draw("same");
  xhi->Draw("same");
  ylo->Draw("same");
  yhi->Draw("same");

  hPvP1Lo->Write();
  hPvP1Hi->Write();

  hPvP2Lo->Write();
  hPvP2Hi->Write();

  cPvP1->Write();
  cPvP2->Write();

  hPvP1->Write();
  hPvP2->Write();
  hMiM->Write();

  hfile.Close();

};
