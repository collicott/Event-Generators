class ComptonGen : public BaseGen {
 protected:
  BasePart Photon, Target, Scatter, Recoil;
  Bool_t Incoh;
  Float_t Ang[3][181], Unp[3][181], Pos[3][181], Min[3][181], Asy[3][181];
  TH2F *hPvP, *hPvPLo, *hPvPHi, *hMiM;

 public:
  ComptonGen(TString, TString, Float_t, Float_t, Float_t);
  ~ComptonGen();
  void Init(Float_t, Float_t);
  Bool_t NewEvent(Float_t);
  void Reset();
  void SaveHists(TString);
};

ComptonGen::ComptonGen(TString name, TString target, Float_t recomass, Float_t beamlo, Float_t beamhi) : BaseGen(name, beamlo, beamhi), Photon("CompPhot"), Target("CompTarg"), Scatter("CompScat"), Recoil("CompReco") {

  cout << ("Constructing generator for "+ProcName) << endl;

  TargName = target;

  if(name.Contains("Incoh")) Incoh = kTRUE;

  else Incoh = kFALSE;
  
  Photon.Mass = 0;
  Target.Mass = recomass;
  Scatter.Mass = 0;
  Recoil.Mass = recomass;

  Int_t ptag[2] = {14,1};

  InitNtuple(2,ptag);
  
  hPvP = new TH2F("hPvP","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvP->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvP->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");

  hPvPLo = new TH2F("hPvPLo","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvPLo->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvPLo->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");
  hPvPLo->SetMarkerColor(2);

  hPvPHi = new TH2F("hPvPHi","Proton vs Photon Angle",181,0,181,181,0,181);
  hPvPHi->GetXaxis()->SetTitle("Photon Recoil Angle (deg)");
  hPvPHi->GetYaxis()->SetTitle("Proton Scatter Angle (deg)");

  hMiM = new TH2F("hMissM","Missing Mass vs Proton K",201,0,201,201,0,201);
  hMiM->GetXaxis()->SetTitle("Proton Kinetic Energy (MeV)");
  hMiM->GetYaxis()->SetTitle("Missing Mass (MeV)");
};

ComptonGen::~ComptonGen(){

  cout << ("Deleting generator for "+ProcName) << endl;

};

void ComptonGen::Init(Float_t PolG, Float_t PolT){

  cout << "Loading Compton data files" << endl;
  
  TString Text[3] = {"par/Comp_lab_", "par/Comp_lab_", "par/Comp_lab_"};
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
      fin>>Ang[i][180-AngN];
      fin>>Unp[i][180-AngN];
      fin>>Pos[i][180-AngN];
      fin>>Min[i][180-AngN];
     
      AngN++;
    }
    fin.close();
    
    for(AngN=0; AngN<181; AngN++){
      Asy[i][AngN]=(PolG*PolT*((Pos[i][AngN]-Min[i][AngN])/(Pos[i][AngN]+Min[i][AngN])));
      for(PhiN=0; PhiN<181; PhiN++){
	val = (Unp[i][AngN]*(1+Asy[i][AngN]*cos(PhiN*kD2R)));
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

Bool_t ComptonGen::NewEvent(Float_t BeamE){

  Reset();
  NewVertex();

  Photon.SetP4Lab(BeamE,BeamE,0,0);
  Target.SetP4Lab(Target.Mass,0,0,0);

  Collision2B(Photon, Target, Scatter, Recoil);

  if(Weight(BeamE,Scatter.Theta,Scatter.Phi)){
    return kFALSE;
  }
  
  if(Incoh){
    SpecModel(Target);
    Collision2B(Photon, Target, Scatter, Recoil);    
  }

  Photon.HistLab();
  Target.HistLab();
  Scatter.HistLab();
  Recoil.HistLab();

  Photon.FillNtuple(var,3);
  Recoil.FillNtuple(var,8);
  Scatter.FillNtuple(var,13);
  
  // Fill ntuple
  h1->Fill(var);

  hPvP->Fill(Scatter.Theta,Recoil.Theta);
  if(Recoil.KEner<40) hPvPLo->Fill(Scatter.Theta,Recoil.Theta);
  else hPvPHi->Fill(Scatter.Theta,Recoil.Theta);

  hMiM->Fill(Recoil.KEner,(Photon.KEner-(Scatter.KEner+Recoil.KEner)));

  return kTRUE;
  
};

void ComptonGen::Reset(){

  Photon.Reset();
  Target.Reset();
  Scatter.Reset();
  Recoil.Reset();
  
};

void ComptonGen::SaveHists(TString file){

  cout << ("Saving "+ProcName+" histograms") << endl;

  TFile hfile(file, "RECREATE", "MC_Hists_File");

  Photon.WriteHists();
  Target.WriteHists();
  Scatter.WriteHists();
  Recoil.WriteHists();

  hPvPLo->Write();
  hPvPHi->Write();

  TCanvas *cPvP = new TCanvas("cPvP", "Proton vs Photon Angle", 800, 1000);

  TLine *xlo = new TLine(20,0,20,180);
  TLine *xhi = new TLine(160,0,160,180);
  TLine *ylo = new TLine(0,20,180,20);
  TLine *yhi = new TLine(0,160,180,160);

  TArrow ar1(30,30,10,10,0.01,"|>");
  TArrow ar2(90,30,90,10,0.01,"|>");
  TArrow ar3(150,30,170,10,0.01,"|>");
  TArrow ar4(30,90,10,90,0.01,"|>");
  TArrow ar6(150,90,170,90,0.01,"|>");
  TArrow ar7(30,150,10,170,0.01,"|>");
  TArrow ar8(90,150,90,170,0.01,"|>");
  TArrow ar9(150,150,170,170,0.01,"|>");

  ar1.SetLineWidth(2);
  ar1.SetLineColor(4);
  ar1.SetFillColor(4);
  ar2.SetLineWidth(2);
  ar2.SetLineColor(4);
  ar2.SetFillColor(4);
  ar3.SetLineWidth(2);
  ar3.SetLineColor(4);
  ar3.SetFillColor(4);
  ar4.SetLineWidth(2);
  ar4.SetLineColor(4);
  ar4.SetFillColor(4);
  ar6.SetLineWidth(2);
  ar6.SetLineColor(4);
  ar6.SetFillColor(4);
  ar7.SetLineWidth(2);
  ar7.SetLineColor(4);
  ar7.SetFillColor(4);
  ar8.SetLineWidth(2);
  ar8.SetLineColor(4);
  ar8.SetFillColor(4);
  ar9.SetLineWidth(2);
  ar9.SetLineColor(4);
  ar9.SetFillColor(4);

  TPaveText pt1(30,30,50,40);
  TPaveText pt2(80,30,100,40);
  TPaveText pt3(130,30,150,40);
  TPaveText pt4(30,85,50,95);
  TPaveText pt5(80,85,100,95);
  TPaveText pt6(130,85,150,95);
  TPaveText pt7(30,140,50,150);
  TPaveText pt8(80,140,100,150);
  TPaveText pt9(130,140,150,150);

  pt1.AddText("Test1");
  pt2.AddText("Test2");
  pt3.AddText("Test3");
  pt4.AddText("Test4");
  pt5.AddText("Test5");
  pt6.AddText("Test6");
  pt7.AddText("Test7");
  pt8.AddText("Test8");
  pt9.AddText("Test9");

  hPvPLo->SetStats(kFALSE);
  hPvPLo->Draw();
  hPvPHi->Draw("same");
  xlo->Draw("same");
  xhi->Draw("same");
  ylo->Draw("same");
  yhi->Draw();
  ar1.Draw();
  ar2.Draw();
  ar3.Draw();
  ar4.Draw();
  ar6.Draw();
  ar7.Draw();
  ar8.Draw();
  ar9.Draw();
  pt1.Draw();
  pt2.Draw();
  pt3.Draw();
  pt4.Draw();
  pt5.Draw();
  pt6.Draw();
  pt7.Draw();
  pt8.Draw();
  pt9.Draw();

  cPvP->Write();

  hPvP->Write();
  hMiM->Write();

  hfile.Close();

};
