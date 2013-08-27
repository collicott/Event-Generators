class BasePart {
 protected:
  TString PartName;
  void FillCM();
  void FillLab();
 public:
  BasePart(TString);
  ~BasePart();
  TLorentzVector P4, P4CM;
  Float_t Mass, Ener, KEner, Mom, Theta, Phi;
  Float_t EnerCM, KEnerCM, MomCM, ThetaCM, PhiCM;
  TH1F *hKEner, *hMom, *hTheta, *hPhi, *hMass;
  void BoostCM(TVector3);
  void BoostLab(TVector3);
  void HistCM();
  void HistLab();
  void Reset();
  void SetP4CM(Float_t, Float_t, Float_t, Float_t);
  void SetP4Lab(Float_t, Float_t, Float_t, Float_t);
  void SetP4CM(Float_t, Float_t);
  void SetP4Lab(Float_t, Float_t);
  TString WhichDet();
  void WriteHists();
  void FillNtuple(Float_t*, Int_t);
};

BasePart::BasePart(TString name){

  PartName = name;

  cout << ("Constructing "+PartName) << endl;

  hKEner = new TH1F(PartName+"KEner",PartName+" Kinetic E (MeV)",1000,0,300);
  hKEner->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  hKEner->GetYaxis()->SetTitle("# of Events");

  hMom = new TH1F(PartName+"Mom",PartName+" Momentum (MeV)",1000,0,500);
  hMom->GetXaxis()->SetTitle("Momentum (MeV)");
  hMom->GetYaxis()->SetTitle("# of Events");

  hTheta = new TH1F(PartName+"Theta",PartName+" Theta (deg)",1000,0,180);
  hTheta->GetXaxis()->SetTitle("Theta (deg)");
  hTheta->GetYaxis()->SetTitle("# of Events");

  hPhi = new TH1F(PartName+"Phi",PartName+" Phi (deg)",1000,-180,180);
  hPhi->GetXaxis()->SetTitle("Phi (deg)");
  hPhi->GetYaxis()->SetTitle("# of Events");

  hMass = new TH1F(PartName+"Mass",PartName+" Mass (MeV)",1000,0,1000);
  hMass->GetXaxis()->SetTitle("Mass (MeV)");
  hMass->GetYaxis()->SetTitle("# of Events");

};

BasePart::~BasePart(){

  cout << ("Deleting "+PartName) << endl;

  delete hKEner;
  delete hMom;
  delete hTheta;
  delete hPhi;
  delete hMass;

};

void BasePart::FillCM(){

  EnerCM = P4CM.E();
  KEnerCM = (EnerCM-Mass);
  MomCM = P4CM.Rho();
  ThetaCM = (P4CM.Theta()*kR2D);
  PhiCM = (P4CM.Phi()*kR2D);

};

void BasePart::FillLab(){

  Ener = P4.E();
  KEner = (Ener-Mass);
  Mom = P4.Rho();
  Theta = (P4.Theta()*kR2D);
  Phi = (P4.Phi()*kR2D);

};

void BasePart::BoostCM(TVector3 v){

  P4CM = P4;
  P4CM.Boost(v);
  FillCM();

};

void BasePart::BoostLab(TVector3 v){

  P4 = P4CM;
  P4.Boost(v);
  FillLab();

};

void BasePart::HistCM(){

  hKEner->Fill(KEnerCM);
  hMom->Fill(MomCM);
  hTheta->Fill(ThetaCM);
  hPhi->Fill(PhiCM);
  hMass->Fill(Mass);

};

void BasePart::HistLab(){

  hKEner->Fill(KEner);
  hMom->Fill(Mom);
  hTheta->Fill(Theta);
  hPhi->Fill(Phi);
  hMass->Fill(Mass);

};

void BasePart::Reset(){

  P4.SetPxPyPzE(0,0,0,0);
  P4CM.SetPxPyPzE(0,0,0,0);
  FillCM();
  FillLab();

};

void BasePart::SetP4CM(Float_t ener,Float_t mom,Float_t theta,Float_t phi){

  EnerCM = ener;
  KEnerCM = (Ener-Mass);
  MomCM = mom;
  ThetaCM = theta;
  PhiCM = phi;
  P4CM.SetPxPyPzE(1,1,1,EnerCM);
  P4CM.SetRho(MomCM);
  P4CM.SetTheta(ThetaCM*kD2R);
  P4CM.SetPhi(PhiCM*kD2R);

};

void BasePart::SetP4Lab(Float_t ener,Float_t mom,Float_t theta,Float_t phi){

  Ener = ener;
  KEner = (Ener-Mass);
  Mom = mom;
  Theta = theta;
  Phi = phi;
  P4.SetPxPyPzE(1,1,1,Ener);
  P4.SetRho(Mom);
  P4.SetTheta(Theta*kD2R);
  P4.SetPhi(Phi*kD2R);

};

void BasePart::SetP4CM(Float_t ener, Float_t mom){

  Float_t theta = acos(-1+2*gRandom->Rndm());
  Float_t phi = kPI*(-1+2*gRandom->Rndm());
  SetP4CM(ener, mom, (theta*kR2D), (phi*kR2D));

};

void BasePart::SetP4Lab(Float_t ener, Float_t mom){

  Float_t theta = acos(-1+2*gRandom->Rndm());
  Float_t phi = kPI*(-1+2*gRandom->Rndm());
  SetP4Lab(ener, mom, (theta*kR2D), (phi*kR2D));

};

TString BasePart::WhichDet(){
  TString Det;
  if(Theta<=20) Det = "TAPS";
  else if(Theta<=160) Det = "CB";
  else Det = "Out";

  return Det;
};

void BasePart::WriteHists(){

  hKEner->Write();
  hMom->Write();
  hTheta->Write();
  hPhi->Write();
  hMass->Write();

};

void BasePart::FillNtuple(Float_t* var, Int_t index){

  var[index++] = sin(Theta*kD2R)*cos(Phi*kD2R);
  var[index++] = sin(Theta*kD2R)*sin(Phi*kD2R);
  var[index++] = cos(Theta*kD2R);
  var[index++] = Mom/1000;
  var[index++] = Ener/1000;

};
