class BaseGen {
 protected:
  TString ProcName;
  TString TargName;
  Float_t BeamLo;
  Float_t BeamMd;
  Float_t BeamHi;
  Float_t BeamSt;
  TH3F *CrossSec;
  TH1F *CrossMax;
  TH1F *CrossTot;
  TLorentzVector ptot;
  TVector3 cm_to_lab, lab_to_cm;
  Float_t vtx_x, vtx_y, vtx_z;
 public:
  BaseGen(TString, Float_t, Float_t);
  ~BaseGen();
  void Collision2B(BasePart&, BasePart&, BasePart&, BasePart&);
  void Decay2B(BasePart&, BasePart&, BasePart&);
  void FillCross();
  void Interp1(TH1F*, Float_t, Float_t);
  void Interp3(TH3F*, Float_t, Float_t);
  void SpecModel(BasePart&);
  void InitNtuple(Int_t, Int_t*);
  void SaveNtuple(TString);
  void NewVertex();
  Float_t Sqr(Float_t x){return(x*x);};
  Bool_t Weight(Float_t, Float_t, Float_t);
  TNtuple *h1;
  Float_t var[100];
};

BaseGen::BaseGen(TString name, Float_t beamlo, Float_t beamhi){

  ProcName = name;
  BeamLo = beamlo;
  BeamMd = ((beamlo+beamhi)/2);
  BeamHi= beamhi;
  BeamSt = ((beamhi-beamlo)/20);

  cout << ("Constructing base generator for "+ProcName) << endl;

  CrossSec = new TH3F(ProcName,ProcName,21,BeamLo,(BeamHi+BeamSt),181,0,181,181,0,181);
  CrossMax = new TH1F(ProcName+"Max",(ProcName+" - max"),21,BeamLo,(BeamHi+BeamSt));
  CrossTot = new TH1F(ProcName+"Tot",(ProcName+" - tot"),21,BeamLo,(BeamHi+BeamSt));
  
};

BaseGen::~BaseGen(){

  cout << ("Deleting base generator for "+ProcName) << endl;

  delete CrossSec;
  delete CrossMax;
  delete CrossTot;

};

void BaseGen::Collision2B(BasePart& qi,BasePart& ki,BasePart& qf,BasePart& kf){

  Float_t ener, mom;
  
  ptot = qi.P4 + ki.P4;
  cm_to_lab = ptot.BoostVector();
  lab_to_cm = -ptot.BoostVector();
  
  qi.BoostCM(lab_to_cm);
  ki.BoostCM(lab_to_cm);
  
  ener = ((ptot.M2()+Sqr(qf.Mass)-Sqr(kf.Mass))/(2*ptot.M()));
  mom = sqrt(Sqr(ener)-Sqr(qf.Mass));

  if(qf.Ener==0){  
    qf.SetP4CM(ener,mom);
    qf.BoostLab(cm_to_lab);
  }

  else{
    qf.BoostCM(lab_to_cm);
    qf.SetP4CM(ener,mom,qf.ThetaCM,qf.PhiCM);
    qf.BoostLab(cm_to_lab);
  }

  ener = (qi.EnerCM + ki.EnerCM - qf.EnerCM);

  kf.SetP4CM(ener,mom,(180-qf.ThetaCM),-(180-qf.PhiCM));
  kf.BoostLab(cm_to_lab);
  
};

void BaseGen::Decay2B(BasePart& k,BasePart& p1,BasePart& p2){

  Float_t ener, mom;
  
  ptot = k.P4;
  cm_to_lab = ptot.BoostVector();
  lab_to_cm = -ptot.BoostVector();
  
  k.BoostCM(lab_to_cm);
  
  ener = ((ptot.M2()+Sqr(p1.Mass)-Sqr(p2.Mass))/(2*ptot.M()));
  mom = sqrt(Sqr(ener)-Sqr(p1.Mass));

  p1.SetP4CM(ener,mom);
  p1.BoostLab(cm_to_lab);

  p2.SetP4CM(ener,mom,(180-p1.ThetaCM),-(180-p1.PhiCM));
  p2.BoostLab(cm_to_lab);

};

void BaseGen::FillCross(){

  Int_t AngN, PhiN;
  Float_t SolAng = 0, SolAngTot = 0, CTotLo = 0, CTotMd = 0, CTotHi = 0;

  for(AngN=0; AngN<180; AngN++){
    SolAng = ((cos(AngN*kD2R)-cos((AngN+1)*kD2R))*kD2R);
    for(PhiN=0; PhiN<180; PhiN++){
      CTotLo += (2*(CrossSec->GetBinContent( 1, (AngN+1), (PhiN+1)))*SolAng);
      CTotMd += (2*(CrossSec->GetBinContent(11, (AngN+1), (PhiN+1)))*SolAng);
      CTotHi += (2*(CrossSec->GetBinContent(21, (AngN+1), (PhiN+1)))*SolAng);
      SolAngTot += (2*SolAng);
    }
  }

  cout << BeamLo << " MeV - " << CTotLo << " nb" << endl;
  cout << BeamMd << " MeV - " << CTotMd << " nb" << endl;
  cout << BeamHi << " MeV - " << CTotHi << " nb" << endl;

  CrossTot->Fill(BeamLo, CTotLo);
  CrossTot->Fill(BeamMd, CTotMd);
  CrossTot->Fill(BeamHi, CTotHi);

  Interp3(CrossSec, BeamLo, BeamMd);
  Interp3(CrossSec, BeamMd, BeamHi);

  Interp1(CrossMax, BeamLo, BeamMd);
  Interp1(CrossMax, BeamMd, BeamHi);

  Interp1(CrossTot, BeamLo, BeamMd);
  Interp1(CrossTot, BeamMd, BeamHi);

};

void BaseGen::Interp1(TH1F* hist, Float_t lo, Float_t hi){

  Float_t y1, y2, x1, x2, m, b, EnrN;

  y1 = hist->GetBinContent(((lo-BeamLo)/BeamSt)+1);
  y2 = hist->GetBinContent(((hi-BeamLo)/BeamSt)+1);
  x1 = lo;
  x2 = hi;
  m = ((y2-y1)/(x2-x1));
  b = (y1-m*x1);
  for(EnrN=(lo+BeamSt); EnrN<hi; EnrN+=BeamSt){
    hist->Fill(EnrN, (m*EnrN+b));
  }

};

void BaseGen::Interp3(TH3F* hist, Float_t lo, Float_t hi){

  Int_t AngN, PhiN;
  Float_t y1, y2, x1, x2, m, b, EnrN;
  
  for(AngN=0; AngN<181; AngN++){
    for(PhiN=0; PhiN<181; PhiN++){
      y1 = hist->GetBinContent((((lo-BeamLo)/BeamSt)+1), (AngN+1), (PhiN+1));
      y2 = hist->GetBinContent((((hi-BeamLo)/BeamSt)+1), (AngN+1), (PhiN+1));
      x1 = lo;
      x2 = hi;
      m = ((y2-y1)/(x2-x1));
      b = (y1-m*x1);
      for(EnrN=(lo+BeamSt); EnrN<hi; EnrN+=BeamSt){
	hist->Fill(EnrN, AngN, PhiN, (m*EnrN+b));
      }
    }
  }

};

void BaseGen::SpecModel(BasePart& ki){

  Float_t kpi, prob, func, ener, N_s, N_p;
  Float_t k_max, A[2], coeff[2], sigma[2], Emin[2], Emax[2], sf_max[2];
  Int_t L;

  if(TargName=="12C"){
    N_s = 0.52; N_p = 1.65;
    k_max = 300.0;
    A[0] = 43.5; A[1] = 13500.0;
    coeff[0] = 0.0; coeff[1] = 2.0;
    sigma[0] = 90.0; sigma[1] = 75.0;
    Emin[0] = 26.0; Emin[1] = 16.0;
    Emax[0] = 50.0; Emax[1] = 26.0;
    sf_max[0] = 0.26; sf_max[1] = 0.93;
  }

  if(TargName=="16O"){
    N_s = 0.52; N_p = 1.65;
    k_max = 300.0;
    A[0] = 43.5; A[1] = 13500.0;
    coeff[0] = 0.0; coeff[1] = 2.0;
    sigma[0] = 90.0; sigma[1] = 75.0;
    Emin[0] = 26.0; Emin[1] = 16.0;
    Emax[0] = 50.0; Emax[1] = 26.0;
    sf_max[0] = 0.26; sf_max[1] = 0.93;
  }

  if((gRandom->Rndm()) < N_s/(N_s+N_p)) L = 0;
  else L = 1;
  
  prob = 1;
  func = 0;

  while(prob>func){
    kpi = k_max*gRandom->Rndm();
    prob = sf_max[L]*gRandom->Rndm();
    func = A[L]*pow((kpi/1000.0),(coeff[L]+2))*exp(-Sqr(kpi/sigma[L])/2);
  }

  ener = (kMP_MEV-(Emin[L]+(Emax[L]-Emin[L])*gRandom->Rndm()));

  ki.SetP4Lab(ener, kpi);

};

Bool_t BaseGen::Weight(Float_t Enr, Float_t Ang, Float_t Phi){
  
  Bool_t Check = kFALSE;

  Float_t val = CrossSec->GetBinContent((((Enr-BeamLo)/BeamSt)+1), (Ang+1), (fabs(Phi)+1));
  Float_t max = CrossMax->GetBinContent(((Enr-BeamLo)/BeamSt)+1);

  if(val <= (max*gRandom->Rndm())) Check = kTRUE;

  return Check;

};

void BaseGen::InitNtuple(Int_t npart, Int_t* ptag) {
  
  Int_t i, j;
  
  TString pstr[] = {"Px", "Py", "Pz", "Pt", "En"};
  TString beam = "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm";
  TString particles;
  TString names;
  
  for ( i = 0; i < npart; i++) {
    for ( j = 0; j < 5; j++) {
      particles.Append( pstr[j]);
      if ( ( i == (npart-1)) && ( j == 4))
	particles.Append( Form( "_l%02d%02d", i+1, ptag[i]));
      else
	particles.Append( Form( "_l%02d%02d:", i+1, ptag[i]));
    }
  }
  
  names = beam + ":" + particles;

  h1 = new TNtuple("h1", "TMCUserGenerator", names);
  
}

void BaseGen::SaveNtuple(TString file){

  TFile hfile(file, "RECREATE", "MC_Ntuple_File");

  h1->Write();

  hfile.Close();
  
};

void BaseGen::NewVertex() {

    vtx_x = 0.5;
    vtx_y = 0.5;

    while(sqrt(Sqr(vtx_x)+Sqr(vtx_y)) > 0.5){
      vtx_x = gRandom->Gaus(0,0.5);
      vtx_y = gRandom->Gaus(0,0.5);
    }

    vtx_z = 2.0*(-0.5 + gRandom->Rndm());

    var[0] = vtx_x;
    var[1] = vtx_y;
    var[2] = vtx_z;

}
