/*
 * EvGenBasic
 *
 * This routine is the event generator for GEANT4.  It outputs a root file
 * with an ntuple with the appropriate variable and names.
 *
 * 27.03.2010		DLH		First Version adapted from EvGenRes
 *
 */

// This ifndef section allows you to use this code from the CINT command
// line or as stand-alone executable code.
#ifndef __CINT__

	#include "TF1.h"
	#include "TMath.h"
	#include "TH1F.h"
	#include "TH2F.h"
	#include "TNtuple.h"
	#include "TFile.h"
	#include "TString.h"
	#include "TRandom.h"
	#include "TVector3.h"
	#include "TLorentzVector.h"
	#include <fstream>

	Double_t Sqr( Double_t);
	Double_t Momentum( Double_t, Double_t);
	int EvGenBasic();
	TString GenNames( Int_t, Int_t*);

	int main()
	{
	  return EvGenBasic();
	}

#endif

// Some physical constants and functions.
#include "physics.h"

//
// EvGenBasic
//
// The main event generator code.
//
int EvGenBasic() 
{
	Int_t i, counts;
	Int_t npart, update;
	Int_t ptag[10];
	
	Double_t tgt_len, vtx_rad, sig_bm;
	Double_t e_lo, e_hi, e_mid, e_step;
	Double_t th_lo, th_hi, th_mid, th_step;
	Double_t ph_lo, ph_hi, ph_mid, ph_step;
	Double_t KE, energy, momentum, theta, phi;
	Double_t vtx_x, vtx_y, vtx_z;
	Double_t pm;

	TString particle, string, name, gnames;

	TVector3 vtx, dircos;
	TLorentzVector p;

	std::cout << "--------" << std::endl;
	std::cout << "EvGenBasic" << std::endl;
	std::cout << "--------" << std::endl;

	// Default parameters
	particle = "g";
	tgt_len = 2.00;
	vtx_rad = 0.50;
	counts = (int) 1e5;
	e_lo = 150;
	e_hi = 200;
	e_step = 10;
	th_lo = 10;
	th_hi = 170;
	th_step = 20;
	ph_lo = -180;
	ph_hi = 180;
	ph_step = 20;

	// Read in parameters from parameter file
	name ="par/EvGenBasic.in";
	ifstream inFile( name);
	if ( !inFile.is_open()) 
	{
		std::cout << "Error opening file ";
		std::cout << name;
		std::cout << std::endl;
		exit( -1);
	}
	while( !inFile.eof()) {
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			string = "Particle: ";
			if ( name.Contains( string))
				particle = name.Remove( 0, string.Length());
			string = "TargetLength: ";
			if ( name.Contains( string)) {
				name.Remove( 0, string.Length());
				tgt_len = name.Atof();
			}
			string = "BeamSpotRadius: ";
			if ( name.Contains( string)) {
				name.Remove( 0, string.Length());
				vtx_rad = name.Atof();
			}
			string = "Throws: ";
			if ( name.Contains( string)) {
				name.Remove( 0, string.Length());
				counts = name.Atoi();
			}
			string = "ParticleEnergy: ";
			if ( name.Contains( string)) {
				name.Remove( 0, string.Length());
				string = name;
				string.Remove( 0 , string.Last(' '));
				e_step = string.Atof();

				string = name;
				string.Remove( string.First(' '));
				e_lo = string.Atof();

				string = name( name.First(' '), name.Last(' ')-name.First(' '));
				e_hi = string.Atof();
			}
			string = "LabTheta: ";
			if ( name.Contains( string)) {
				name.Remove( 0, string.Length());
				string = name;
				string.Remove( 0 , string.Last(' '));
				th_step = string.Atof();

				string = name;
				string.Remove( string.First(' '));
				th_lo = string.Atof();

				string = name( name.First(' '), name.Last(' ')-name.First(' '));
				th_hi = string.Atof();
			}
			string = "LabPhi: ";
			if ( name.Contains( string)) {
				name.Remove( 0, string.Length());
				string = name;
				string.Remove( 0 , string.Last(' '));
				ph_step = string.Atof();

				string = name;
				string.Remove( string.First(' '));
				ph_lo = string.Atof();

				string = name( name.First(' '), name.Last(' ')-name.First(' '));
				ph_hi = string.Atof();
			}
		}
	}
	inFile.close();

	if ( ( particle != "g") && ( particle != "p")) {
		std::cout << "Particle \"" << particle << "\" not supported";
		std::cout << std::endl;
		exit( -1);
	}

	name = "Particle = " + particle;
	std::cout << name << std::endl;
	name = Form( "Target Length = %4.2f cm", tgt_len);
	std::cout << name << std::endl;
	name = Form( "Beam Spot Radius = %4.2f cm", vtx_rad);
	std::cout << name << std::endl;
	name = Form( "Throws = %d", counts);
	std::cout << name << std::endl;
	name = Form( "Particle Energy Range = %5.1f - %5.1f MeV in %4.1f MeV steps",
			e_lo, e_hi, e_step);
	std::cout << name << std::endl;
	name = Form( "Particle Theta Range = %5.1f - %5.1f deg in %4.1f deg steps",
			th_lo, th_hi, th_step);
	std::cout << name << std::endl;
	name = Form( "Particle Phi Range = %6.1f - %5.1f deg in %4.1f deg steps",
			ph_lo, ph_hi, ph_step);
	std::cout << name << std::endl;

	update = counts/20;

	// Particle ID #s
	// They are standard geant numbers

	// One particle
	npart = 1;

	// Default particle is a photon.
	ptag[0] = 1;
	pm = 0;

	// If proton
	if ( particle == "p") {
		ptag[0] = 14;
		pm = kMP_MEV/1000;
	}

	// Array for filling ntuple
	Float_t var[100];

	//	Width of gaussian for beam profile on target in cm
	sig_bm = 0.5;

	// Set the seed for the random number generator
	gRandom->SetSeed();

	//	Generate GEANT name string for ntuple
	gnames = GenNames( npart, ptag);

	// Loop over energy or angle
	for ( e_mid = e_lo; e_mid <= e_hi; e_mid += e_step) {

		name = Form( " Particle Energy = %5.1f MeV", e_mid);
		std::cout << name << std::endl;
		KE = e_mid/1000;

		for ( th_mid = th_lo; th_mid <= th_hi; th_mid += th_step) {

			name = Form( " Particle Theta = %5.1f MeV", th_mid);
			std::cout << name << std::endl;
			theta = th_mid*kD2R;

			for ( ph_mid = ph_lo; ph_mid <= ph_hi; ph_mid += ph_step) {

				name = Form( " Particle Phi = %5.1f MeV", ph_mid);
				std::cout << name << std::endl;
				phi = ph_mid*kD2R;

				// Filename
				name = Form( "out/basic/%s_%d_%d_%d_in.root",
						(const char*) particle, (int) e_mid, (int) th_mid,
						(int) ph_mid);
				TFile hfile( name, "RECREATE", "MC_Ntuple_File");

				//	Create ntuple for kinematic variables
				//	It is absolutely necessary for Geant4!
				TNtuple *h1 = new TNtuple( "h1", "TMCUserGenerator", gnames);

				// These histograms are only for debugging.  You can comment them
				// out if you want.
				TH1F *h2 = new TH1F( "h2", "Particle KE (MeV)", 1000, 0, 1000);
				TH1F *h3 = new TH1F( "h3", "Particle Momentum (MeV/c)", 2000, 0,
						2000);
				TH1F *h4 = new TH1F( "h4", "Particle #theta (deg)", 180, 0, 180);
				TH1F *h5 = new TH1F( "h5", "Particle #phi (deg)", 360, -180, 180);

				for ( i = 1; i <= counts; i++)
				{
					if ( i && (i%update) == 0)
						std::cout << "     events analysed: " << i << std::endl;

					vtx_z = tgt_len*(-0.5 + gRandom->Rndm());
					while ( sqrt( Sqr( vtx_x = gRandom->Gaus(0,sig_bm))
									+ Sqr( vtx_y = gRandom->Gaus(0,sig_bm)))
									> vtx_rad);
					vtx.SetXYZ( vtx_x, vtx_y, vtx_z);

					// Interaction vertex position.
					var[0] = vtx.X();
					var[1] = vtx.Y();
					var[2] = vtx.Z();

					// Incident photon beam, in this case turned off.
					var[3] = 0;
					var[4] = 0;
					var[5] = 0;
					var[6] = 0;
					var[7] = 0;

					// Direction cosines
					dircos.SetXYZ( sin( theta)*cos( phi), sin( theta)*sin( phi),
							cos( theta));

					// Energy and momentum
					energy = KE + pm;
					momentum = Momentum( energy, pm);

					// Particle variables
					var[8] = dircos.X();
					var[9] = dircos.Y();
					var[10] = dircos.Z();
					var[11] = momentum;
					var[12] = energy;

					// Fill ntuple
					h1->Fill( var);

					// Fill histograms with quantities (in MeV)
					h2->Fill( e_mid);
					h3->Fill( momentum*1000);
					h4->Fill( th_mid);
					h5->Fill( ph_mid);
				}

				// Write histograms to file
				hfile.Write();

				// This isn't really necessary, but can be used for debugging.
//				h1->Print();

				// Close file
				hfile.Close();
			}
		}
	}

	return 0;
}

// Calculates relativistic momentum from mass and energy.
Double_t Momentum( Double_t en, Double_t m)
{
	if ( en >= m) return( sqrt( en*en - m*m));
	else return( -1);
}

// Calculates the square.
Double_t Sqr( Double_t x)
{
	return( x*x);
}

TString GenNames( Int_t npart, Int_t* ptag) {

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

	return( names);
}
