//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(CarmTracking* _carm)
: G4VUserPrimaryGeneratorAction(), peak_energy(80), filter_thickness(1.0), frameNo(-1), isFirst(true)
, eIntensity(0), tubeVoltage(0), tubeCurrent(0), monitorDAP(0), monitorPower(true), monitorTime(1)
{
	fPrimary = new G4ParticleGun();
	fMessenger = new PrimaryMessenger(this);
	source = G4ThreeVector(0,0,-810) * mm;
	cosTheta = cos(22*deg);
	fPrimary->SetParticlePosition(rot * source);
	carm = _carm;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fPrimary;
	delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	if(isFirst) {
		SetSourceEnergy();
		isFirst = false;
	}

	G4double rand_energy = G4UniformRand();

	if (rand_energy == 1)	rand_energy = pdf_sort.rbegin()->second;
	else {
		for (auto itr : cdf_sort) {
			if (rand_energy < itr.first) {
				rand_energy = itr.second;
				break;
			}
		}
	}

	fPrimary->SetParticleEnergy(rand_energy*keV);
	fPrimary->SetParticleMomentumDirection(SampleRectangularBeamDirection());
	fPrimary->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SetSourceEnergy()
{
//	G4String fileName("E" + to_string(peak_energy) + "_Al" + to_string(filter_thickness).substr(0,3) + ".spec");
	peak_energy = carm->GetTubeVoltage(frameNo);
	G4String fileName("E" + to_string(peak_energy) + ".spec");
	G4String spectra("./spectra/" + fileName);
	G4cout << "Read x-ray spectra: " << spectra << G4endl;
	ifstream ifs(spectra);
	if(!ifs.is_open()) { G4cerr << "X-ray spectra file was not opened" << G4endl; exit(1); }

	G4double sum(0);
	G4String dump;
	while(getline(ifs,dump)) {
		stringstream ss(dump);
		ss >> dump;
		if (dump == "Energy[keV]") {
			while (getline(ifs,dump)) {
				G4double energy, intensity;
				stringstream ss2(dump);
				ss2 >> energy >> intensity;
				sum += energy * intensity;
				pdf[energy] = energy*intensity;
			}
		}
	}

	for (auto &itr : pdf) {
		itr.second /= sum;
		pdf_sort[itr.second] = itr.first;
	}

	G4double cdf(0);
	for (auto itr : pdf_sort) {
		cdf += itr.first;
		cdf_sort[cdf] = itr.second;
	}

	eIntensity = sum;
	monitorPower = carm->GetMonitorPower(frameNo);
	monitorTime  = carm->GetMonitorTime(frameNo);
	monitorDAP   = carm->GetMonitorDAP(frameNo);
	tubeVoltage  = carm->GetTubeVoltage(frameNo);
	tubeCurrent  = carm->GetTubeCurrent(frameNo);

	ifs.close();
}


G4ThreeVector PrimaryGeneratorAction::SampleRectangularBeamDirection()
{
	G4ThreeVector ref(0,0,0);
	G4ThreeVector centerVec = ref - source;
	G4double a = fabs(centerVec.z()) * tan(8.216*deg);
	G4double b = fabs(centerVec.z()) * tan(10.565*deg);
	G4double rand1 = G4UniformRand();
	G4double rand2 = G4UniformRand();
	G4double x = 2*a*rand1 - a;
	G4double y = 2*b*rand2 - b;
	G4double z = centerVec.z();
	G4ThreeVector rec_source(x,y,z);

	return rot * rec_source;
}
