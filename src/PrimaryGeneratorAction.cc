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

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "../include/PrimaryGeneratorAction.hh"
#include "G4RandomDirection.hh"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <tuple>
#include "G4RotationMatrix.hh"
#include <iomanip>
using namespace std;


PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(), peak_energy("100"), filter_thickness("0.2"),
  	  	  	  	  	  	  	  	   l_arm_pos(0),l_arm_rot(0),
								   c_arm_rot(0),c_arm_ang(0), isfirst(true)
{
	fMessenger = new PrimaryMessenger(this);

	sourcePosition = G4ThreeVector(0.0*mm, -810*mm, 0.0*mm);
//	SetSourceEnergy();
//	SetSourcePosition();

	fPrimary = new G4ParticleGun();
	fPrimary->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("geantino"));
//	fPrimary->SetParticlePosition(sourcePosition);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fPrimary;
	delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	if(isfirst)
	{
		SetSourceEnergy();
		SetSourcePosition();
		sourcePosition += G4ThreeVector(0, 0, l_arm_pos);
		fPrimary->SetParticlePosition(sourcePosition);

		isfirst = false;
	}

	double rand_energy = G4UniformRand();

	if (rand_energy == 1)
		rand_energy = pdf_sort.rbegin()->second;
	else
	{
		for (auto itr : cdf_sort)
		{
			if (rand_energy < itr.first)
			{
				rand_energy = itr.second;
				break;
			}
		}
	}

	fPrimary->SetParticleEnergy(rand_energy*keV);
	G4ThreeVector u = G4ThreeVector(1,0,0);
	G4ThreeVector v = G4ThreeVector(0, cos(M_PI*0.5), -sin(M_PI*0.5));
	G4ThreeVector w = G4ThreeVector(0, sin(M_PI*0.5),  cos(M_PI*0.5));
	G4RotationMatrix rot = G4RotationMatrix(u,v,w);
	fPrimary->SetParticleMomentumDirection(sourceRotM*rot*G4RandomDirection(cos(22*deg))); // 22 deg (anode angle = 11 deg)

	//G4cout << "sampling_energy=" << rand_energy << G4endl;
	fPrimary->GeneratePrimaryVertex(anEvent);

}

void PrimaryGeneratorAction::SetSourceEnergy()
{
//	cout << peak_energy << endl;
//	cout << filter_thickness << endl;
//	stringstream stream;
//	stream << fixed << setprecision(1) << filter_thickness;
//	G4String ss = stream.str();

	G4String fileName;
	fileName = "E" + peak_energy + "_Al" + filter_thickness + ".spec";
	G4String spectra_folder = "./spectra/";

	spectra_folder += fileName;
	G4cout << spectra_folder << G4endl;
	std::ifstream ifs(spectra_folder);
	if(!ifs.is_open())
		G4cout << "X-ray spectra file was not opened" << G4endl;

	double energy(0);
	double intensity(0);

	double sum(0);
	while(!ifs.eof())
	{
		ifs >> energy >> intensity;
		sum += energy * intensity;
		pdf.insert(make_pair(energy, energy*intensity));
	}

	for (auto &itr : pdf)
	{
		itr.second /= sum;
		pdf_sort[itr.second] = itr.first;
	}

	double cdf(0);
	for (auto itr : pdf_sort)
	{
		cdf += itr.first;
		cdf_sort.insert(make_pair(cdf, itr.second));
	}

	ifs.close();
}

G4ThreeVector PrimaryGeneratorAction::SetSourcePosition()
{
	// +:counterclockwise, -:clockwise (criteria : from the top-view of the reference axis)
	// R_y => +:counterclockwise(~90), -:clockwise(~90)
	// R_z => +:LAO(~120/~90), -:RAO(~185/~90)
	// R_x => +:caudal(~90/~120), -:cranial(~90/~185)
	sourceRotM.rotateX(c_arm_ang).rotateZ(c_arm_rot).rotateY(l_arm_rot);
//	sourcePosition += G4ThreeVector(0, 0, l_arm_pos);

	return sourcePosition = sourceRotM * sourcePosition; // R = R_x * R_z * R_y => XZY Euler Angle
}
