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
// TETPrimaryMessenger.cc
// \author Haegin Han
//

#include "PrimaryMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ParticleGun.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
	fSourceDef = new G4UIdirectory("/source/");

	fsrc_energy    = new G4UIcmdWithAString("/source/energy", this);
	fsrc_filter    = new G4UIcmdWithAString("/source/filter", this);
	fsrc_l_arm_pos = new G4UIcmdWithADoubleAndUnit("/source/lpos", this); fsrc_l_arm_pos->SetDefaultUnit("mm");
	fsrc_l_arm_rot = new G4UIcmdWithADoubleAndUnit("/source/lrot", this); fsrc_l_arm_rot->SetDefaultUnit("deg");
	fsrc_c_arm_rot = new G4UIcmdWithADoubleAndUnit("/source/crot", this); fsrc_c_arm_rot->SetDefaultUnit("deg");
	fsrc_c_arm_ang = new G4UIcmdWithADoubleAndUnit("/source/cang", this); fsrc_c_arm_ang->SetDefaultUnit("deg");

//	fExternalDir = new G4UIdirectory("/external/");
//	fBeamDirCmd = new G4UIcmdWithAString("/external/dir", this);
//	fBeamDirCmd->SetCandidates("AP PA LLAT RLAT ROT ISO");
//
//	fInternalDir      = new G4UIdirectory("/internal/");
//	fSourceOrganCmd   = new G4UIcmdWithAString("/internal/source", this);
//	fSurfaceSourceCmd = new G4UIcmdWithAString("/internal/surface", this);
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fSourceDef;
	delete fsrc_energy;
	delete fsrc_filter;
	delete fsrc_l_arm_pos;
	delete fsrc_l_arm_rot;
	delete fsrc_c_arm_rot;
	delete fsrc_c_arm_ang;

//	delete fExternalDir;
//	delete fBeamDirCmd;
//	delete fInternalDir;
//	delete fSourceOrganCmd;
}

//void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4int _energy, G4double _thickness,
//								   G4double _l_pos, G4double _l_rot, G4double _c_rot, G4double _c_ang)
void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	     if(command == fsrc_energy)		fPrimary->GetPeakEnergy(newValue);
	else if(command == fsrc_filter)     fPrimary->GetFilterThickness(newValue);
	else if(command == fsrc_l_arm_pos)  fPrimary->GetLarmPosition(fsrc_l_arm_pos->GetNewDoubleValue(newValue));
	else if(command == fsrc_l_arm_rot)  fPrimary->GetLarmRotation(fsrc_l_arm_rot->GetNewDoubleValue(newValue));
	else if(command == fsrc_c_arm_rot)  fPrimary->GetCarmRotation(fsrc_c_arm_rot->GetNewDoubleValue(newValue));
	else if(command == fsrc_c_arm_ang)  fPrimary->GetCarmAngulation(fsrc_c_arm_ang->GetNewDoubleValue(newValue));

}

