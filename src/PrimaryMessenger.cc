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

#include "PrimaryMessenger.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary), fDetector(0), fDetectorDef(0), fdet_frame(0)
{
	G4cout << "PrimaryMessenger()" << G4endl;
	fSourceDef = new G4UIdirectory("/source/");

	fsrc_frame     = new G4UIcmdWithAnInteger("/source/frame", this);
	fsrc_energy    = new G4UIcmdWithAnInteger("/source/energy", this);
	fsrc_filter    = new G4UIcmdWithADouble("/source/filter", this);
	
}

PrimaryMessenger::PrimaryMessenger(DetectorConstruction* _detector)
:G4UImessenger(), fPrimary(0), fDetector(_detector), fSourceDef(0), fsrc_frame(0), fsrc_energy(0), fsrc_filter(0)
{
	G4cout << "~PrimaryMessenger()" << G4endl;
	fDetectorDef = new G4UIdirectory("/detector/");
	fdet_frame = new G4UIcmdWithAnInteger("/detector/frame", this);
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fSourceDef;
	delete fsrc_frame;
	delete fsrc_energy;
	delete fsrc_filter;

	delete fDetectorDef;
	delete fdet_frame;
}

void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
		 if(command == fsrc_frame )	 	{ fPrimary->SetFrameNo(fsrc_frame->GetNewIntValue(newValue));
		 	 	 	 	 	 	 	 	  fPrimary->SetSource(); }
	else if(command == fsrc_energy)	    { fPrimary->SetPeakEnergy(fsrc_energy->GetNewIntValue(newValue));
		 	 	 	 	 	 	 	 	  fPrimary->SetSourceEnergy(); }
	else if(command == fsrc_filter)       fPrimary->SetFilterThickness(fsrc_filter->GetNewDoubleValue(newValue));
	else if(command == fdet_frame)      { fDetector->SetFrameNo(fsrc_frame->GetNewIntValue(newValue));
			 	 	 	 	 	 	 	  fDetector->SetDAPMeter(); }
}
