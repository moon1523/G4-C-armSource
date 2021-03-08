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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"


RunAction::RunAction(TETModelImport* _tetData)
: G4UserRunAction(), tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile("out.txt")
, primaryEnergy(-1.), beamArea(-1.)
{
	if (!isMaster) return;
}


RunAction::~RunAction()
{

}


G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new Run(tetData);
	return fRun;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{


}


void RunAction::EndOfRunAction(const G4Run* aRun)
{
	G4int nofEvents = aRun->GetNumberOfEvent();
	if ( nofEvents == 0 ) return;

	const DetectorConstruction* detectorConstruction
	= static_cast<const DetectorConstruction*>
	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());


	const Run* runE = static_cast<const Run*>(aRun);
	std::vector<G4VPhysicalVolume*> phyVec = detectorConstruction->GetScoringPV();
		for (G4VPhysicalVolume* phy:phyVec) {
		G4double mass = phy->GetLogicalVolume()->GetMass();
		G4int idx = phy->GetCopyNo();
		massMap[idx] += mass;
	}

	std::map<G4int, G4double> dose;

	auto edep = runE->GetfEdep();
	for (auto itr:edep) {
		dose[itr.first] = itr.second/massMap[itr.first]/nofEvents;
	}

	if (isMaster) {
		std::ofstream ofs("result.txt");
		for (auto itr:dose) {
			ofs << itr.first << "\t" << dose[itr.first] / (joule/kg) << G4endl;
		}
	}


}
