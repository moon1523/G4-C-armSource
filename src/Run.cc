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
// TETRun.cc
// \file   MRCP_GEANT4/External/src/TETRun.cc
// \author Haegin Han
//


#include "Run.hh"

Run::Run(TETModelImport* _tetData)
:G4Run(), tetData(_tetData), fCollID_skinTet(-1), eIntensity(0),
 monitorPower(0), monitorTime(0), monitorDAP(0),
 tubeVoltage(0), tubeCurrent(0)
{
	G4cout << "Run()" << G4endl;
	fCollID_skinTet = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/eDep");
}

Run::~Run()
{
	G4cout << "~Run()" << G4endl;
	edepMap.clear();
}

void Run::RecordEvent(const G4Event* event)
{
	//Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	G4THitsMap<G4double>* evtMap =
	static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_skinTet));

	auto hitsMap = *evtMap->GetMap();

	for(auto itr:hitsMap)
	{
		G4double edep = *(itr.second);
		G4double edepSquare = (*itr.second) * (*itr.second);
		edepMap[itr.first].first += edep;
		edepMap[itr.first].second += edepSquare;
	}

	G4Run::RecordEvent(event);
}

void Run::Merge(const G4Run* run)
{
	const Run* localRun = static_cast<const Run*>(run);

	eIntensity = localRun->eIntensity;

	monitorPower = localRun->monitorPower;
	monitorTime  = localRun->monitorTime;
	monitorDAP   = localRun->monitorDAP;
	tubeVoltage  = localRun->tubeVoltage;
	tubeCurrent  = localRun->tubeCurrent;

	auto localMap = localRun->edepMap;
	for (auto itr:localMap) {
		edepMap[itr.first].first  += itr.second.first;
		edepMap[itr.first].second += itr.second.second;
	}
	G4Run::Merge(run);
}

