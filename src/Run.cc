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
:G4Run(), tetData(_tetData), fCollID_skinTet(-1)
{
	fCollID_skinTet = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/eDep");
}

Run::~Run()
{

}

void Run::RecordEvent(const G4Event* event)
{
	//Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	//Energy in crystals : identify 'good events'
	//
	G4THitsMap<G4double>* evtMap =
	static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_skinTet));

	auto hitsMap = *evtMap->GetMap();
	for(auto itr:hitsMap)
	{

		G4double edep = *(itr.second);
		fEdep[itr.first] += edep;
	}

	std::ofstream ofs("hitsMap.txt",std::ios::app);
	for(auto itr:hitsMap) {
		G4double edep = *(itr.second);
		fEdep[itr.first] += edep;
		ofs << itr.first << " " << *(itr.second) << " " << edep << G4endl;
	}

	G4Run::RecordEvent(event);
}

void Run::Merge(const G4Run* run)
{
	const Run* localRun = static_cast<const Run*>(run);

	auto localEdep = localRun->fEdep;
	for(auto itr:localEdep)
		fEdep[itr.first] += itr.second;

	G4Run::Merge(run);
}

