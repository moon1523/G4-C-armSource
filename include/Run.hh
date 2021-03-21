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
// TETRun.hh
// \file   MRCP_GEANT4/External/include/TETRun.hh
// \author Haegin Han
//

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"
#include "TETModelImport.hh"
#include "PrimaryGeneratorAction.hh"


class Run : public G4Run
{
public:
	Run(TETModelImport* _tetData);
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
	virtual void Merge(const G4Run*);

	map<G4int, pair<G4double, G4double>> GetEdepMap() const { return edepMap; }

	G4double GetEIntensty()    { return eIntensity;   }
	G4bool   GetMonitorPower() { return monitorPower; }
	G4double GetMonitorTime()  { return monitorTime;  }
	G4double GetMonitorDAP()   { return monitorDAP;   }
	G4int    GetTubeVoltage()  { return tubeVoltage;  }
	G4double GetTubeCurrent()  { return tubeCurrent;  }


	void SetPrimary(G4double _eIntensity,
			G4bool _monitorPower, G4double _monitorTime, G4double _monitorDAP,
			G4int _tubeVoltage, G4double _tubeCurrent)
	{
		eIntensity = _eIntensity;

		monitorPower = _monitorPower;
		monitorTime  = _monitorTime;
		monitorDAP   = _monitorDAP;
		tubeVoltage  = _tubeVoltage;
		tubeCurrent  = _tubeCurrent;
	}

private:
	TETModelImport* tetData;
	G4int 			fCollID_skinTet;
	map<G4int, pair<G4double, G4double>> edepMap;

	G4double eIntensity;

	G4bool   monitorPower;
	G4double monitorTime;
	G4double monitorDAP;
	G4int    tubeVoltage;
	G4double tubeCurrent;
};
#endif

    
