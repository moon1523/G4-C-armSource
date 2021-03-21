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
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "G4Timer.hh"

#include <fstream>


class Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:
    RunAction(TETModelImport* _tetData, G4String _output, G4Timer* initTimer);
    virtual ~RunAction();

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void PrintPLY(G4String fileName, std::vector<G4double> vertexWeight);



  private:
    TETModelImport* tetData;
    G4String        outputFile;
    G4Timer*        initTimer;
    G4Timer*        runTimer;

    Run*            fRun;
    G4int 			numOfEvent;

    G4bool          monitorPower;
    G4double        monitorTime;
    G4double        monitorDAP;
    G4int           tubeVoltage;
    G4double        tubeCurrent;
    G4double        eIntensity;


    G4double   ratioDAP;
    map<G4int, G4double> massMap;
    map<G4int, pair<G4double,G4double>> doses;
    map<G4int, pair<G4double,G4double>> dosesRunSum;

    ofstream ofs;
    vector<G4double> vertexDose2WeightMerge;

	vector<G4ThreeVector> vertVec;
    vector<vector<G4int>> faceVec;


};

#endif

