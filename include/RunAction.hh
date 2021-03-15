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

#include <fstream>


class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:
    RunAction(TETModelImport* _tetData);
    virtual ~RunAction();

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void PrintPLY(G4String fileName, std::vector<G4double> vertexWeight);
    void PrintSUMPLY(G4String fileName);


  private:
    TETModelImport* tetData;
    Run*            fRun;
    G4int           numOfEvent;
    G4int           runID;
    G4String        outputFile;

    G4String primaryParticle;
    G4String primarySourceName;
    G4double primaryEnergy;
    G4double beamArea;


    std::map<G4int, G4double> massMap;
    std::map<G4int, std::pair<G4int, G4double>> massMap2;
    std::map<G4int, std::pair<G4double,G4double>> doses;
    std::map<G4int, G4String> nameMap;

    std::ofstream ofs;
    std::vector<G4double> vertexDose2WeightMerge;

	std::vector<G4ThreeVector> vertVec;
    std::vector<std::vector<G4int>> faceVec;
};

#endif

