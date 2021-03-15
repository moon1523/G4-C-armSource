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
//  Author: Sungho Moon (2021.03.05)
//  
/// \file MRCGS0407.cc
/// \brief main source



#include "G4UImanager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "ActionInitialization.hh"

#include "TETModelImport.hh"
#include "CarmTracking.hh"

#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

void PrintUsage(){
  G4cerr<< "Usage: ./TetCal -m [MACRO] -o [OUTPUT] -p [phantom name]"  <<G4endl;
  G4cerr<< "Example: ./TetCal -m sample.in -o run.out -p ./phantoms/00M" <<G4endl;
}

int main(int argc,char** argv)
{
  // Read the arguments for batch mode
	//
    G4Timer* initTimer = new G4Timer;
    initTimer->Start();
    G4String macro;
    G4String output;
    G4String phantomName;
    G4UIExecutive* ui = 0;


  for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[i+1];
			i++;
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[i+1];
			i++;
		}
		// switch for MRCP-AF phantom
		else if ( G4String(argv[i]) == "-p" ) {
			phantomName = argv[i+1];
			i++;
		}
		else {
			PrintUsage();
			return 1;
		}
	}

  // print usage when there are more than six arguments
	if ( argc>7 || phantomName.empty()){
		PrintUsage();
		return 1;
	}
  
  // Detect interactive mode (if no macro file name) and define UI session
	//
	if ( !macro.size() ) {
		ui = new G4UIExecutive(argc, argv, "qt");
//		G4cerr<<"ERROR: Interactive mode is not available. Please provide macro file."<<G4endl;
//		return 1;
	}

  // Choose the Random engine
	//
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set a class to import phantom data, tracking result
  //
  TETModelImport* tetData = new TETModelImport(phantomName, ui);
  CarmTracking* carm = new CarmTracking();

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction(tetData, carm));
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization(tetData, carm));
  
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}
