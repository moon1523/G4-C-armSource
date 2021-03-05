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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "../include/DetectorConstruction.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "ChamberParameterisation.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
	G4Material* AIR = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	G4Material* WATER = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
	// World
	G4VSolid* sol_World = new G4Box("World", 1.0*m, 1.0*m, 1.0*m);
	G4LogicalVolume* lv_World = new G4LogicalVolume(sol_World, AIR, "World");
	G4VPhysicalVolume* pv_World =
	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), lv_World, "World", 0, false, 0);

	G4VSolid* sol_Tet = new G4Tet("ee", G4ThreeVector(0,0,10*cm), G4ThreeVector(0,10*cm,0), G4ThreeVector(5*cm,-5*cm,0),G4ThreeVector(-5*cm,-5*cm,0));
	//G4VSolid* sol_Tet1 = new G4Tet("ee", G4ThreeVector(0,0,30*cm), G4ThreeVector(0,10*cm,20*cm), G4ThreeVector(5*cm,-5*cm,20*cm),G4ThreeVector(-5*cm,-5*cm,20*cm));
	G4LogicalVolume* lv_Tet = new G4LogicalVolume(sol_Tet, WATER, "Tet");
	//G4LogicalVolume* lv_Tet1 = new G4LogicalVolume(sol_Tet1, WATER, "Tet1");

	//	G4PVPlacement* phy_Tet = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), lv_Tet, "Tet", lv_World, false, 1);
	//G4VPVParameterisation* chamberParam = new ChamberParameterisation();
	//G4PVParameterised* pv_Tet = new G4PVParameterised("para",lv_Tet, lv_World, kXAxis, 3, chamberParam);
	//G4PVParameterised* pv_Tet1 = new G4PVParameterised("para1",lv_Tet1, lv_World, kXAxis, 3, chamberParam); //dksehlsek



	// Visualization
	G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	va_World->SetForceWireframe(true);
	lv_World->SetVisAttributes(va_World);

	G4VisAttributes* va_Tet = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1));
	va_Tet->SetForceSolid(true);
	lv_Tet->SetVisAttributes(va_Tet);


	return pv_World;
}

void DetectorConstruction::ConstructSDandField()
{
	G4MultiFunctionalDetector* fMFD = new G4MultiFunctionalDetector("detector");
	fMFD->RegisterPrimitive(new G4PSEnergyDeposit("TetDose"));
	G4SDManager::GetSDMpointer()->AddNewDetector(fMFD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
