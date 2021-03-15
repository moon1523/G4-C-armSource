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

#include "DetectorConstruction.hh"

using namespace std;

DetectorConstruction::DetectorConstruction(TETModelImport* _tetData, CarmTracking* _carm)
:worldLogical(0) ,worldPhysical(0), container_logic(0), tetData(_tetData),
 carm(_carm), tetLogic(0), frameNo(0),
 lv_DAP(0), pv_DAP(0),
 lv_pic(0), pv_pic(0),
 a(30*cm * tan(11*deg))
{
	G4cout << "DetectorConstruction() !!!" << G4endl;
	// initialisation of the variables for phantom information
	phantomSize     = tetData -> GetPhantomSize();
	phantomBoxMin   = tetData -> GetPhantomBoxMin();
	phantomBoxMax   = tetData -> GetPhantomBoxMax();
	nOfTetrahedrons = tetData -> GetNumTetrahedron();

	fMessenger = new PrimaryMessenger(this);
	rot_DAP = new G4RotationMatrix();
}

DetectorConstruction::~DetectorConstruction()
{
	delete tetData;
	delete fMessenger;
	delete rot_DAP;
	delete carm;
	delete pv_DAP;
	G4cout << "~DetectorConstruction() !!!" << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
	SetupWorldGeometry();
	ConstructPhantom();
	PrintPhantomInformation();

	return worldPhysical;
}

void DetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 10*10*10 m3)
	//
	G4double worldXYZ = 10. * m;
	G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

	G4VSolid* worldSolid
	  = new G4Box("worldSolid", worldXYZ/2, worldXYZ/2, worldXYZ/2);

	worldLogical
	  = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");

	worldPhysical
	  = new G4PVPlacement(0,G4ThreeVector(), worldLogical,"worldPhysical", 0, false,0,false);


	G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	va_World->SetForceWireframe(true);
	worldLogical->SetVisAttributes(va_World);

	// Define the phantom container (10-cm margins from the bounding box of phantom)
	//
	G4Box* containerSolid = new G4Box("phantomBox", phantomSize.x()/2 + 1.*cm,
										            phantomSize.y()/2 + 1.*cm,
										            phantomSize.z()/2 + 1.*cm);

	container_logic = new G4LogicalVolume(containerSolid, vacuum, "phantomLogical");
	G4VisAttributes* container_vis = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	container_vis->SetForceWireframe(true);
	container_logic->SetVisAttributes(container_vis);


	G4RotationMatrix* rotM = new G4RotationMatrix();
	rotM->rotateY(192.6*deg);
	new G4PVPlacement(rotM, G4ThreeVector(374.91301*mm,-77.09209*mm,2480.90340*mm), container_logic, "PhantomPhysical",
			          worldLogical, false, 0);
	container_logic->SetOptimisation(TRUE);
	container_logic->SetSmartless( 0.5 ); // for optimization (default=2)

	G4Box* sol_DAP = new G4Box("DAP_meter", a, a, 0.5*1/a);
	lv_DAP = new G4LogicalVolume(sol_DAP, vacuum, "lv_DAP");
	lv_DAP->SetVisAttributes(G4Colour(1.,1.,0.));
	pv_DAP = new G4PVPlacement(0, G4ThreeVector(0,0,-810*mm), lv_DAP, "pv_DAP", worldLogical, false, 1000);
//	scoringPV.push_back(pv_DAP);
}

void DetectorConstruction::ConstructPhantom()
{
//	// View
//	G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
//	lv_pic = new G4LogicalVolume(tetData->GetPicTess(), water, "lv_pic");
//	lv_pic->SetVisAttributes(new G4VisAttributes(G4Colour(1.000000, 0.752941, 0.627451)));
//	new G4PVPlacement(0, G4ThreeVector(), lv_pic, "pv_pic", container_logic, false, 999);

	// Skin 50-100 um ( 3 Tet in 1 Facet Layer )
	vector<vector<G4int>> innerFaces = tetData->GetInnerFaces();
	vector<vector<G4int>> outerFaces = tetData->GetOuterFaces();
	G4int skinCount(0);
	G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
	for (int i=0;i<innerFaces.size();i++) {
		G4VSolid* skin_Tet1 = new G4Tet("SkinTet1_Solid",
										tetData->GetAVertex(outerFaces[skinCount][0])-center,
										tetData->GetAVertex(outerFaces[skinCount][1])-center,
										tetData->GetAVertex(outerFaces[skinCount][2])-center,
										tetData->GetAVertex(innerFaces[skinCount][1])-center);
		G4VSolid* skin_Tet2 = new G4Tet("SkinTet2_Solid",
										tetData->GetAVertex(outerFaces[skinCount][0])-center,
										tetData->GetAVertex(innerFaces[skinCount][1])-center,
										tetData->GetAVertex(outerFaces[skinCount][2])-center,
										tetData->GetAVertex(innerFaces[skinCount][0])-center);
		G4VSolid* skin_Tet3 = new G4Tet("SkinTet3_Solid",
										tetData->GetAVertex(innerFaces[skinCount][0])-center,
										tetData->GetAVertex(innerFaces[skinCount][2])-center,
										tetData->GetAVertex(outerFaces[skinCount][2])-center,
										tetData->GetAVertex(innerFaces[skinCount][1])-center);
		G4LogicalVolume* lv_Tet1 = new G4LogicalVolume(skin_Tet1, tetData->GetMaterial(126), "lv_Tet1");
		G4LogicalVolume* lv_Tet2 = new G4LogicalVolume(skin_Tet2, tetData->GetMaterial(126), "lv_Tet2");
		G4LogicalVolume* lv_Tet3 = new G4LogicalVolume(skin_Tet3, tetData->GetMaterial(126), "lv_Tet3");

		G4PVPlacement* pv_Tet1 = new G4PVPlacement(0, G4ThreeVector(), lv_Tet1, "pv_Tet1", container_logic, false, 12600+skinCount);
		G4PVPlacement* pv_Tet2 = new G4PVPlacement(0, G4ThreeVector(), lv_Tet2, "pv_Tet2", container_logic, false, 12600+skinCount);
		G4PVPlacement* pv_Tet3 = new G4PVPlacement(0, G4ThreeVector(), lv_Tet3, "pv_Tet3", container_logic, false, 12600+skinCount);
		scoringPV.push_back(pv_Tet1);
		scoringPV.push_back(pv_Tet2);
		scoringPV.push_back(pv_Tet3);
		skinCount++;
	}

	// Phantom
	for (G4int i=0; i<tetData->GetNumTetrahedron(); i++) {
		G4Material* mat_Tet = tetData->GetMaterial(tetData->GetMaterialIndex(i));
		G4VSolid* sol_Tet = tetData->GetTetrahedron(i);
		G4LogicalVolume* lv_Tet = new G4LogicalVolume(sol_Tet, mat_Tet, "lv_Tet");

		if (tetData->GetMaterialIndex(i) == 126) {
			continue;
		} else {
			new G4PVPlacement(0, G4ThreeVector(), lv_Tet, "pv_Tet",
													  container_logic, false, tetData->GetMaterialIndex(i));
		}
	}
}

void DetectorConstruction::ConstructSDandField()
{
	G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector("PhantomSD");
	MFDet->RegisterPrimitive(new G4PSEnergyDeposit("eDep"));
	G4SDManager::GetSDMpointer()->AddNewDetector(MFDet);

	for (auto itr:scoringPV)
		SetSensitiveDetector(itr->GetLogicalVolume(), MFDet);
	// DAP meter
	SetSensitiveDetector(lv_DAP, MFDet);

	G4cout << "SCORING SIZE:" << scoringPV.size() << G4endl;
}

void DetectorConstruction::PrintPhantomInformation()
{
	// print brief information on the imported phantom
	G4cout<< G4endl;
	G4cout.precision(3);
	G4cout<<"   Phantom name               "<<tetData->GetPhantomName() << " TET phantom"<<G4endl;
	G4cout<<"   Phantom size               "<<phantomSize.x()<<" * "<<phantomSize.y()<<" * "<<phantomSize.z()<<" mm3"<<G4endl;
	G4cout<<"   Phantom box position (min) "<<phantomBoxMin.x()<<" mm, "<<phantomBoxMin.y()<<" mm, "<<phantomBoxMin.z()<<" mm"<<G4endl;
	G4cout<<"   Phantom box position (max) "<<phantomBoxMax.x()<<" mm, "<<phantomBoxMax.y()<<" mm, "<<phantomBoxMax.z()<<" mm"<<G4endl;
	G4cout<<"   Number of tetrahedrons     "<<nOfTetrahedrons<<G4endl<<G4endl;
}
