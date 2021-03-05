/*
 * ChamberParameterisation.cc
 *
 *  Created on: May 18, 2020
 *      Author: sungho
 */

#include "ChamberParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"

ChamberParameterisation::ChamberParameterisation()
{
	WATER = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
	eventID = 0;
}

ChamberParameterisation::~ChamberParameterisation()
{
}

void ChamberParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
	G4double Xposition = 0;
	G4double Yposition = 0;
	G4double Zposition = 0;

	G4ThreeVector origin(Xposition, Yposition, Zposition);

	physVol->SetTranslation(origin);
	physVol->SetRotation(0);
}

//G4Material* ChamberParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable *parentTouch=0)
//{
//	return WATER;
//}
//
//void ChamberParameterisation::ComputeDimensions(G4Tet& tet, const G4int copyNo, const G4VPhysicalVolume* physVol) const
//{
//
////	tet.SetVertices(G4ThreeVector(0,0,10*cm) + copyNo*G4ThreeVector(10*cm, 0, 0),
////						        G4ThreeVector(0,10*cm,0)  + copyNo*G4ThreeVector(10*cm, 0, 0),
////								G4ThreeVector(5*cm,-5*cm,0) + copyNo*G4ThreeVector(10*cm, 0, 0),
////								G4ThreeVector(-5*cm,-5*cm,0) + copyNo*G4ThreeVector(10*cm, 0, 0));
//
//
//}

G4VSolid* ChamberParameterisation::ComputeSolid(const G4int copyNo, G4VPhysicalVolume* physVol)
{
	if(G4RunManager::GetRunManager()->GetCurrentEvent())
		eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	//eventID=1;
	auto myTet = dynamic_cast<G4Tet*>(physVol->GetLogicalVolume()->GetSolid());
	myTet->SetVertices(G4ThreeVector(0,0,10*cm)     + (copyNo+1)*G4ThreeVector(10*cm, 0, 0) + eventID * G4ThreeVector(10*cm, 0, 0),
	        		   G4ThreeVector(0,10*cm,0)     + (copyNo+1)*G4ThreeVector(10*cm, 0, 0) + eventID * G4ThreeVector(10*cm, 0, 0),
	        		   G4ThreeVector(5*cm,-5*cm,0)  + (copyNo+1)*G4ThreeVector(10*cm, 0, 0) + eventID * G4ThreeVector(10*cm, 0, 0),
	        		   G4ThreeVector(-5*cm,-5*cm,0) + (copyNo+1)*G4ThreeVector(10*cm, 0, 0) + eventID * G4ThreeVector(10*cm, 0, 0));


	return myTet;
}


