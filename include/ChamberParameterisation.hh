/*
 * ChamberParameterisation.hh
 *
 *  Created on: May 18, 2020
 *      Author: sungho
 */

#ifndef INCLUDE_CHAMBERPARAMETERISATION_HH_
#define INCLUDE_CHAMBERPARAMETERISATION_HH_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4Tet.hh"
#include "G4NistManager.hh"

class ChamberParameterisation : public G4VPVParameterisation
{
public:
	ChamberParameterisation();
	virtual ~ChamberParameterisation();

	virtual void ComputeTransformation // position, rotation
		(const G4int copyNo, G4VPhysicalVolume* physVol) const;

//	virtual G4Material* ComputeMaterial // mat, sensitivity, visAtt
//		(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable *parentTouch=0);
//
//	virtual void ComputeDimensions //size
//		(G4Tet& tet, const G4int copyNo, const G4VPhysicalVolume* physVol) const;

	virtual G4VSolid* ComputeSolid // shape
		(const G4int copyNo, G4VPhysicalVolume* physVol);

private:
	G4Material* WATER;
	G4int eventID;
};

#endif /* INCLUDE_CHAMBERPARAMETERISATION_HH_ */
