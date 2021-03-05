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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"
#include <fstream>
#include <map>
#include "PrimaryMessenger.hh"

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

//    void SetSourceEnergy(G4int _peak_energy, G4double _filter_thickness);
//    G4ThreeVector SetSourcePosition(G4double _l_arm_pos, G4double _l_arm_rot,
//    					   	   	    G4double _c_arm_rot,
//								   	G4double _c_arm_ang);
    void     GetPeakEnergy(G4String _peak_energy)           { peak_energy      = _peak_energy;      }
    void     GetFilterThickness(G4String _filter_thickness) { filter_thickness = _filter_thickness; }
    void     GetLarmPosition(G4double _l_arm_pos)           { l_arm_pos        = _l_arm_pos;        }
    void     GetLarmRotation(G4double _l_arm_rot)           { l_arm_rot        = _l_arm_rot;        }
    void     GetCarmRotation(G4double _c_arm_rot)           { c_arm_rot        = _c_arm_rot;        }
    void     GetCarmAngulation(G4double _c_arm_ang)         { c_arm_ang        = _c_arm_ang;        }

    void SetSourceEnergy();
    G4ThreeVector SetSourcePosition();

  private:
    //G4GeneralParticleSource* fPrimary;
    G4ParticleGun* fPrimary;
    PrimaryMessenger* fMessenger;

    std::map<G4double, G4double> pdf;
	std::map<G4double, G4double, std::greater<G4double>> pdf_sort;
	std::map<G4double, G4double> cdf_sort;

	G4String peak_energy;
	G4String filter_thickness;
	G4double l_arm_pos;
	G4double l_arm_rot;
	G4double c_arm_rot;
	G4double c_arm_ang;

	G4ThreeVector sourcePosition;
	G4RotationMatrix sourceRotM;
	G4bool isfirst;
};

#endif
