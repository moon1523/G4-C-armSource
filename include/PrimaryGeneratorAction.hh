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

#ifndef PrimaryGeneratorAction_hh
#define PrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryMessenger.hh"
#include "G4RotationMatrix.hh"
#include "G4RandomDirection.hh"

#include "CarmTracking.hh"

#include <map>
#include <algorithm>

using namespace std;

class PrimaryMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(CarmTracking* carm);
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    void SetAngle(G4double angle)             { cosTheta = cos(angle); }
    void SetDefaultSource(G4ThreeVector pos)  { source = pos; }

    void SetSource() {
    	rot = carm->GetRotationMatrix(frameNo);
    	G4ThreeVector trans = carm->GetTranslationMatrix(frameNo);
    	fPrimary->SetParticlePosition(rot*source + trans);
    }

//    void SetSource(G4RotationMatrix _rot, G4ThreeVector trans) {
//        rot = _rot;
//        fPrimary->SetParticlePosition(rot*source + trans);
//    }

    G4ThreeVector SampleADirection(){
        return rot*G4RandomDirection(cosTheta);
    }
    G4ThreeVector SampleRectangularBeamDirection();

    void  SetPeakEnergy(G4int _peak_energy) { peak_energy      = _peak_energy;      }
    void  SetFilterThickness(G4double _filter_thickness) { filter_thickness = _filter_thickness; }
    void  SetFrameNo(G4int _frameNo) { frameNo = _frameNo; }
    
    void          SetSourceEnergy();
    G4ThreeVector SetSourcePosition();

    G4ParticleGun* GetParticleGun()     const { return fPrimary; }
    G4double       GetEIntensity()      const { return eIntensity; }

    G4bool         GetMonitorPower()    const { return monitorPower; }
    G4double       GetMonitorTime()     const { return monitorTime;  }
    G4double       GetMonitorDAP()      const { return monitorDAP;   }
    G4int          GetTubeVoltage()     const { return tubeVoltage;  }
    G4double       GetTubeCurrent()     const { return tubeCurrent;  }

  private:
    G4ParticleGun*    fPrimary;
    PrimaryMessenger* fMessenger;
    G4double          cosTheta;
    G4RotationMatrix  rot;
    G4ThreeVector     source;
    G4int             frameNo;

    map<G4double, G4double> pdf;
    map<G4double, G4double, greater<G4double>> pdf_sort;
    map<G4double, G4double> cdf_sort;
    G4double eIntensity;

    G4int peak_energy;
	G4double filter_thickness;

	G4bool   monitorPower;
	G4double monitorTime;
	G4double monitorDAP;
	G4int    tubeVoltage;
	G4double tubeCurrent;

    CarmTracking* carm;
	G4bool isFirst;
};

#endif
