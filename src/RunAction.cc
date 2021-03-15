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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"


RunAction::RunAction(TETModelImport* _tetData)
: G4UserRunAction(), tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile("out.txt")
, primaryEnergy(-1.), beamArea(-1.)
{
	if (!isMaster) return;
	G4cout << "RunAction()" << G4endl;
	ofs.open("result.txt");
	for(G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge.push_back(0);
	}

}


RunAction::~RunAction()
{
	if(!isMaster) return;
	G4cout << "~RunAction() !!!" << G4endl;

	ofstream vert("vertexDoseMerge.txt");
	for (auto itr:vertexDose2WeightMerge) {
		vert << itr / (joule/kg) << G4endl;
	}

	G4double sumvertDose(0);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		sumvertDose += vertexDose2WeightMerge[i];
	}
	vert << "sumvertDose" << G4endl;
	vert << sumvertDose << G4endl;
	vert << "vertexWeightMerge" << G4endl;
	G4double weightMin(DBL_MAX);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] /= sumvertDose;
		if ((vertexDose2WeightMerge[i] > 0) && (vertexDose2WeightMerge[i] < weightMin))
			weightMin = vertexDose2WeightMerge[i];
		vert << vertexDose2WeightMerge[i] << G4endl;
	}

	vert << "weightMin" << G4endl;
	vert << weightMin << G4endl;

	G4double logMax(-DBL_MAX);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] += weightMin;
		vertexDose2WeightMerge[i] /= weightMin;
		vertexDose2WeightMerge[i] = log10(vertexDose2WeightMerge[i]);
		if (logMax < vertexDose2WeightMerge[i]) logMax = vertexDose2WeightMerge[i];
	}
	ofstream fout("sum.txt");
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] /= logMax;
		fout << vertexDose2WeightMerge[i] << G4endl;
	}

	PrintPLY("sum.ply", vertexDose2WeightMerge);

	ofs.close();
}


G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new Run(tetData);
	return fRun;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{


}


void RunAction::EndOfRunAction(const G4Run* aRun)
{
	if(!isMaster) return;

	G4int nofEvents = aRun->GetNumberOfEvent();
	if ( nofEvents == 0 ) return;

	const DetectorConstruction* detectorConstruction
	= static_cast<const DetectorConstruction*>
	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());


	const Run* runE = static_cast<const Run*>(aRun);
	std::vector<G4VPhysicalVolume*> phyVec = detectorConstruction->GetScoringPV();

	for (G4VPhysicalVolume* phy:phyVec) {
			G4double mass = phy->GetLogicalVolume()->GetMass();
			G4int idx = phy->GetCopyNo();
			massMap[idx] += mass;
	}

	std::map<G4int, G4double> dose;

	auto edep = runE->GetfEdep();
	for (auto itr:edep) {
		dose[itr.first] = itr.second/massMap[itr.first]/nofEvents;
	}

	G4double areaDAP = detectorConstruction->GetDAParea();
	ofs << "Frame " << aRun->GetRunID() << G4endl;
	ofs << "SkinFaceID \t Dose(J/kg)" << G4endl;
	for (auto itr:dose) {
		if (itr.first == 1000) { ofs << "DAP (Gy): " << itr.second * areaDAP / (joule/kg) / (cm*cm) << G4endl; continue; }
		ofs << itr.first-12600 << "\t" << itr.second / (joule/kg) << G4endl;
	}
	ofs << G4endl;

	vector<G4double> faceDose;
	vector<G4double> vertexDose2Weight;
	vertVec = tetData->GetOuterVec();
	faceVec = tetData->GetOuterOgFaces();

	for (G4int i=0; i<tetData->GetSkinFaceSize(); i++) {
		G4bool pass(false);
		for (auto itr:dose) {
			if (itr.first == 1000) continue;
			if ((itr.first-12600) == i) { faceDose.push_back(itr.second); pass = true; break;}
		}
		if (pass) continue;
		faceDose.push_back(0);
	}

	G4double sumVert(0);
	for (G4int i=0; i<tetData->GetSkinVertexSize(); i++) {
		G4double sumAdjacent(0);
		for (G4int j=0; j<faceVec.size(); j++) {
			if (i == faceVec[j][0] || i == faceVec[j][1] || i == faceVec[j][2]) {
				sumAdjacent += faceDose[j];
			}
		}
		sumVert += sumAdjacent;
		vertexDose2Weight.push_back(sumAdjacent);
	}

	G4double weightMin(DBL_MAX);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] += vertexDose2Weight[i];
		vertexDose2Weight[i] /= sumVert;
		if ((vertexDose2Weight[i] > 0) && (vertexDose2Weight[i] < weightMin)) {
			weightMin = vertexDose2Weight[i];
		}

	}

	G4double logMax(-DBL_MAX);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2Weight[i] += weightMin;
		vertexDose2Weight[i] /= weightMin;
		vertexDose2Weight[i] = log10(vertexDose2Weight[i]);
		if (logMax < vertexDose2Weight[i]) logMax = vertexDose2Weight[i];
	}

	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2Weight[i] /= logMax;
	}

	G4String fileName = "frame" + to_string(aRun->GetRunID()) + ".ply";
	PrintPLY(fileName, vertexDose2Weight);
}

void RunAction::PrintPLY(G4String fileName, vector<G4double> vertexWeight)
{
	vector<tuple<G4int,G4int,G4int>> tetColor;

	for (auto v:vertexWeight) {
		G4int rr, gg, bb;
		if (v < 0.5) {
			rr = 0;
			gg = floor(255 * v * 2 + 0.5);
			bb = floor(255 * (0.5 - v) * 2 + 0.5);
		} else {
			rr = floor(255 * v * 2 + 0.5);
			gg = floor(255 * (1-v) * 2 + 0.5);
			bb = 0;
		}
		tetColor.push_back(make_tuple(rr,gg,bb));
	}

	ofstream ofs_ply("./SkinDistVertex/" + fileName);
	ofs_ply << "ply" << G4endl;
	ofs_ply << "format ascii 1.0" << G4endl;
	ofs_ply << "comment exported by rapidform" << G4endl;
	ofs_ply << "element vertex " << vertVec.size() << G4endl;
	ofs_ply << "property float x" << G4endl;
	ofs_ply << "property float y" << G4endl;
	ofs_ply << "property float z" << G4endl;
	ofs_ply << "property uchar red" << G4endl;
	ofs_ply << "property uchar green" << G4endl;
	ofs_ply << "property uchar blue" << G4endl;
	ofs_ply << "element face " << faceVec.size() << G4endl;
	ofs_ply << "property list uchar int vertex_index" << G4endl;
	ofs_ply << "end_header" << G4endl;

	for (G4int i=0; i<vertVec.size(); i++) {
		ofs_ply << vertVec[i].x()/cm << " " << vertVec[i].y()/cm << " " << vertVec[i].z()/cm <<
				" " << get<0>(tetColor[i]) << " " << get<1>(tetColor[i]) << " " << get<2>(tetColor[i]) << G4endl;
	}
	for (G4int i=0; i<faceVec.size(); i++) {
			ofs_ply << "3 " << faceVec[i][0] << " " << faceVec[i][1] << " " << faceVec[i][2] << G4endl;
	}
}
