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
	ofs.open("result.txt");
}


RunAction::~RunAction()
{
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

	if (isMaster) {
		ofs << "Frame " << aRun->GetRunID() << G4endl;
		ofs << "SkinFaceID \t Dose(J/kg)" << G4endl;
		G4double doseSum(0);
		for (auto itr:dose) {
			doseSum += itr.second;
			ofs << itr.first-12600 << "\t" << itr.second / (joule/kg) << G4endl;
		}
		ofs << G4endl;

		G4cout << doseSum << G4endl;
		doseSum = 1/doseSum;

		vector<G4double> faceDose;
		vector<G4double> vertexDose;
		vector<G4double> vertexWeight;
		vector<vector<G4int>> faces = tetData->GetOuterOgFaces();
		ofstream vert("vertexDoseWeight.txt");

		for (G4int i=0; i<tetData->GetSkinVertexSize(); i++) {
			G4double sumAdjacent(0);
			for (auto f:faces) {
				if (i == f[0] || i == f[1] || i == f[2]) {
					sumAdjacent += faceDose[i];
				}
			}
			vertexDose.push_back(sumAdjacent);
			vertexWeight.push_back(sumAdjacent*doseSum);
			vert << i << " " << vertexDose[i] << " " << vertexWeight[i] << G4endl;
		}

		G4String fileName = "frame" + to_string(aRun->GetRunID()) + ".ply";
		PrintPLY(fileName, vertexWeight);


//		for (G4int i=0; i<tetData->GetSkinVertexSize(); i++) {
//			G4bool pass(false);
//			for (auto itr:dose) {
//				if ((itr.first-12600) == i) { faceDose.push_back(itr.second); pass = true; break;}
//			}
//			if (pass) continue;
//			faceDose.push_back(0);
//		}
//
//		ofstream faceD("faceDose.txt");
//		for (auto itr:faceDose) faceD << itr << G4endl;
//
//		G4String fileName = "frame" + to_string(aRun->GetRunID()) + ".ply";
//		auto vertexWeight = SkinDoseDistNorm(faceDose);
//		PrintPLY(fileName, vertexWeight);
//
//
//		auto doseDist = SkinDoseDistNorm(dose);
//		PrintPLY(fileName, doseDist);
	}
}

vector<G4double> RunAction::SkinDoseDistNorm(vector<G4double> faceDose)
{
	vector<G4double> vertexWeight;
	vector<vector<G4int>> faces = tetData->GetOuterOgFaces();

	G4double sumAll(0);
	for (auto itr:faceDose) {
		sumAll += itr;
	}
	G4cout << "dose Sum: " << sumAll << G4endl;

	for(G4int i=0; i<tetData->GetSkinVertexSize(); i++) {
		G4double sumAdjcent(0);
		for (auto f:faces) {
			if (i == f[0] || i == f[1] || i == f[2]) {
				sumAdjcent += faceDose[i];
			}
		}
		vertexWeight.push_back(sumAdjcent/sumAll);
	}

	ofstream ofs("vertexWeight.txt");
	for (auto v:vertexWeight) {
		ofs << v << G4endl;
	}
	return vertexWeight;
}

void RunAction::PrintPLY(G4String fileName, vector<G4double> vertexWeight)
{
	vector<G4ThreeVector> vertVec = tetData->GetOuterVec();
	vector<vector<G4int>> faceVec = tetData->GetOuterOgFaces();
	vector<tuple<G4int,G4int,G4int>> tetColor;

	ofstream rgb("rgb_vertex.txt");
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

	G4int idx(0);
	for (auto col:tetColor) {
		rgb << idx++ << " " << get<0>(col) << " " << get<1>(col) << " " << get<2>(col) << G4endl;
	}

	system("mkdir SkinDistVertex");
	ofstream ofs("./SkinDistVertex/" + fileName);
	ofs << "ply" << G4endl;
	ofs << "format ascii 1.0" << G4endl;
	ofs << "comment exported by rapidform" << G4endl;
	ofs << "element vertex " << vertVec.size() << G4endl;
	ofs << "property float x" << G4endl;
	ofs << "property float y" << G4endl;
	ofs << "property float z" << G4endl;
	ofs << "property uchar red" << G4endl;
	ofs << "property uchar green" << G4endl;
	ofs << "property uchar blue" << G4endl;
	ofs << "element face " << faceVec.size() << G4endl;
	ofs << "property list uchar int vertex_index" << G4endl;
	ofs << "end_header" << G4endl;

	for (G4int i=0; i<vertVec.size(); i++) {
		ofs << vertVec[i].x() << " " << vertVec[i].y() << " " << vertVec[i].z() <<
				" " << get<0>(tetColor[i]) << " " << get<1>(tetColor[i]) << " " << get<2>(tetColor[i]) << G4endl;
	}
	for (G4int i=0; i<faceVec.size(); i++) {
			ofs << "3 " << faceVec[i][0] << " " << faceVec[i][1] << " " << faceVec[i][2] << G4endl;
	}

}

map<G4int, G4double> RunAction::SkinDoseDistNorm(map<G4int, G4double> dose)
{
	G4double sum(0);
	map<G4int, G4double> doseDist;
	for (auto itr:dose) {
		sum += itr.second;
	}
	G4cout << "dose Sum: " << sum << G4endl;
	for (auto itr:dose) {
		doseDist[itr.first-12600] = itr.second/sum;
		G4cout << itr.first-12600 << " " << itr.second/sum << G4endl;
	}

	return doseDist;
}

void RunAction::PrintPLY(G4String fileName, map<G4int, G4double> doseDist)
{
	vector<G4ThreeVector> vertVec = tetData->GetOuterVec();
	vector<vector<G4int>> faceVec = tetData->GetOuterOgFaces();
	map<G4int, tuple<G4int,G4int,G4int>> tetColor;

	ofstream dist("doseDist.txt");
	for (auto itr:doseDist) {
		dist << itr.first << " " << itr.second << G4endl;
	}

	ofstream rgb("rgb.txt");
	for (auto f:doseDist) {
		G4int rr, gg, bb;
		if (f.second < 0.5) {
			rr = 0;
			gg = floor(255 * f.second * 2 + 0.5);
			bb = floor(255 * (0.5 - f.second) * 2 + 0.5);
		} else {
			rr = floor(255 * (f.second) * 2 + 0.5);
			gg = floor(255 * (1-f.second) * 2 + 0.5);
			bb = 0;
		}
		tetColor[f.first] = make_tuple(rr,gg,bb);
	}

	for (auto col:tetColor) {
		rgb << col.first << " " << get<0>(col.second) << " " << get<1>(col.second) << " " << get<2>(col.second) << G4endl;
	}

	system("mkdir SkinDist");
	ofstream ofs("./SkinDist/" + fileName);
	ofs << "ply" << G4endl;
	ofs << "format ascii 1.0" << G4endl;
	ofs << "comment exported by rapidform" << G4endl;
	ofs << "element vertex " << vertVec.size() << G4endl;
	ofs << "property float x" << G4endl;
	ofs << "property float y" << G4endl;
	ofs << "property float z" << G4endl;
	ofs << "property uchar red" << G4endl;
	ofs << "property uchar green" << G4endl;
	ofs << "property uchar blue" << G4endl;
	ofs << "element face " << faceVec.size() << G4endl;
	ofs << "property list uchar int vertex_index" << G4endl;
	ofs << "end_header" << G4endl;

	for (G4int i=0; i<vertVec.size(); i++) {
		ofs << vertVec[i].x() << " " << vertVec[i].y() << " " << vertVec[i].z() << G4endl;
	}


	for (G4int i=0; i<faceVec.size(); i++) {
		G4bool pass(false);
		for (auto itr:tetColor) {
			if (itr.first == i) {
				pass = true;
				ofs << "3 " << faceVec[i][0] << " " << faceVec[i][1] << " " << faceVec[i][2] << " " <<
						get<0>(itr.second) << " " << get<1>(itr.second) << " " << get<2>(itr.second) << G4endl;
				break;
			}
		}
		if (pass) { continue; }
		else {
			ofs << "3 " << faceVec[i][0] << " " << faceVec[i][1] << " " << faceVec[i][2] << G4endl;
//					" " << 0 << " " << 0 << " " << 0 << G4endl;
		}
	}

}
