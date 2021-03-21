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


RunAction::RunAction(TETModelImport* _tetData, G4String _output, G4Timer* _init)
: G4UserRunAction(), tetData(_tetData), outputFile(_output), initTimer(_init),
  fRun(0), ratioDAP(1.0),
  eIntensity(0), numOfEvent(0),
  monitorPower(0), monitorTime(0), monitorDAP(0),
  tubeVoltage(0), tubeCurrent(0)
{
	if (!isMaster) return;
	G4cout << "RunAction()" << G4endl;

	runTimer = new G4Timer;

	ofs.open(outputFile);
	for(G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge.push_back(0);
		dosesRunSum[i] = make_pair(0,0);
	}
}


RunAction::~RunAction()
{
	if(!isMaster) return;
	G4cout << "~RunAction() !!!" << G4endl;

	G4double maxDose(0), maxRelativeErr; G4int maxID(0);
	for (auto itr:dosesRunSum) {
		if (maxDose < itr.second.first) {
			maxDose = itr.second.first;
			maxRelativeErr = itr.second.second;
			maxID = itr.first;
		}
	}

	ofs << "Maximum skin dose: " <<  maxID << "\t" << scientific << maxDose/(joule/kg) << "\t" << maxRelativeErr << G4endl;

	G4double sumvertDose(0);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		sumvertDose += vertexDose2WeightMerge[i];
//		if (maximumDose < vertexDose2WeightMerge[i]) maximumDose = vertexDose2WeightMerge[i];
	}



	G4double weightMin(DBL_MAX);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] /= sumvertDose;
		if ((vertexDose2WeightMerge[i] > 0) && (vertexDose2WeightMerge[i] < weightMin))
			weightMin = vertexDose2WeightMerge[i];
	}

	G4double logMax(-DBL_MAX);
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] += weightMin;
		vertexDose2WeightMerge[i] /= weightMin;
		vertexDose2WeightMerge[i] = log10(vertexDose2WeightMerge[i]);
		if (logMax < vertexDose2WeightMerge[i]) logMax = vertexDose2WeightMerge[i];
	}
	for (G4int i=0;i<tetData->GetSkinVertexSize(); i++) {
		vertexDose2WeightMerge[i] /= logMax;
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
	G4cout << "BeginOfRunAction()" << G4endl;
	numOfEvent = aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(int(numOfEvent*0.1));

	if(isMaster) {
		initTimer->Stop();
		runTimer->Start();
	}

	const PrimaryGeneratorAction* primary =
			dynamic_cast<const PrimaryGeneratorAction*>
			(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	if(!primary) return;

	monitorPower = primary->GetMonitorPower();
	monitorTime  = primary->GetMonitorTime();
	monitorDAP   = primary->GetMonitorDAP();
	tubeVoltage  = primary->GetTubeVoltage();
	tubeCurrent  = primary->GetTubeCurrent();
	eIntensity   = primary->GetEIntensity();
	fRun->SetPrimary(eIntensity, monitorPower, monitorTime, monitorDAP, tubeVoltage, tubeCurrent);
}


void RunAction::EndOfRunAction(const G4Run* aRun)
{
	if(!isMaster) return;
	runTimer->Stop();
	G4cout << "EndOfRunAction()" << G4endl;


	const DetectorConstruction* detectorConstruction =
			static_cast<const DetectorConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());


	const Run* runE = static_cast<const Run*>(aRun);
	std::vector<G4VPhysicalVolume*> phyVec = detectorConstruction->GetScoringPV();

	G4double massSum(0);
	for (G4VPhysicalVolume* phy:phyVec) {
		G4double mass = phy->GetLogicalVolume()->GetMass();
		G4int idx = phy->GetCopyNo();
		massMap[idx] += mass;
		massSum += mass;
	}
	map<G4int, G4double> dosechk;

	// Call Source Information
	eIntensity   = fRun->GetEIntensty(); // #/cm2-mAs
	monitorPower = fRun->GetMonitorPower();
	monitorTime  = fRun->GetMonitorTime();
	monitorDAP   = fRun->GetMonitorDAP();
	tubeVoltage  = fRun->GetTubeVoltage();
	tubeCurrent  = fRun->GetTubeCurrent();


	G4double h = 100; // 1 meter
	G4double a = h * tan(8.216*deg);
	G4double b = h * tan(10.565*deg);
	G4double area1m = 4*a*b; // 1077.202 cm2
	G4double massDAP = detectorConstruction->GetDAPPhysicalVolume()->GetLogicalVolume()->GetMass(); //
	G4double areaDAP = detectorConstruction->GetDAParea(); // 96.94822 cm2
	G4double intensity = eIntensity * area1m * tubeCurrent;

	// Set Doses
	doses.clear();
	auto edepMap = runE->GetEdepMap();
	G4double doseSum(0); G4int psdID(0); G4double psdDose(0), psdrelError(0);
	for (auto itr:edepMap) {
		G4double meanDose   = itr.second.first / numOfEvent;
		G4double squareDose = itr.second.second / numOfEvent;
		G4double variance   = ((squareDose) - (meanDose * meanDose)) / numOfEvent;
		G4double relativeE  = sqrt(variance)/meanDose;

		if (itr.first == 1000) { // DAP
			doses[itr.first] = make_pair(meanDose/massDAP*intensity*areaDAP, relativeE);
			continue;
		}
		doses[itr.first] = make_pair(meanDose/massMap[itr.first]*intensity, relativeE);
		dosechk[itr.first] = meanDose/massMap[itr.first];
		doseSum += doses[itr.first].first;
		if (psdDose < doses[itr.first].first) {
			psdDose = doses[itr.first].first;
			psdrelError = doses[itr.first].second;
			psdID = itr.first - 12600;
		}
		dosesRunSum[itr.first].first  += meanDose/massMap[itr.first]*intensity;
		dosesRunSum[itr.first].second += relativeE;
	}
	doseSum /= massSum;
	ratioDAP = monitorDAP / (doses[1000].first / (joule/kg) / (cm*cm));

	for (auto itr:dosesRunSum) {
		itr.second.first *= ratioDAP;
	}

	// Print Results
	ofs << "# Frame " << detectorConstruction->GetFrameNo() << " ============================================" << G4endl;
	ofs << "NPS                          : " << numOfEvent << G4endl;
	ofs << "init/run Time                : " << initTimer->GetRealElapsed() << "\t" << runTimer->GetRealElapsed() << G4endl;
	ofs << "Intensity [#/s]              : " << scientific << intensity << G4endl;
	ofs << "Power                        : " << monitorPower << G4endl;
	ofs << "Time                         : " << monitorTime << G4endl;
	ofs << "Tube Voltage/Current [kVp,mA]: " << tubeVoltage << "\t" << tubeCurrent << G4endl;
	ofs << "DAP ratio                    : " << ratioDAP << G4endl;
	ofs << "Monitored DAP_rate [Gy-cm2/s]: " << monitorDAP << G4endl;
	ofs << "Estimated DAP_rate [Gy-cm2/s]: " << doses[1000].first / (joule/kg) / (cm*cm) << G4endl;
	ofs << "PSD rate [FaceID,Gy/s,relE]  : " << psdID << "\t" << psdDose * ratioDAP / (joule/kg) << "\t" << psdrelError << G4endl;
	ofs << "Skin dose_rate [Gy-cm2/s]    : " << doseSum * ratioDAP / (joule/kg) << G4endl;

	ofs << "-- Skin Dose Distribution ----------------------------" << G4endl;
	ofs << "SkinFaceID\tDose rate (Gy/s)\tRel.error" << G4endl;
	for (auto itr:doses) {
		if (itr.first == 1000) continue;
		ofs <<  setw(10) << fixed <<  itr.first-12600 << "\t"
			<<  setw(15) << scientific << itr.second.first * ratioDAP / (joule/kg) << "\t"
			<<  setw(10) << fixed << itr.second.second << G4endl;
	}

	ofs << "## Dose Check " << detectorConstruction->GetFrameNo() << G4endl;
	for (auto itr:dosechk) {
		if (itr.first == 1000) { ofs << itr.first 		  << " " << scientific << itr.second / (joule/kg) << G4endl; continue; }
							     ofs << itr.first - 12600 << " " << scientific << itr.second / (joule/kg) << G4endl;
	}
	ofs << G4endl;


	// Set vertex weight to print PLY.
	vector<G4double> faceDose;
	vector<G4double> vertexDose2Weight;
	vertVec = tetData->GetOuterVec();
	faceVec = tetData->GetOuterOgFaces();

	for (G4int i=0; i<tetData->GetSkinFaceSize(); i++) {
		G4bool pass(false);
		for (auto itr:doses) {
			if (itr.first == 1000) continue;
			if ((itr.first-12600) == i) { faceDose.push_back(itr.second.first); pass = true; break;}
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
	G4String fileName = "frame" + to_string(detectorConstruction->GetFrameNo()) + ".ply";
	PrintPLY(fileName, vertexDose2Weight);

	initTimer->Start();
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

	ofstream ofs_ply("./skinPLY/" + fileName);
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
