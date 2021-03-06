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
// TETModelImport.cc
// \file   MRCP_GEANT4/External/src/TETModelImport.cc
// \author Haegin Han
//

#include "TETModelImport.hh"

TETModelImport::TETModelImport(G4String _phantomName, G4UIExecutive* ui)
:doseOrganized(false)
{
	G4cout << "TETModelImport()" << G4endl;
	// set phantom name
	phantomName = _phantomName;

	G4cout << "================================================================================"<<G4endl;
	G4cout << "\t" << phantomName << " was implemented in this CODE!!   "<< G4endl;
	G4cout << "================================================================================"<<G4endl;

	G4String eleFile      =  phantomName + ".ele";
	G4String nodeFile     =  phantomName + ".node";
	G4String materialFile =  phantomName + ".material";
	G4String doseFile     =  phantomName + ".dose";
	G4String boneFile     =  phantomName + ".RBMnBS";
	G4String drfFile      =  phantomName + ".DRF";
	G4String mtlFile      =  phantomName + ".mtl";

	// read dose file (*.dose) -if there is any
	DoseRead(doseFile);
	// read phantom data files (*. ele, *.node)
    DataRead(eleFile, nodeFile);
	// read material file (*.material)
	MaterialRead(materialFile);
	// read bone file (*.RBMnBS)
	RBMBSRead(boneFile);
	// read bone file (*.DRF)
	DRFRead(drfFile);
	// read colour data file (colour.dat) if this is interactive mode
	if(ui) ColourRead(mtlFile);
	// print the summary of phantom information
	PrintMaterialInfomation();
}


void TETModelImport::DoseRead(G4String doseFile){
	//read dose file : PLEASE be careful not to include dose ID 0
	std::ifstream ifs(doseFile);
	if(!ifs.is_open()) return;
	doseOrganized = true;

	G4String aLine;
	while(!ifs.eof()){
		getline(ifs, aLine);
		if(aLine.empty()) break;

		std::stringstream ss(aLine);
		G4int doseID; ss>>doseID;
		G4String name; ss>>name; doseName[doseID] = name;
		G4int organID;
		while(ss>>organID){
			if(organ2dose.find(organID)==organ2dose.end()) organ2dose[organID] = {doseID};
			else	                                       organ2dose[organID].push_back(doseID);
		}
	}
	ifs.close();

}

void TETModelImport::DataRead(G4String eleFile, G4String nodeFile)
{
	using namespace std;

	G4String tempStr;
	G4int tempInt;

	// Read *.node file
	//
	std::ifstream ifpNode;

	ifpNode.open(nodeFile);
	if(!ifpNode.is_open()) {
		// exception for the case when there is no *.node file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + nodeFile ).c_str());
	}
	G4cout << "  Opening TETGEN node (vertex points: x y z) file '" << nodeFile << "'" <<G4endl;

	G4int numVertex;
	G4double xPos, yPos, zPos;
	G4double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
	G4double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

	ifpNode >> numVertex >> tempInt >> tempInt >> tempInt;

	for(G4int i=0; i<numVertex; i++)
	{
		ifpNode >> tempInt >> xPos >> yPos >> zPos;

		// set the unit
		xPos*=cm;
		yPos*=cm;
		zPos*=cm;

		// save the node data as the form of std::vector<G4ThreeVector>
		vertexVector.push_back(G4ThreeVector(xPos, yPos, zPos));

		// to get the information of the bounding box of phantom
		if (xPos < xMin) xMin = xPos;
		if (xPos > xMax) xMax = xPos;
		if (yPos < yMin) yMin = yPos;
		if (yPos > yMax) yMax = yPos;
		if (zPos < zMin) zMin = zPos;
		if (zPos > zMax) zMax = zPos;
	}

	// set the variables for the bounding box and phantom size
	boundingBox_Min = G4ThreeVector(xMin,yMin,zMin);
	boundingBox_Max = G4ThreeVector(xMax,yMax,zMax);
	G4ThreeVector center = (boundingBox_Max+boundingBox_Min)*0.5;
	phantomSize = G4ThreeVector(xMax-xMin,yMax-yMin,zMax-zMin);

	ifpNode.close();

	// Read *.ele file
	//
	std::ifstream ifpEle;

	ifpEle.open(eleFile);
	if(!ifpEle.is_open()) {
		// exception for the case when there is no *.ele file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + eleFile ).c_str());
	}
	G4cout << "  Opening TETGEN elements (tetrahedron with node No.) file '" << eleFile << "'" <<G4endl;

	G4int numEle;
	ifpEle >> numEle  >> tempInt >> tempInt;

	for(G4int i=0; i<numEle; i++)
	{
		ifpEle >> tempInt;
		G4int* ele = new G4int[4];
		for(G4int j=0;j<4;j++){
			ifpEle >> tempInt;
			ele[j]=tempInt;
		}
		eleVector.push_back(ele);
		ifpEle >> tempInt;
		materialVector.push_back(tempInt);

		// save the element (tetrahedron) data as the form of std::vector<G4Tet*>
		tetVector.push_back(new G4Tet("Tet_Solid",
							   		  vertexVector[ele[0]]-center,
									  vertexVector[ele[1]]-center,
									  vertexVector[ele[2]]-center,
									  vertexVector[ele[3]]-center));

		// calculate the total volume and the number of tetrahedrons for each organ
		std::map<G4int, G4double>::iterator FindIter = volumeMap.find(materialVector[i]);

		if(FindIter!=volumeMap.end()){
			FindIter->second += tetVector[i]->GetCubicVolume();
			numTetMap[materialVector[i]]++;
		}
		else {
			volumeMap[materialVector[i]] = tetVector[i]->GetCubicVolume();
			numTetMap[materialVector[i]] = 1;
		}
	}

	ifstream ifs_inner("skin100.obj");
	ifstream ifs_outer("skin50.obj");
	if (!ifs_inner.is_open()) { G4cerr << "skin100.obj was not opened" << G4endl; exit(1); }
	if (!ifs_outer.is_open()) { G4cerr << "skin50.obj was not opened" << G4endl; exit(1); }
	G4String dump;
	G4cout << "Reading skin100.obj file..." << flush;
	while (getline(ifs_inner, dump)) {
		stringstream ss(dump);
		ss >> dump;
		if (dump == "v") {
			double x, y, z;
			ss >> x >> y >> z; G4ThreeVector point(x*cm,y*cm,z*cm);
			innerVec.push_back(point);
		}
		if (dump == "f") {
			G4int a, b, c;
			ss >> a >> b >> c;
			innerFaces.push_back({a-1, b-1, c-1});
			innerOgFaces.push_back({a-1, b-1, c-1});
		}
	} ifs_inner.close();
	G4cout << "done" << G4endl;

	G4cout << "Reading skin50.obj file..." << flush;
	pictess = new G4TessellatedSolid();
	while (getline(ifs_outer, dump)) {
		stringstream ss(dump);
		ss >> dump;
		if (dump == "v") {
			double x, y, z;
			ss >> x >> y >> z; G4ThreeVector point(x*cm,y*cm,z*cm);
			outerVec.push_back(point);
		}
		if (dump == "f") {
			G4int a, b, c;
			ss >> a >> b >> c;
			outerFaces.push_back({a-1, b-1, c-1});
			outerOgFaces.push_back({a-1, b-1, c-1});
			pictess->AddFacet(new G4TriangularFacet(outerVec[a-1], outerVec[b-1], outerVec[c-1], ABSOLUTE));
		}
	} pictess->SetSolidClosed(true); ifs_outer.close();
	G4cout << "done" << G4endl;


	map<IJK, vector<G4int>> grid;
	for (G4int i=0;i<vertexVector.size();i++) {
		grid[IJK(floor(vertexVector[i].x()),floor(vertexVector[i].y()),floor(vertexVector[i].z()))].push_back(i);
	}

	double epsilon(1e-5);

	for (auto i:innerVec) {
		IJK inner_ijk = IJK(floor(i.x()), floor(i.y()), floor(i.z()));
		for (auto n:grid[inner_ijk]) {
			if (fabs(i.x()-vertexVector[n].x()) > epsilon) continue;
			if (fabs(i.y()-vertexVector[n].y()) > epsilon) continue;
			if (fabs(i.z()-vertexVector[n].z()) > epsilon) continue;
			innerNodes.push_back(n);
			break;
		}
	} assert(innerNodes.size()==innerVec.size());

	for (auto o:outerVec) {
		IJK outer_ijk = IJK(floor(o.x()), floor(o.y()), floor(o.z()));
		for (auto n:grid[outer_ijk]) {
			if (fabs(o.x()-vertexVector[n].x()) > epsilon) continue;
			if (fabs(o.y()-vertexVector[n].y()) > epsilon) continue;
			if (fabs(o.z()-vertexVector[n].z()) > epsilon) continue;
			outerNodes.push_back(n);
			break;
		}
	} assert(outerNodes.size()==outerVec.size());

	for(vector<G4int> &aFace:innerFaces) {
		aFace = {innerNodes[aFace[0]],innerNodes[aFace[1]],innerNodes[aFace[2]]};
	}
	for(vector<G4int> &aFace:outerFaces) {
		aFace = {outerNodes[aFace[0]],outerNodes[aFace[1]],outerNodes[aFace[2]]};
	}

	ifpEle.close();
}

void TETModelImport::MaterialRead(G4String materialFile)
{
	// Read material file (*.material)
	//
	std::ifstream ifpMat;

	ifpMat.open(materialFile);
	if(!ifpMat.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + materialFile ).c_str());
	}

	G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4String MaterialName;
	G4double density;

	while(!ifpMat.eof())
	{
		ifpMat >> read_data;                   //ex) 'C' RBM
		ifpMat >> MaterialName;                //ex)  C 'RBM'
		ifpMat >> read_data;
		density = std::atof(read_data);        //ex) 1.30
		ifpMat >> read_data;                   //ex) g/cm3
		ifpMat >> read_data;
        if(G4String(read_data).empty()) continue;
		token = std::strtok(read_data,"m");
		G4int matID = std::atoi(token);        //ex) m'10'
        materialIndex.push_back(matID);
		organNameMap[matID]= MaterialName;
		densityMap[matID] = density*g/cm3;

		for(G4int i=0 ;  ; i++)
		{
			ifpMat >> read_data;
			if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

			zaid = (G4int)(std::atoi(read_data)/1000);
			ifpMat >> read_data;
			fraction = -1.0 * std::atof(read_data);
			materialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
		}
	}
	ifpMat.close();

	// Construct materials for each organ
	//
    G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
    G4NistManager* nistManager = G4NistManager::Instance();

	for(G4int i=0;i<(G4int)materialIndex.size();i++){
		G4int idx = materialIndex[i];
		G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
			if(materialIndexMap[idx][j].first==1) mat->AddElement(elH, materialIndexMap[idx][j].second);
			else mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
		}
		materialMap[idx]=mat;
		massMap[idx]=densityMap[idx]*volumeMap[idx];
	}

	if(DoseWasOrganized()){
		for(auto dm:doseName){
			doseMassMap[dm.first] = 0;
		}
		for(auto od:organ2dose){
			for(auto doseID:od.second)
				doseMassMap[doseID] += massMap[od.first];
		}
	}
}

void TETModelImport::RBMBSRead(G4String bonefile){
	std::ifstream ifs(bonefile);
	if(!ifs.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::RBMBSRead","",JustWarning,
				G4String("      There is no " + bonefile ).c_str());
		return;
	}
	G4int idx;
	G4double rbm, bs;
	while(ifs>>idx>>rbm>>bs){
        if(rbmRatio.find(idx)!=rbmRatio.end()) {
            G4cerr<<idx<<" is duplicated in RBMBS file.."<<G4endl;
            exit(0);
        }
		rbmRatio[idx]=rbm;
		bsRatio[idx]=bs;
	}
}

void TETModelImport::DRFRead(G4String DRFfile){
	std::ifstream ifp;
	ifp.open(DRFfile.c_str());

	if(!ifp.is_open()) {
		G4cerr << DRFfile << " not found!!" << G4endl;
		return;
	}

	G4int ID;
    G4double DRF;
    while (!ifp.eof()) {
        G4String dump;
        getline(ifp, dump);
        std::stringstream ss(dump); dump.clear();
        ss >> dump;
        if(dump.empty()) continue;
        ID = atoi(dump.c_str());
        if(rbmDRF.find(ID)!=rbmDRF.end()) {
            G4cerr<<ID<<" is duplicated in DRF file.."<<G4endl;
            exit(0);
        }
        rbmDRF[ID]={};
        bsDRF[ID]={};
    	for (int j=0; j<25; j++) {
            ss >> DRF;
    		rbmDRF[ID].push_back(DRF);
    	}
    	for (int j=0; j<25; j++) {
            ss >> DRF;
    		bsDRF[ID].push_back(DRF);
    	}
    }
    ifp.close();
}

void TETModelImport::ColourRead(G4String mtlFile)
{
	// Read colour data file (colour.dat)
	//
	std::ifstream ifpColour;

	ifpColour.open(mtlFile);
	if(!ifpColour.is_open()) {
		// exception for the case when there is no colour.dat file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String(mtlFile + " file was not found ").c_str());
	}

	G4cout << "  Opening colour data file "+mtlFile <<G4endl;

	G4int organID;
	G4ThreeVector rgb;
	G4String dump;
	while( ifpColour >> dump ){
		if(dump=="newmtl"){
			ifpColour>>dump;
			organID=GetID(dump);
		}
		else if(dump=="Kd"){
			ifpColour>>rgb;
			colourMap[organID] = G4Colour(rgb.getX(), rgb.getY(), rgb.getZ(), 0.);
		}
	}
	ifpColour.close();
	colourMap[160].SetAlpha(0.1);
	colourMap[140].SetAlpha(0.1);
	colourMap[141].SetAlpha(0.1);
	colourMap[142].SetAlpha(0.1);
}

void TETModelImport::PrintMaterialInfomation()
{
	// Print the overall information for each organ
	//
	G4cout << G4endl
		   << std::setw(9)  << "Organ ID"
		   << std::setw(11) << "# of Tet"
		   << std::setw(11) << "vol [cm3]"
		   << std::setw(11) << "d [g/cm3]"
		   << std::setw(11) << "mass [g]"
		   << "\t" << "organ/tissue"<< G4endl ;
	G4cout << "--------------------------------------------------------------------------------"<<G4endl;

	std::map<G4int, G4Material*>::iterator matIter;
	G4cout<<std::setiosflags(std::ios::fixed);
	G4cout.precision(3);
	for(matIter=materialMap.begin(); matIter!=materialMap.end();matIter++)
	{
		G4int idx = matIter->first;

		G4cout << std::setw(9)  << idx                         // organ ID
			   << std::setw(11) << numTetMap[idx]              // # of tetrahedrons
			   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
			   << std::setw(11) << materialMap[idx]
			                       ->GetDensity()/(g/cm3)      // organ density
			   << std::setw(11) << massMap[idx]/g              // organ mass
			   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
	}
}
