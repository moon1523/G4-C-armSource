#include "CarmTracking.hh"

CarmTracking::CarmTracking()
{
    ReadTrackingResults("source_xyz.txt");
}

CarmTracking::~CarmTracking()
{

}

void CarmTracking::ReadTrackingResults(string fileName)
{
    G4cout << "Read " << fileName << "..." << flush;
    ifstream ifs(fileName);
    if(!ifs.is_open()) { G4cerr << "Tracking result file is not opened" << G4endl; exit(1); }

    string dump;
    while(getline(ifs, dump)) {
        stringstream ss(dump);
        ss >> dump;
        G4int frameNo, classID, templateID;
        G4double a,b,c,d,
                 e,f,g,h,
                 i,j,k,l,
                 m,n,o,p;
        G4ThreeVector isocenter, source;
        G4bool  power;
        G4double time;
        G4double rateDAP, tubeCurrent; G4int tubeVoltage;
        if (dump == "Frame")     {
            ss >> frameNo;
            frameVec.push_back(frameNo);
        }
        if (dump == "Index")     { 
            ss >> classID >> templateID; 
            indexVec.push_back(make_pair(classID,templateID)); 
        }
        if (dump == "AffineT")   {
             ifs >> a >> b >> c >> d
                 >> e >> f >> g >> h
                 >> i >> j >> k >> l
                 >> m >> n >> o >> p;

            G4ThreeVector rot1(a,e,i); G4ThreeVector rot2(b,f,j); G4ThreeVector rot3(c,g,k);

			rotVec.push_back(G4RotationMatrix(rot1, rot2, rot3));
			rot_col1.push_back(rot1);
			rot_col2.push_back(rot2);
			rot_col3.push_back(rot3);
        }
        if (dump == "Isocenter") { 
            ss >> isocenter; 
            isocenterVec.push_back(isocenter);
			transVec.push_back(isocenter);
        }
        if (dump == "Source")    { 
            ss >> source;    
            sourceVec.push_back(source); 
        }
        if (dump == "Power") {
        	ss >> power;
        	monitorPowerVec.push_back(power);
        }
        if (dump == "Time") {
        	ss >> time;
        	monitorTimeVec.push_back(time);
        }
        if (dump == "DAP_rate") {
        	ss >> rateDAP;
        	monitorDAPVec.push_back(rateDAP);
        }
        if (dump == "T_Voltage") {
        	ss >> tubeVoltage;
        	tubeVoltageVec.push_back(tubeVoltage);
        }
        if (dump == "T_Current") {
        	ss >> tubeCurrent;
        	tubeCurrentVec.push_back(tubeCurrent);
        }
    }
    ifs.close();

    G4cout << "done" << G4endl;
}
