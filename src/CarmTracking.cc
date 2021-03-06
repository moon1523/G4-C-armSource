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
            transVec.push_back(G4ThreeVector(d,h,l));
        }
        if (dump == "Isocenter") { 
            ss >> isocenter; 
            isocenterVec.push_back(isocenter); 
        }
        if (dump == "Source")    { 
            ss >> source;    
            sourceVec.push_back(source); 
        }
    }
    ifs.close();

    G4cout << "done" << G4endl;
}