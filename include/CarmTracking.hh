#ifndef CarmTracking_h
#define CarmTracking_h

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

class CarmTracking
{
    public:
        CarmTracking();
        virtual ~CarmTracking();

        void ReadTrackingResults(string fileName);

        pair<G4int, G4int> GetTemplateIndex(G4int frameNo) { return indexVec[frameNo]; }
        G4RotationMatrix GetRotationMatrix(G4int frameNo)  { return rotVec[frameNo]; }
        G4ThreeVector GetTranslationMatrix(G4int frameNo)  { return transVec[frameNo]; }
        G4ThreeVector GetIsocenterPosition(G4int frameNo)  { return isocenterVec[frameNo]; }
        G4ThreeVector GetSourcePosition(G4int frameNo)     { return sourceVec[frameNo]; }
        G4ThreeVector GetRotMColumn1(G4int frameNo)        { return rot_col1[frameNo]; }
        G4ThreeVector GetRotMColumn2(G4int frameNo)        { return rot_col2[frameNo]; }
        G4ThreeVector GetRotMColumn3(G4int frameNo)        { return rot_col3[frameNo]; }
        
    private:
        vector<G4int>             frameVec;
        vector<pair<G4int,G4int>> indexVec;
        vector<G4RotationMatrix>  rotVec;
        vector<G4ThreeVector>     transVec;
        vector<G4ThreeVector>     isocenterVec;
        vector<G4ThreeVector>     sourceVec;
        
        vector<G4ThreeVector>     rot_col1;
        vector<G4ThreeVector>     rot_col2;
        vector<G4ThreeVector>     rot_col3;

        
};

#endif
