/*
 * SourceGenerator.hh
 *
 *  Created on: Jun 3, 2020
 *      Author: sungho
 */

#ifndef INCLUDE_SOURCEGENERATOR_HH_
#define INCLUDE_SOURCEGENERATOR_HH_

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class SourceGenerator
{
public:
	virtual ~SourceGenerator() {}
	virtual void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir)
private:

};







#endif /* INCLUDE_SOURCEGENERATOR_HH_ */
