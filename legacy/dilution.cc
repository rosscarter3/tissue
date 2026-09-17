//
// Filename     : dilution.cc
// Description  : Classes describing dilution of moleculecular concentrations due to growth
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//
#include <vector>
#include <cmath>
#include"dilution.h"
#include"baseReaction.h"


Dilution::FromVertexDerivs::
FromVertexDerivs(std::vector<double> &paraValue,
		 std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size()) {
    std::cerr << "Dilution::FromVertexDerivs::Dilution::FromVertexDerivs() "
	      << "Uses no parameters." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() < 1) {
    std::cerr << "Dilution::FromVertexDerivs::Dilution::FromVertexDerivs() "
	      << "List of concentration variable index must be given in "
	      << "first level." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("Dilution::FromVertexDerivs");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  setParameterId(tmp);
}

void Dilution::FromVertexDerivs::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension;
  dimension = vertexData[0].size();
  assert(dimension==2);
  
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    double area = cell.calculateVolume(vertexData,1);
    
    double areaDerivs=0.0;
    for( size_t k=0 ; k<cell.numVertex() ; ++k ) {
      size_t vI = cell.vertex(k)->index();
      size_t vIPlus = cell.vertex((k+1)%(cell.numVertex()))->index();
      areaDerivs += vertexData[vIPlus][1]*vertexDerivs[vI][0] - 
	vertexData[vI][1]*vertexDerivs[vIPlus][0] -
	vertexData[vIPlus][0]*vertexDerivs[vI][1] +
	vertexData[vI][0]*vertexDerivs[vIPlus][1];
    }
    areaDerivs *= 0.5;
    double fac = areaDerivs/area;
    for (size_t k=0; k<numVariableIndex(0); ++k)
      cellDerivs[n][variableIndex(0,k)] -= cellData[n][variableIndex(0,k)]*fac;
  }
}


void Dilution::FromVertexDerivs::derivsWithAbs(Tissue &T,
					       DataMatrix &cellData,
					       DataMatrix &wallData,
					       DataMatrix &vertexData,
					       DataMatrix &cellDerivs,
					       DataMatrix &wallDerivs,
					       DataMatrix &vertexDerivs,
					       DataMatrix &sdydtCell,
					       DataMatrix &sdydtWall,
					       DataMatrix &sdydtVertex ) 
{
  size_t dimension;
  dimension = vertexData[0].size();
  assert(dimension==2);
  
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    double area = cell.calculateVolume(vertexData,1);
    
    double areaDerivs=0.0;
    for( size_t k=0 ; k<cell.numVertex() ; ++k ) {
      size_t vI = cell.vertex(k)->index();
      size_t vIPlus = cell.vertex((k+1)%(cell.numVertex()))->index();
      areaDerivs += vertexData[vIPlus][1]*vertexDerivs[vI][0] - 
	vertexData[vI][1]*vertexDerivs[vIPlus][0] -
	vertexData[vIPlus][0]*vertexDerivs[vI][1] +
	vertexData[vI][0]*vertexDerivs[vIPlus][1];
    }
    areaDerivs *= 0.5;
    double fac = areaDerivs/area;
    for (size_t k=0; k<numVariableIndex(0); ++k)
      cellDerivs[n][variableIndex(0,k)] -= cellData[n][variableIndex(0,k)]*
	fac;
  }
}

