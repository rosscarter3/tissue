//
// Filename     : initiation.cc
// Description  : Classes describing initiation rules for a tissue simulation
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : August 2019
// Revision     : $Id:$
//
#include"tissue.h"
#include "baseReaction.h"
#include "initiation.h"
#include<cstdlib>

namespace Initiation {

  RandomBoolean::
  RandomBoolean(std::vector<double> &paraValue,
       std::vector< std::vector<size_t> >
       &indValue )
  {
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size() != 1 && paraValue.size() != 2) {
      std::cerr << "Initiation::RandomBoolean::"
                << "RandomBoolean() "
                << "Uses one parameter p (probability to set cell variable to one)." << std::endl
		<< "A second (long) parameter can be given to set the seed of the random generator." << std::endl
		<< "(Otherwise the seed is randomized)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[0]<0.0 || paraValue[0]>1.0) {
      std::cerr << "Initiation::RandomBoolean::"
                << "RandomBoolean() "
                << "First parameter a probability and have to be in [0:1]." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Initiation::RandomBoolean::"
                << "RandomBoolean() "
                << "Index for cell variable to be initiated given." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Initiation::RandomBoolean");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "p";
    if (paraValue.size()==2)
      tmp[1] = "seed";
    setParameterId( tmp );
  }

  void RandomBoolean::
  derivs(Tissue &T,
         DataMatrix &cellData,
         DataMatrix &wallData,
         DataMatrix &vertexData,
         DataMatrix &cellDerivs,
         DataMatrix &wallDerivs,
         DataMatrix &vertexDerivs )
  {
    // nothing
  }

  void RandomBoolean::
  initiate(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs)
  {
    // Randomize or initiate with given seed
    long int idum=0;
    if (numParameter()==1)
      idum = myRandom::ran3Randomize();
    else
      idum = long(parameter(1));
    myRandom::sran3(idum);

    //Do the initiation for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double prob = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      cellData[cellI][cIndex] = 0.0;
      if (myRandom::ran3()<prob)
	cellData[cellI][cIndex] = 1.0;
    }
  }
  
} // end namespace Initiation
  
