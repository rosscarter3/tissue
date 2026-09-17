//
// Filename     : degradation.cc
// Description  : Classes describing molecular degradation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#include"tissue.h"
#include"baseReaction.h"
#include"degradation.h"
#include<cmath>

namespace Degradation {
  One::
  One(std::vector<double> &paraValue, 
      std::vector< std::vector<size_t> > 
      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::One::"
		<< "One() "
		<< "Uses one parameter k_d (constant degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Degradation::One::"
		<< "One() "
		<< "Index for variable to be updated (degraded) given." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::One");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_d";
    setParameterId( tmp );
  }

  void One::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double k_d = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] -= k_d * cellData[cellI][cIndex];
    }
  }

  void One::
  derivsWithAbs(Tissue &T,
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
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {  
      double value = k_c * cellData[cellI][cIndex];
      
      cellDerivs[cellI][cIndex]  -= value;
      sdydtCell[cellI][cIndex]  += value;
    }
  }

  Two::
  Two(std::vector<double> &paraValue, 
      std::vector< std::vector<size_t> > 
      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::Two::"
		<< "Two() "
		<< "Uses one parameter k_d (degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Degradation::Two::"
		<< "Two() "
		<< "Index for variable to be updated given in first row and "
		<< "index for degradation-dependent variable in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::Two");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_d";
    setParameterId( tmp );
  }

  void Two::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    size_t xIndex = variableIndex(1,0);
    double k_d = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] -= k_d * cellData[cellI][xIndex] * cellData[cellI][cIndex]; 
    }
  }

  N::
  N(std::vector<double> &paraValue, 
    std::vector< std::vector<size_t> > 
    &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::N::"
		<< "N() "
		<< "Uses one parameter k_d (degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "Degradation::N::"
		<< "N() "
		<< "Index for variable to be updated given in first row and "
		<< "indices for degradation-dependent variables in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::N");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_d";
    setParameterId( tmp );
  }

  void Degradation::N::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double k_d = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      
      double product = k_d;
      for( size_t i=0 ; i<numVariableIndex(1) ; i++ ) { // for each rate-limiting variable
	product *= cellData[cellI][variableIndex(1,i)]; // multiply rate with its concentration
      }
      cellDerivs[cellI][cIndex] -= product*cellData[cellI][cIndex]; 
    }
  }

  Hill::Hill(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=3 ) {
      std::cerr << "Degradation::Hill::Hill() "
		<< "Uses three parameters k_deg K_hill n_hill." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Degradation::Hill::Hill() "
		<< "One variable index in first layer for variable to be updated "
		<< "and one in the second that activate degradation must be given." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Set the variable values
    setId("Degradation::Hill");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_deg";
    tmp[1] = "K_hill";
    tmp[2] = "n_hill";
    setParameterId( tmp );
  }

  void Hill::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {  
    size_t degradedIndex = variableIndex(0,0);
    size_t degraderIndex = variableIndex(1,0);
    double Kpow = std::pow(parameter(1),parameter(2));

    for (size_t cellI=0; cellI<T.numCell(); ++cellI) {
      double fac = parameter(0)*cellData[cellI][degradedIndex];
      double yPow = std::pow(cellData[cellI][degraderIndex],parameter(2));
      cellDerivs[cellI][degradedIndex] -= fac*yPow/(Kpow+yPow);
    }
  }
  
  HillN::HillN(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( indValue.size() != 3 || indValue[0].size() != 1 ) {
      std::cerr << "Degradation::HillN::HillN() "
		<< "Three levels of variable indices are used, "
		<< "the first for the updated molecule, 2nd for activators and 3rd for repressors"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( 1 + 2*(indValue[1].size()+indValue[2].size()) != paraValue.size() ) {
      std::cerr << "Degradation::HillN::HillN() "
		<< "Number of parameters does not agree with number of "
		<< "activators/repressors." << std::endl
		<< indValue[1].size() << " activators, " 
		<< indValue[2].size() << " repressors, and "
		<< paraValue.size() << " parameters\n"
		<< "One plus pairs of parameters (V_max,K_half1,n_Hill1,"
		<< "K2,n2,...) must be given." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    setId("Degradation::HillN");
    setParameter(paraValue);  
    setVariableIndex(indValue);  
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "d_max";
    size_t count=1;
    for (size_t i=0; i<indValue[1].size(); ++i) {
      tmp[count++] = "K_A";
      tmp[count++] = "n_A";
    }
    for (size_t i=0; i<indValue[2].size(); ++i) {
      tmp[count++] = "K_R";
      tmp[count++] = "n_R";
    }
    setParameterId( tmp );
  }

  void HillN::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {      
    // Do the update for each cell
    size_t numCells = T.numCell();
    size_t cIndex = variableIndex(0,0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      double contribution=parameter(0);
      size_t parameterIndex=1;
      // Activator contributions
      for( size_t i=0 ; i<numVariableIndex(1) ; i++ ) {
	double c = std::pow(cellData[cellI][variableIndex(1,i)], 
			    parameter(parameterIndex+1));
	contribution *= c
	  / ( std::pow(parameter(parameterIndex),parameter(parameterIndex+1)) + c );
	parameterIndex+=2;
      }
      // Repressor contributions
      for( size_t i=0 ; i<numVariableIndex(2) ; i++ ) {
	double c = std::pow(parameter(parameterIndex),parameter(parameterIndex+1));
	contribution *= c /
	  ( c + std::pow(cellData[cellI][variableIndex(2,i)],
			 parameter(parameterIndex+1)) );
	parameterIndex+=2;
      }   
      cellDerivs[cellI][cIndex] -= contribution*cellData[cellI][cIndex]; 
    }
  }

  TwoGeometric::
  TwoGeometric(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::TwoGeometric::"
		<< "TwoGeometric() "
		<< "Uses one parameter k_d (degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Degradation::TwoGeometric::"
		<< "TwoGeometric() "
		<< "Index for variable to be updated given in first row and "
		<< "index for degradation-dependent variable in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::TwoGeometric");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_d";
    setParameterId( tmp );
  }

  void TwoGeometric::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    size_t xIndex = variableIndex(1,0);
    double k_d = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      double cellVolume = T.cell(cellI).calculateVolume(vertexData);
      cellDerivs[cellI][cIndex] -= cellVolume* k_d * cellData[cellI][xIndex] * cellData[cellI][cIndex]; 
    }
  }

  OneWall::
  OneWall(std::vector<double> &paraValue, 
	  std::vector< std::vector<size_t> > 
	  &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::OneWall::"
		<< "OneWall() "
		<< "Uses one parameter k_cw (constant degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Degradation::OneWall::"
		<< "OneWall() "
		<< "Index for wall variable to be updated (degraded) given." << std::endl;
    exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::OneWall");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_cw";
    setParameterId( tmp );
  }

  void OneWall::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t wIndex = variableIndex(0,0);
    double k_cw = parameter(0);
    for (size_t k=0; k<numWalls; ++k) {
      wallDerivs[k][wIndex] -= k_cw * wallData[k][wIndex];
      //std::cerr << k << "\t" << wIndex << "\t" << wallDerivs[k][wIndex] << "\n";
    }
  }

  void OneWall::
  derivsWithAbs(Tissue &T,
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
    size_t numWalls = T.numWall();
    size_t wIndex = variableIndex(0,0);
    double k_cw = parameter(0);
    for (size_t k=0; k<numWalls; ++k) {
      double value = k_cw * wallData[k][wIndex];
      wallDerivs[k][wIndex] -= value;
      sdydtWall[k][wIndex] += value;
    }
  }

  OneBoundary::
  OneBoundary(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::OneBoundary::"
		<< "OneBoundary() "
		<< "Uses one parameter k_d (constant degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Degradation::OneBoundary::"
		<< "DegradationOneBoundary() "
		<< "Index for variable to be updated (degraded) given." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::OneBoundary");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_d";
    setParameterId( tmp );
  }

  void OneBoundary::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double k_d = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      
      if( T.cell(cellI).isNeighbor(T.background()) ){
	double value = k_d*cellData[cellI][cIndex];
	cellDerivs[cellI][cIndex]  -= value;
      }
    }
  }

  void OneBoundary::
  derivsWithAbs(Tissue &T,
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
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) { 
      
      if( T.cell(cellI).isNeighbor(T.background()) ){
	double value = k_c*cellData[cellI][cIndex];
	cellDerivs[cellI][cIndex]  -= value;
	sdydtCell[cellI][cIndex]  += value;
      }
    }
  }

  OneFromList::
  OneFromList(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Degradation::OneFromList::"
		<< "OneFromList() "
		<< "Uses one parameter k_d (constant degradation rate)." << std::endl;
      exit(EXIT_FAILURE);
  }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() < 1) {
      std::cerr << "Degradation::OneFromList::"
		<< "OneFromList() "
		<< "Two rows of indices needed." << std::endl
		<< "Index for variable to be updated (degraded) given in first row." << std::endl
		<< "Cell indices for cells to be updated (degraded) given in second row." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Degradation::OneFromList");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_d";
    setParameterId( tmp );
    
    numCellI = indValue[1].size();
  }

  void OneFromList::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t cIndex = variableIndex(0,0);
    double k_d = parameter(0);
    //For each cell in list of cell indices
    for (size_t cellI = 0; cellI < numCellI; ++cellI) {      
      cellDerivs[variableIndex(1,cellI)][cIndex] -= k_d * cellData[variableIndex(1,cellI)][cIndex];
    }
  }

  void OneFromList::
  derivsWithAbs(Tissue &T,
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
    //Do the update for each cell
    size_t numCells = T.numCell();
  
    size_t cIndex = variableIndex(0,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {  
      double value = k_c * cellData[cellI][cIndex];
      
      cellDerivs[cellI][cIndex]  -= value;
      sdydtCell[cellI][cIndex]  += value;
    }
  }
} // end namespace Degradation
