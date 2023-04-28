//
// Filename     : creation.cc
// Description  : Classes describing molecular production/creation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#include"tissue.h"
#include"baseReaction.h"
#include"creation.h"
#include<cmath>

namespace Creation {
  Zero::
  Zero(std::vector<double> &paraValue, 
       std::vector< std::vector<size_t> > 
       &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Creation::Zero::"
		<< "Zero() "
		<< "Uses one parameter k_c (constant production rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Creation::Zero::"
		<< "Zero() "
		<< "Index for variable to be updated given." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Creation::Zero");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_c";
    setParameterId( tmp );
  }
  
  void Zero::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] += k_c;
    }
  }
  
  void Zero::
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
      cellDerivs[cellI][cIndex] += k_c;
      sdydtCell[cellI][cIndex] += k_c;
    }
  }

  One::
  One(std::vector<double> &paraValue, 
      std::vector< std::vector<size_t> > 
      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Creation::One::"
	      << "One() "
		<< "Uses one parameter k_c (linear production rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Creation::One::"
		<< "One() "
		<< "Index for variable to be updated given in first row and "
		<< "index for production-dependent variable in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Creation::One");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_c";
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
    size_t xIndex = variableIndex(1,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] += k_c * cellData[cellI][xIndex];
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
    size_t xIndex = variableIndex(1,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      double value = k_c * cellData[cellI][xIndex];
      cellDerivs[cellI][cIndex] += value;
      sdydtCell[cellI][cIndex] += value;    
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
      std::cerr << "Creation::Two::"
		<< "Two() "
		<< "Uses one parameter k_c (linear production rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 2 ) {
      std::cerr << "Creation::Two::"
		<< "Two() "
		<< "One index for variable to be updated given in first row and "
		<< "Two indices for production-dependent variables in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Creation::Two");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_c";
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
    size_t yIndex = variableIndex(1,1);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] += k_c * cellData[cellI][xIndex] * cellData[cellI][yIndex];
    }
  }

  Three::
  Three(std::vector<double> &paraValue, 
      std::vector< std::vector<size_t> > 
      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Creation::Three::"
		<< "Three() "
		<< "Uses one parameter k_c (linear production rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 3 ) {
      std::cerr << "Creation::Three::"
		<< "Three() "
		<< "One index for variable to be updated given in first row and "
		<< "three indices for production-dependent variables in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Creation::Three");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_c";
    setParameterId( tmp );
  }

  void Three::
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
    size_t yIndex = variableIndex(1,1);
    size_t zIndex = variableIndex(1,2);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] += k_c * cellData[cellI][xIndex] * cellData[cellI][yIndex] *
	cellData[cellI][zIndex];
    }
  }

  SpatialSphere::
  SpatialSphere(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=4 ) {
      std::cerr << "Creation::SpatialSphere::SpatialSphere() "
		<< "Uses four parameters V_max R(K_Hill) n_Hill and R_sign."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Creation::SpatialSphere::"
		<< "SpatialSphere() "
		<< "Index for variable to be updated given." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Sign should be -/+1
    if( paraValue[3] != -1 && paraValue[3] != 1 ) {
      std::cerr << "Creation::SpatialSphere::SpatialSphere() "
		<< "R_sign should be +/-1 to set production inside/outside R"
		<< std::endl;
      exit(EXIT_FAILURE);
    }	
    // Set the variable values
    setId("Creation::SpatialSphere");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_max";
    tmp[1] = "R (K_Hill)";
    tmp[2] = "n_Hill";
    tmp[3] = "R_sign";  
    setParameterId( tmp );
  }

  void SpatialSphere::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double powK_ = std::pow(parameter(1),parameter(2));
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      
      //Calculate cell center from vertices positions
      std::vector<double> cellCenter;
      cellCenter = T.cell(cellI).positionFromVertex(vertexData);
      assert( cellCenter.size() == vertexData[0].size() );
      double r=0.0;
      for( size_t d=0 ; d<cellCenter.size() ; ++d )
	r += cellCenter[d]*cellCenter[d];
      r = std::sqrt(r);
      
      double powR = std::pow(r,parameter(2));
      
      if (parameter(3)>0.0)
	cellDerivs[cellI][cIndex] += parameter(0)*powR/(powK_+powR);
      else
	cellDerivs[cellI][cIndex] += parameter(0)*powK_/(powK_+powR);
    }
  }

  SpatialCylinder::
  SpatialCylinder(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=4 ) {
      std::cerr << "Creation::SpatialCylinder::SpatialCylinder() "
		<< "Uses four parameters V_max R(K_Hill) n_Hill and R_sign."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Creation::SpatialCylinder::"
		<< "SpatialCylinder() "
		<< "Index for variable to be updated given." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Sign should be -/+1
    if( paraValue[3] != -1 && paraValue[3] != 1 ) {
      std::cerr << "Creation::SpatialCylinder::SpatialCylinder() "
		<< "R_sign should be +/-1 to set production inside/outside R"
		<< std::endl;
      exit(EXIT_FAILURE);
    }	
    // Set the variable values
    setId("Creation::SpatialCylinder");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_max";
    tmp[1] = "R (K_Hill)";
    tmp[2] = "n_Hill";
    tmp[3] = "R_sign";  
    setParameterId( tmp );
  }

  void SpatialCylinder::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double powK_ = std::pow(parameter(1),parameter(2));
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      
      //Calculate cell center from vertices positions
      std::vector<double> cellCenter;
      cellCenter = T.cell(cellI).positionFromVertex(vertexData);
      assert( cellCenter.size() == vertexData[0].size() );
      assert( cellCenter.size() > 1 );
      double r=0.0;
      for( size_t d=0 ; d<2 ; ++d )
	r += cellCenter[d]*cellCenter[d];
      r = std::sqrt(r);
      
      double powR = std::pow(r,parameter(2));
      
      if (parameter(3)>0.0)
	cellDerivs[cellI][cIndex] += parameter(0)*powR/(powK_+powR);
      else
	cellDerivs[cellI][cIndex] += parameter(0)*powK_/(powK_+powR);
    }
  }
  
  SpatialRing::
  SpatialRing(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=5 ) {
      std::cerr << "Creation::SpatialRing::SpatialRing() "
		<< "Uses five parameters V_max R(K_Hill) r_ring n_Hill and R_sign"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Creation::SpatialRing::"
		<< "SpatialRing() "
		<< "Index for variable to be updated given." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Sign should be -/+1
    if( paraValue[4] != -1 && paraValue[4] != 1 ) {
      std::cerr << "Creation::SpatialRing::SpatialRing() "
		<< "R_sign should be +/-1 to set production inside/outside R"
		<< std::endl;
      exit(EXIT_FAILURE);
    }	
    // Set the variable values
    setId("Creation::SpatialRing");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_max";
    tmp[1] = "R (K_Hill)";
    tmp[2] = "r_ring";
    tmp[3] = "n_Hill";
    tmp[4] = "R_sign";  
    setParameterId( tmp );
  }

  void SpatialRing::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    double powK_ = std::pow(parameter(1),parameter(3));
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      
      //Calculate cell center from vertices positions
      std::vector<double> cellCenter;
      cellCenter = T.cell(cellI).positionFromVertex(vertexData);
      assert( cellCenter.size() == vertexData[0].size() );
      double r=0.0;
      for( size_t d=0 ; d<cellCenter.size() ; ++d )
	r += cellCenter[d]*cellCenter[d];
      r = std::sqrt(r);
      
      r = std::sqrt( (r - parameter(2)) * (r - parameter(2))  ); // we want the distance to the ring
      
      double powR = std::pow(r,parameter(3));
      
      if (parameter(4)>0.0)
	cellDerivs[cellI][cIndex] += parameter(0)*powR/(powK_+powR);
      else
	cellDerivs[cellI][cIndex] += parameter(0)*powK_/(powK_+powR);
    }
  }

  SpatialCoordinate::
  SpatialCoordinate(std::vector<double> &paraValue, 
		    std::vector< std::vector<size_t> > 
		    &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=4 ) {
      std::cerr << "Creation::SpatialCoordinate::SpatialCoordinate() "
		<< "Uses four parameters V_max X(K_Hill) n_Hill and X_sign"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Creation::SpatialCoordinate::"
		<< "SpatialCoordinate() "
		<< "Two levels of indices used: index for variable to be updated given as first index, "
		<< "index for the spatial coordinate to use given as the second index" << std::endl;
      exit(EXIT_FAILURE);
    }
    // Sign should be -/+1
    if( paraValue[3] != -1 && paraValue[3] != 1 ) {
      std::cerr << "Creation::SpatialCoordinate::SpatialCoordinate() "
		<< "X_sign should be +/-1 to set production inside/outside X"
		<< std::endl;
      exit(EXIT_FAILURE);
    }	
    // Set the variable values
    setId("Creation::SpatialCoordinate");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_max";
    tmp[1] = "X (K_Hill)";
    tmp[2] = "n_Hill";
    tmp[3] = "X_sign";  
    setParameterId( tmp );
  }

  void SpatialCoordinate::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    size_t xIndex = variableIndex(1,0);
    double powK_ = std::pow(parameter(1),parameter(2));
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      
      //Calculate cell center from vertices positions
      std::vector<double> cellCenter;
      cellCenter = T.cell(cellI).positionFromVertex(vertexData);
      assert( cellCenter.size() == vertexData[0].size() );
      
      double powX = std::pow(cellCenter[xIndex],parameter(2));
      
      if (parameter(3)>0.0)
	cellDerivs[cellI][cIndex] += parameter(0)*powX/(powK_+powX);
      else
	cellDerivs[cellI][cIndex] += parameter(0)*powK_/(powK_+powX);
    }
  }
  
  SpatialPlane::
  SpatialPlane(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=3 ) {
      std::cerr << "Creation::SpatialPlane::SpatialPlane() "
		<< "Uses three parameters V_max X(K_Hill) and X_sign."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Creation::SpatialPlane::"
		<< "SpatialPlane() "
		<< "Two levels of indices used: index for variable to be updated given as first index, "
		<< "index for the spatial coordinate to use given as the second index" << std::endl;
      exit(EXIT_FAILURE);
    }
    // Sign should be -/+1
    if( paraValue[2] != -1 && paraValue[2] != 1 ) {
      std::cerr << "Creation::SpatialCoordinate::SpatialCoordinate() "
		<< "X_sign should be +/-1 to set production inside/outside X"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
	
    // Set the variable values
    setId("Creation::SpatialCoordinate");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_max";
    tmp[1] = "X (K_Hill)";
    tmp[2] = "n_Hill";
    setParameterId( tmp );
  }

  void SpatialPlane::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    
    size_t cIndex = variableIndex(0,0);
    size_t xIndex = variableIndex(1,0);
    
    double X = parameter(1);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      //Calculate cell center from vertices positions
      std::vector<double> cellCenter;
      cellCenter = T.cell(cellI).positionFromVertex(vertexData);
      assert( cellCenter.size() == vertexData[0].size() );
      double x_coord = cellCenter[xIndex];
      double deriv_add=0;

      if (parameter(2)<0.0)
	deriv_add = ((x_coord <= X) ? parameter(0) : 0);
      else
	deriv_add = ((x_coord >= X) ? parameter(0) : 0);
      
      cellDerivs[cellI][cIndex] += deriv_add;
    }
  }

  FromList::
  FromList(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=1 && paraValue.size()!=2 ) {
      std::cerr << "Creation::FromList::FromList() "
		<< "Uses one or two parameter(s)  k_c, the constant production rate"
		<< "and a flag(equal to 1) if constant -amount- creation is needed"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() < 1 ) {
      std::cerr << "Creation::FromList::"
		<< "FromList() "
		<< "Two levels of indices used: index for variable to be updated given "
		<< "as first index, a list of cell indices given in the second level" 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    setId("Creation::FromList");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_c";
    if(numParameter()>1)
      tmp[1] = "number_flag";
    setParameterId(tmp);
    proCells=indValue[1].size();
  }
  
  void FromList::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    size_t cIndex = variableIndex(0,0);
    if(numParameter()>1 && parameter(1)==1) {
      //Add const number of molecules to the listed cells
      for (size_t cellI = 0; cellI < proCells; ++cellI) 
	cellDerivs[variableIndex(1,cellI)][cIndex] += parameter(0) /
	  T.cell(cellI).calculateVolume(vertexData);
    }
    else {
      //Add const conc to the listed cells 
      for (size_t cellI = 0; cellI < proCells; ++cellI) 
	cellDerivs[variableIndex(1,cellI)][cIndex] += parameter(0);
    }
  }
  
  OneGeometric::
  OneGeometric(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Creation::OneGeometric::"
		<< "OneGeometric() "
		<< "Uses one parameter k_c (linear production rate)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
      std::cerr << "Creation::OneGeometric::"
		<< "OneGeometric() "
		<< "Index for variable to be updated given in first row and "
		<< "index for production-dependent variable in 2nd." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Creation::OneGeometric");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "k_c";
    setParameterId( tmp );
  }

  void OneGeometric::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    size_t cIndex = variableIndex(0,0);
    size_t xIndex = variableIndex(1,0);
    double k_c = parameter(0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      double cellVolume = T.cell(cellI).calculateVolume(vertexData);	
      cellDerivs[cellI][cIndex] += cellVolume*k_c*cellData[cellI][xIndex];
    }
  }

  Sinus::Sinus(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=3 ) {
      std::cerr << "Creation::Sinus::Sinus() "
		<< "Uses three parameters amplitude, period, and phase." << std::endl;
      exit(EXIT_FAILURE);
  }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Creation::Sinus::"
		<< "Sinus() "
		<< "Index for variable to be updated given." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    setId("Creation::Sinus");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "amplitude";
    tmp[1] = "frequency";
    tmp[2] = "phase";
  }

  void Sinus::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    //Do the update for each cell
    size_t numCells = T.numCell();
    size_t cIndex = variableIndex(0,0);
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {      
      cellDerivs[cellI][cIndex] += parameter(0)*
	( 1.0 + std::sin(6.28*(time_/parameter(1) + parameter(2) ) ) );};
  }

  void Sinus::
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
      double value = k_c*
	( 1.0 + std::sin(6.28*(time_/parameter(1) + parameter(2) ) ) );
      
      cellDerivs[cellI][cIndex]  += value;
      sdydtCell[cellI][cIndex]  += value;
    }
  }

  void Sinus::
  update(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &walldata,
	 DataMatrix &vertexData,
	 double h) 
  {
    time_+=h;
  }

  void Sinus::
  initiate(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &walldata,
	   DataMatrix &vertexData,
	   DataMatrix &cellderivs, 
	   DataMatrix &wallderivs,
	   DataMatrix &vertexDerivs)
  {
    time_=0;
  }
} // end namespace Creation
