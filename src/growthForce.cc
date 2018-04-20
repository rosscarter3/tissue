//
// Filename     : growthForce.cc
// Description  : Classes describing growth updates by applying forces on vertices
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : April 2018
// Revision     : $Id:$
//

#include "growthForce.h"
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "tissue.h"

namespace GrowthForce {
  Radial::
  Radial(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
      std::cerr << "GrowthForce::Radial::"
		<< "Radial() "
		<< "Uses two parameters k_growth and r_pow (0,1)" << std::endl;
      exit(EXIT_FAILURE);
    }  
    if( indValue.size() != 0 ) {
      std::cerr << "GrowthForce::Radial::"
		<< "Radial() "
		<< "No variable index is used." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::Radial");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[0] = "r_pow";
    setParameterId( tmp );
  }

  void Radial::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numVertices = T.numVertex();
    size_t dimension=vertexData[0].size();
    
    for( size_t i=0 ; i<numVertices ; ++i ) {
      double fac=parameter(0);
      if( parameter(1)==0.0 ) {
	double r=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  r += vertexData[i][d]*vertexData[i][d];
	if( r>0.0 )
	  r = std::sqrt(r);
	if( r>0.0 )
	  fac /= r;
	else
	  fac=0.0;
      }
      for( size_t d=0 ; d<dimension ; ++d )
	vertexDerivs[i][d] += fac*vertexData[i][d];
    }
  }
  
  void Radial::derivsWithAbs(Tissue &T,
			     DataMatrix &cellData,
			     DataMatrix &wallData,
			     DataMatrix &vertexData,
			     DataMatrix &cellDerivs,
			     DataMatrix &wallDerivs,
			     DataMatrix &vertexDerivs,
			     DataMatrix &sdydtCell,
			     DataMatrix &sdydtWall,
			     DataMatrix &sdydtVertex )  {
    return derivs(T,cellData,wallData,vertexData,cellDerivs,wallDerivs,vertexDerivs);
  }

  namespace CenterTriangulation {
    Radial::
    Radial(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {      
      // Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
	std::cerr << "GrowthForce::CenterTriangulation::Radial::"
		  << "Radial() " << std::endl
		  << "Uses two parameters k_growth and r_pow (0,1)" << std::endl;
	exit(EXIT_FAILURE);
      }  
      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
	std::cerr << "GrowthForce::CenterTriangulation::Radial::"
		  << "Radial() " << std::endl
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in first level. "
		  << "See Documentation for namespace CenterTriangulation."
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      // Set the variable values
      //
      setId("GrowthForce::CenterTriangulation::Radial");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      // Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[0] = "r_pow";
      setParameterId( tmp );
    }

    void Radial::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numVertices = T.numVertex();
      size_t numCells = T.numCell();
      size_t dimension=vertexData[0].size();
      
      // Move vertices / apply force to vertices
      for( size_t i=0 ; i<numVertices ; ++i ) {
	double fac=parameter(0);
	if( parameter(1)==0.0 ) {
	  double r=0.0;
	  for( size_t d=0 ; d<dimension ; ++d )
	    r += vertexData[i][d]*vertexData[i][d];
	  if( r>0.0 )
	    r = std::sqrt(r);
	  if( r>0.0 )
	    fac /= r;
	  else
	    fac=0.0;
	}
	for( size_t d=0 ; d<dimension ; ++d )
	  vertexDerivs[i][d] += fac*vertexData[i][d];
      }
      // Move / apply force to vertices defined in cell centers 
      for( size_t i=0 ; i<numCells ; ++i ) {
	double fac=parameter(0);
	if( parameter(1)==0.0 ) {
	  double r=0.0;
	  for( size_t d=variableIndex(0,0) ; d<variableIndex(0,0)+dimension ; ++d )
	    r += cellData[i][d]*cellData[i][d];
	  if( r>0.0 )
	    r = std::sqrt(r);
	  if( r>0.0 )
	    fac /= r;
	  else
	    fac=0.0;
	}
	for( size_t d=variableIndex(0,0) ; d<variableIndex(0,0)+dimension ; ++d )
	  cellDerivs[i][d] += fac*cellData[i][d];
      }
    }
  } // end namespace CenterTriangulation
} // end namespace GrowthForce
