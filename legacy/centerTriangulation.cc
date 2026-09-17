//
// Filename     : centerTriangulation.cc
// Description  : Classes describing updates for tissues with cells storing a central point (and internal edges)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2012
// Revision     : $Id:$
//
#include <vector>
#include "baseReaction.h"
#include "centerTriangulation.h"
#include "tissue.h"

#include <ctime>
using namespace std;

namespace CenterTriangulation {

  Initiate::
  Initiate(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > &indValue)
  {
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=0 && paraValue.size()!=1 && paraValue.size()!=2 ) {
      std::cerr << "CenterTriangulation::Initiate::Initiate() "
		<< "Uses zero or one or two parameter(s). " << std::endl
		<< "If zero parameters given it will initiate a single edge CenterTriangulation "
		<< "from scratch if the initial tissue has no CenterTriangulation and "
		<< "using the information from the initial tissue if CT exist." << std::endl
		<< "If one parameter is provided it is a flag equal to one for overriding existing "
		<< "CenterTriangulation information (and initiate CT from scratch always)." << std::endl
                << "If two parameters are provided the second should be equal to zero for initiation"
		<< "of single internal edges (as default); and one for initiation of two "
                << "independent edges in the CT (e.g. used for growth dependent reactions). It uses"
		<< "a central triangulation stored in tissue if available, or will create from scratch if"
		<< " the initial tissue does not have a CT or if p[0]=1 overiding current information."
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // Warn for obselete use of p[1]=2 (now replaced with p[0]=0,p[1]=1)
    if (paraValue.size()==2 && paraValue[1]==2) {
      std::cerr << "CenterTriangulation::Initiate::Initiate() "
		<< "The use of p[1]=2 no longer needed (or allowed). "
		<< "Use combination p[0]=0 (no overriding) and p[1]=1 (double edges) "
		<< "[and provide a CT init file] for using single edge CenterTriangulation "
		<< "input information to "
		<< "create a double-edge CenterTriangulation." << std::endl;
      exit(EXIT_FAILURE);
    }
      
    if( indValue.size()!=1 || indValue[0].size()!=1 ) { 
      std::cerr << "CenterTriangulation::Initiate::Initiate() "
		<< "Start of internal cell variables index given in first level "
		<< "(=cell.numVariable)."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Set the variable values
    setId("CenterTriangulation::Initiate");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    if (paraValue.size())
      tmp[0] = "overRide_flag";
    if (paraValue.size()==2)
      tmp[1] = "doubleEdge_flag";
  
    setParameterId( tmp );    
  }
  
  void Initiate::
  initiate(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs)
  {
    size_t dimension=3; //Only implemented for 3D models
    assert (dimension==vertexData[0].size());
    size_t numVariable = T.cell(0).numVariable();
    assert (numVariable==cellData[0].size());

    // Create the new variables
    if (variableIndex(0,0) != numVariable) {
      std::cerr << "CenterTriangulation::Initiate::initiate() "
		<< "Wrong index given as start index for additional variables."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    size_t numCell = cellData.size();
    std::vector<double> com(dimension);    
    assert (numCell==T.numCell());
    //
    // CenterTriangulation created from scratch if no CT provided in init file or override flag set
    //
    if (!T.cell(0).numCenterPosition() || (numParameter() && parameter(0)==1) ) {
      //
      // Double edge CT created from scratch if doubleEdge_flag set
      //
      if (numParameter()==2 && parameter(1)==1) {
	for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {	
          size_t numInternalWall = T.cell(cellIndex).numVertex();
          cellData[cellIndex].resize(numVariable+dimension+3*numInternalWall);
          cellDerivs[cellIndex].resize(numVariable+dimension+3*numInternalWall);
          com = T.cell(cellIndex).positionFromVertex(vertexData);
          // Set center position to com of the cell
          for (size_t d=0; d<dimension; ++d)
            cellData[cellIndex][numVariable+d] = com[d];    
          // Set internal wall lengths to the distance btw com and the vertex
          for (size_t k=0; k<numInternalWall; ++k) {// internal walls
            Vertex *tmpVertex = T.cell(cellIndex).vertex(k); 
            size_t vertexIndex = tmpVertex->index();
            double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
                                         (com[0]-vertexData[vertexIndex][0])+
                                         (com[1]-vertexData[vertexIndex][1])*
                                         (com[1]-vertexData[vertexIndex][1])+
                                         (com[2]-vertexData[vertexIndex][2])*
                                         (com[2]-vertexData[vertexIndex][2]) );   
            cellData[cellIndex][numVariable+dimension+2*k] = distance;
            cellData[cellIndex][numVariable+dimension+2*k+1] = distance;
          }
          for (size_t k=0; k<numInternalWall; ++k) { // growing component of main walls zero at initiation
            size_t numInternalWall = T.cell(cellIndex).numVertex();
            cellData[cellIndex][numVariable+dimension+2*numInternalWall+k] =0;
          }	  
	}
	std::cerr << "CenterTriangulation::Initiate::initiate() "
		  << "Initiating double edge CT from scratch." << std::endl; 
      }
      //
      // Single edge CT created from scratch
      //
      else {
	for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {	
	  size_t numInternalWall = T.cell(cellIndex).numVertex();
	  cellData[cellIndex].resize(numVariable+dimension+numInternalWall);
	  cellDerivs[cellIndex].resize(numVariable+dimension+numInternalWall);
	  com = T.cell(cellIndex).positionFromVertex(vertexData);
	  // Set center position to com of the cell
	  for (size_t d=0; d<dimension; ++d)
	    cellData[cellIndex][numVariable+d] = com[d];    
	  // Set internal wall lengths to the distance btw com and the vertex
	  for (size_t k=0; k<numInternalWall; ++k) {
            Vertex *tmpVertex = T.cell(cellIndex).vertex(k); 
            size_t vertexIndex = tmpVertex->index();
            double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
                                         (com[0]-vertexData[vertexIndex][0])+
                                         (com[1]-vertexData[vertexIndex][1])*
                                         (com[1]-vertexData[vertexIndex][1])+
                                         (com[2]-vertexData[vertexIndex][2])*
                                         (com[2]-vertexData[vertexIndex][2]) );   
            cellData[cellIndex][numVariable+dimension+k] = distance;            
          }
        }
	std::cerr << "CenterTriangulation::Initiate::initiate() "
		  << "Initiating single edge CT from scratch." << std::endl; 
      }
    }
    //
    // Else use information from CT in init file to initiate CT
    //
    else if (numParameter()==2 && parameter(1)==1) { // Double edge from init file
      // Copy central position and edge length (twice) data from cells in Tissue
      for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
	// Do the check that CT exist in init (should not be needed)
	if (!T.cell(cellIndex).numCenterPosition()) {
	  std::cerr << "CenterTriangulation::Initiate::initiate() ERROR: "
		    << "Original tissue is expected "
		    << "to be of the form CenterTriangulation (with single edge elements). Did you "
		    << "forget to provide -centerTri_init to simulator?" << std::endl;
	  std::cerr << "Error triggered for cell " << cellIndex << std::endl;
	  std::exit(EXIT_FAILURE);
	}
        size_t numInternalWall = T.cell(cellIndex).numVertex();
        com = T.cell(cellIndex).centerPosition();
        
        cellData[cellIndex].resize(numVariable+dimension+3*numInternalWall);
        cellDerivs[cellIndex].resize(numVariable+dimension+3*numInternalWall);
	
        for (size_t k=0; k<dimension; ++k) {// center position
          cellData[cellIndex][numVariable+k] = com[k];
        }
	
        for (size_t k=0; k<numInternalWall; ++k) {// internal walls
          cellData[cellIndex][numVariable+dimension+2*k] = T.cell(cellIndex).edgeLength(k);
          cellData[cellIndex][numVariable+dimension+2*k+1] =T.cell(cellIndex).edgeLength(k); 
        }
        for (size_t k=0; k<numInternalWall; ++k) { // growing component of main walls zero at initiation
          cellData[cellIndex][numVariable+dimension+2*numInternalWall+k] =0;
        }
      }
      std::cerr << "CenterTriangulation::Initiate::initiate() "
		<< "Initiating double edge CT from init file." << std::endl; 
    }
    else { // Single edge from init file      
      // Copy central position and edge length data from cells in Tissue
      for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
	// Do the check that CT exist in init (should not be needed)
	if (!T.cell(cellIndex).numCenterPosition()) {
	  std::cerr << "CenterTriangulation::Initiate::initiate() ERROR: "
		    << "Original tissue is expected "
		    << "to be of the form CenterTriangulation (with single edge elements). Did you "
		    << "forget to provide -centerTri_init to simulator?" << std::endl;
	  std::cerr << "Error triggered for cell " << cellIndex << std::endl;
	  std::exit(EXIT_FAILURE);
	}
	size_t numInternalWall = T.cell(cellIndex).numVertex();
	cellData[cellIndex].resize(numVariable+dimension+numInternalWall);
	cellDerivs[cellIndex].resize(numVariable+dimension+numInternalWall);
	com = T.cell(cellIndex).centerPosition();
	
	// Set center position to com of the cell
	for (size_t d=0; d<dimension; ++d)
	  cellData[cellIndex][numVariable+d] = com[d];    
	// Set internal wall lengths to the distance btw com and the vertex
	for (size_t k=0; k<numInternalWall; ++k) {
	  cellData[cellIndex][numVariable+dimension+k] = T.cell(cellIndex).edgeLength(k);
	}  
      }
      std::cerr << "CenterTriangulation::Initiate::initiate() "
		<< "Initiating single edge CT from init file." << std::endl; 
    }
  }
  
  void Initiate::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    // Does nothing
  }
  
  
}
