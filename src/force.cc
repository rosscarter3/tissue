//
// Filename     : force.cc
// Description  : Classes describing updates due to external forces compiled in namespace Force
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : April 2018
// Revision     : $Id:$
//
#include "force.h"
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "tissue.h"


namespace Force {
  InfiniteWall::InfiniteWall(
			     std::vector<double> &paraValue,
			     std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 3) {
      std::cerr << "Force::InfiniteWall::"
		<< "InfiniteWall() "
		<< "Uses three parameters k_spring threshold and direction "
		<< "(-1 -> less than)" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[2] != 1.0 && paraValue[2] != -1.0) {
      std::cerr << "Force::InfiniteWall::"
		<< "InfiniteWall() "
		<< "direction (third parameter) need to be 1 (greater than) "
		<< "or -1 (less than)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
      std::cerr << "Force::InfiniteWall::"
		<< "InfiniteWall() "
		<< "Pos (coordinate) index in first level." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::InfiniteWall");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_spring";
    tmp[1] = "threshold";
    tmp[2] = "direction";
    setParameterId(tmp);
  }
  
  void InfiniteWall::derivs(Tissue &T, DataMatrix &cellData,
			    DataMatrix &wallData, DataMatrix &vertexData,
			    DataMatrix &cellDerivs, DataMatrix &wallDerivs,
			    DataMatrix &vertexDerivs) {
    // Check the cancelation for every vertex
    size_t numVertices = T.numVertex();
    size_t posIndex = variableIndex(0, 0);
    // size_t dimension = vertexData[0].size();
    assert(posIndex < vertexData[0].size());
    
    for (size_t i = 0; i < numVertices; ++i) {
      if (parameter(2) > 0 && vertexData[i][posIndex] > parameter(1)) {
	vertexDerivs[i][posIndex] -=
	  parameter(0) * (vertexData[i][posIndex] - parameter(1));
      } else if (parameter(2) < 0 && vertexData[i][posIndex] < parameter(1)) {
	vertexDerivs[i][posIndex] -=
	  parameter(0) * (vertexData[i][posIndex] - parameter(1));
      }
    }
  }
  
  EpidermalCoordinate::EpidermalCoordinate(
					   std::vector<double> &paraValue,
					   std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
      std::cerr << "Force::EpidermalCoordinate::"
		<< "EpidermalCoordinate() "
		<< "Uses two parameters k_strength direction_flag "
		<< "(-1 -> down)" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 1.0 && paraValue[1] != -1.0) {
      std::cerr << "Force::EpidermalCoordinate::"
		<< "EpidermalCoordinate() "
		<< "direction (second parameter) need to be 1 (pos) "
		<< "or -1 (neg direction)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
      std::cerr << "Force::EpidermalCoordinate::"
		<< "EpidermalCoordinate() "
		<< "Pos (coordinate) index in first level." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::EpidermalCoordinate");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_spring";
    tmp[1] = "direction";
    setParameterId(tmp);
  }
  
  void EpidermalCoordinate::derivs(Tissue &T, DataMatrix &cellData,
			       DataMatrix &wallData, DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs) {
    // Check the cancelation for every vertex
    size_t numVertices = T.numVertex();
    size_t posIndex = variableIndex(0, 0);
    // size_t dimension = vertexData[0].size();
    assert(posIndex < vertexData[0].size());
    
    for (size_t i = 0; i < numVertices; ++i) {
      int epidermisFlag = 0;
      for (size_t k = 0; k < T.vertex(i).numWall(); ++k) {
	if (T.vertex(i).wall(k)->cell1() == T.background() ||
	    T.vertex(i).wall(k)->cell1() == T.background()) {
	  epidermisFlag++;
	}
      }
      if (epidermisFlag) {
	vertexDerivs[i][posIndex] += parameter(0) * parameter(1);
      }
    }
  }
  
  EpidermalRadial::EpidermalRadial(
				   std::vector<double> &paraValue,
				   std::vector<std::vector<size_t>> &indValue) {
    if (paraValue.size() != 1) {
      std::cerr << "Force::EpidermalRadial::EpidermalRadial() "
		<< "Uses one parameter: Force" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (indValue.size() != 0) {
      std::cerr << "Force::EpidermalRadial::EpidermalRadial() "
		<< "No indices are given." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    setId("Force::EpidermalRadial");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "Force";
    
    setParameterId(tmp);
  }
  
  void EpidermalRadial::derivs(Tissue &T, DataMatrix &cellData,
			       DataMatrix &wallData, DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs) {
    for (size_t n = 0; n < T.numVertex(); ++n) {
      Vertex vertex = T.vertex(n);
      
      // If vertex is not in the epidermal layer then skip.
      bool isEpidermalVertex = false;
      for (size_t i = 0; i < vertex.numWall(); ++i) {
	Wall *wall = vertex.wall(i);
	if ((wall->cell1()->index() == (size_t)-1) ||
	    (wall->cell2()->index() == (size_t)-1)) {
	  isEpidermalVertex = true;
	  break;
	}
      }
      if (isEpidermalVertex == false) {
	continue;
      }
      
      double x = vertex.position(0);
      double y = vertex.position(1);
      double A = std::sqrt(x * x + y * y);
      if (A == 0) continue;
      x /= A;
      y /= A;
      
      // 		std::cerr << "Vertex " << vertex.index() << std::endl;
      // 		std::cerr << " x = " << x << std::endl;
      // 		std::cerr << " y = " << y << std::endl;
      
      vertexDerivs[vertex.index()][0] += -parameter(0) * x;
      vertexDerivs[vertex.index()][1] += -parameter(0) * y;
    }
  }
} //end namespace Force
