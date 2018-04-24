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
  
  Cylinder::Cylinder(std::vector<double> &paraValue,
		     std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
      std::cerr << "Force::Cylinder::"
		<< "Cylinder() "
		<< "Uses two parameters K_force direction(-1 -> inwards)"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 1.0 && paraValue[1] != -1.0) {
      std::cerr << "Force::Cylinder::"
		<< "Cylinder() "
		<< "direction (second parameter) needs to be 1 (outward) "
		<< "or -1 (inwards)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
      std::cerr << "Force::Cylinder::"
		<< "Cylinder() "
		<< "No indices used." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::Cylinder");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    tmp[1] = "direction";
    setParameterId(tmp);
  }

  void Cylinder::
  derivs(Tissue &T, DataMatrix &cellData,
	 DataMatrix &wallData, DataMatrix &vertexData,
	 DataMatrix &cellDerivs, DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs) {
    double coeff = parameter(0) * parameter(1);
    size_t dimension = vertexData[0].size();
    size_t lastPosIndex = dimension - 1;
    // For each vertex
    for (size_t i = 0; i < T.numVertex(); ++i) {
      // On cylinder
      double norm = 0.0;
      for (size_t d = 0; d < lastPosIndex; d++)
	norm += vertexData[i][d] * vertexData[i][d];
      if (norm > 0.0)
	norm = 1.0 / std::sqrt(norm);
      else
	norm = 0.0;
      for (size_t d = 0; d < lastPosIndex; d++)
	vertexDerivs[i][d] += coeff * norm * vertexData[i][d];
    }
  }
  
  SphereCylinder::SphereCylinder(
				 std::vector<double> &paraValue,
				 std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
      std::cerr << "Force::SphereCylinder::"
		<< "SphereCylinder() "
		<< "Uses two parameters K_force direction (-1 -> inwards)"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 1.0 && paraValue[1] != -1.0) {
      std::cerr << "Force::SphereCylinder::"
		<< "SphereCylinder() "
		<< "direction (second parameter) needs to be 1 (outward) "
		<< "or -1 (inwards)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
      std::cerr << "Force::SphereCylinder::"
		<< "SphereCylinder() "
		<< "Uses no indices." << std::endl;
      exit(EXIT_FAILURE);
  }
    // Set the variable values
    //
    setId("Force::SphereCylinder");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    tmp[1] = "direction";
    setParameterId(tmp);
  }
  
  void SphereCylinder::derivs(Tissue &T, DataMatrix &cellData,
			      DataMatrix &wallData, DataMatrix &vertexData,
			      DataMatrix &cellDerivs, DataMatrix &wallDerivs,
			      DataMatrix &vertexDerivs) {
    double coeff = parameter(0) * parameter(1);
    size_t dimension = vertexData[0].size();
    size_t lastPosIndex = dimension - 1;
    // For each vertex
    for (size_t i = 0; i < T.numVertex(); ++i) {
      if (vertexData[i][lastPosIndex] > 0.0) {
	// On sphere
	double norm = 0.0;
	for (size_t d = 0; d < dimension; d++)
	  norm += vertexData[i][d] * vertexData[i][d];
	if (norm > 0.0)
	  norm = 1.0 / std::sqrt(norm);
	else
	  norm = 0.0;
	for (size_t d = 0; d < dimension; d++)
	  vertexDerivs[i][d] += coeff * norm * vertexData[i][d];
      } else {
	// On cylinder
	double norm = 0.0;
	for (size_t d = 0; d < dimension - 1; d++)
	  norm += vertexData[i][d] * vertexData[i][d];
	if (norm > 0.0)
	  norm = 1.0 / std::sqrt(norm);
	else
	  norm = 0.0;
	for (size_t d = 0; d < dimension - 1; d++)
	  vertexDerivs[i][d] += coeff * norm * vertexData[i][d];
      }
    }
  }
  
  SphereCylinderRadius::
  SphereCylinderRadius(
		       std::vector<double> &paraValue,
		       std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 3) {
      std::cerr << "Force::SphereCylinderRadius::"
		<< "SphereCylinderRadius() "
		<< "Uses three parameters F_out, F_in and radius." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[0] < 0.0 || paraValue[1] < 0.0 || paraValue[2] < 0.0) {
      std::cerr << "Force::SphereCylinderRadius::"
		<< "SphereCylinderRadius() "
		<< "Parameters need to be positive." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
      std::cerr << "Force::SphereCylinderRadius::"
		<< "SphereCylinderRadius() "
		<< "No indices used." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::SphereCylinderRadius");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "F_out";
    tmp[1] = "F_in";
    tmp[2] = "R";
    setParameterId(tmp);
  }

  void SphereCylinderRadius::
  derivs(Tissue &T, DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    size_t lastPosIndex = dimension - 1;
    double F_out = parameter(0);
    double F_in = parameter(1);
    double R = parameter(2);
    // For each vertex
    for (size_t i = 0; i < T.numVertex(); ++i) {
      if (vertexData[i][lastPosIndex] > 0.0) {
	// On sphere
	double norm = 0.0;
	for (size_t d = 0; d < dimension; d++)
	  norm += vertexData[i][d] * vertexData[i][d];
	norm = std::sqrt(norm);
	if (norm < R) {
	  double coeff = F_out * (R - norm) / norm;
	  for (size_t d = 0; d < dimension; d++)
	    vertexDerivs[i][d] += coeff * vertexData[i][d];
	} else {
	  double coeff = F_in * (R - norm) / norm;
	  for (size_t d = 0; d < dimension; d++)
	    vertexDerivs[i][d] += coeff * vertexData[i][d];
	}
      } else {
	// On cylinder
	double norm = 0.0;
	for (size_t d = 0; d < dimension - 1; d++)
	  norm += vertexData[i][d] * vertexData[i][d];
	norm = std::sqrt(norm);
	if (norm < R) {
	  double coeff = F_out * (R - norm) / norm;
	  for (size_t d = 0; d < dimension - 1; d++)
	    vertexDerivs[i][d] += coeff * vertexData[i][d];
	} else {
	  double coeff = F_in * (R - norm) / norm;
	  for (size_t d = 0; d < dimension - 1; d++)
	    vertexDerivs[i][d] += coeff * vertexData[i][d];
	}
      }
    }
  }
  
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
  
  
  EpidermalRadial::
  EpidermalRadial(
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
      double R = std::sqrt(x * x + y * y);
      if (R == 0) continue;
      x /= R;
      y /= R;
      
      // std::cerr << "Vertex " << vertex.index() << std::endl;
      // std::cerr << " x = " << x << std::endl;
      // std::cerr << " y = " << y << std::endl;
      
      vertexDerivs[vertex.index()][0] += parameter(0) * x;
      vertexDerivs[vertex.index()][1] += parameter(0) * y;
    }
  }

  IndexRadial::IndexRadial(std::vector<double> &paraValue,
			   std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
      std::cerr << "Force::IndexRadial::"
		<< "IndexRadial() "
		<< "Uses two parameters K_force direction_flag(-1 -> inwards)"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 1.0 && paraValue[1] != -1.0) {
      std::cerr << "Force::IndexRadial::"
		<< "IndexRadial() "
		<< "direction_flag (second parameter) needs to be 1 (outward) "
		<< "or -1 (inwards)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() < 1) {
      std::cerr << "Force::IndexRadial::"
		<< "IndexRadial() "
		<< "Vertex indices to be updated in first level." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::IndexRadial");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    tmp[1] = "direction_flag";
    setParameterId(tmp);
  }
  
  void IndexRadial::derivs(Tissue &T, DataMatrix &cellData,
			   DataMatrix &wallData,
			   DataMatrix &vertexData,
			   DataMatrix &cellDerivs,
			   DataMatrix &wallDerivs,
			   DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[variableIndex(0, 0)].size();
    double coeff = parameter(0) * parameter(1);
    // For each vertex in list
    for (size_t k = 0; k < numVariableIndex(0); ++k) {
      size_t i = variableIndex(0, k);
      for (size_t d = 0; d < dimension; d++)
	vertexDerivs[i][d] += coeff * vertexData[i][d];
    }
  }
  
  CellIndexRadial::CellIndexRadial(std::vector<double> &paraValue,
				   std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
      std::cerr << "Force::CellIndexRadial::"
		<< "CellIndexRadial() "
		<< "Uses two parameters K_force direction_flag (1 -> outwards, "
		<< "-1 -> inwards)"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 1.0 && paraValue[1] != -1.0) {
      std::cerr << "Force::CellIndexRadial::"
		<< "CellIndexRadial() "
		<< "direction_flag (second parameter) needs to be 1 (outward) "
		<< "or -1 (inwards)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() < 1) {
      std::cerr << "Force::CellIndexRadial::"
		<< "CellIndexRadial() "
		<< "Cell indices for which vertices are to be updated to be "
		<< "provided in first level." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::CellIndexRadial");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    tmp[1] = "direction_flag";
    setParameterId(tmp);
  }
  
  void CellIndexRadial::derivs(Tissue &T, DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs) {
    double coeff = parameter(0) * parameter(1);
    size_t dimension = vertexData[0].size();
    // For each vertex in cell list
    for (size_t k = 0; k < numVariableIndex(0); ++k) {
      size_t cellI = variableIndex(0, k);
      for (size_t l = 0; l < T.cell(cellI).numVertex(); ++l) {
	size_t i = T.cell(cellI).vertex(l)->index();
	for (size_t d = 0; d < dimension; d++)
	  vertexDerivs[i][d] += coeff * vertexData[i][d];
      }
    }
  }
  
  Axial::Axial(
	       std::vector<double> &paraValue,
	       std::vector<std::vector<size_t>> &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 4) {
      std::cerr << "Force::Axial::"
		<< "Axial()"
		<< "Applies forces F axially (in z) to the vertices in the regions "
		<< "between [z0+a,z0+a+d] and [z0-a,z0-a-d]"
		<< "uses 4 parameters: z0, a, d and F" << std::endl
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
      std::cerr << "Force::Axial::"
		<< "Axial() "
		<< "No variable indices used." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Force::Axial");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "z_0";  // the center of the growing region
    tmp[1] = "a";   // growing region threshold interval
    tmp[2] = "d";   // width of force application interval
    tmp[3] = "F";   // Force for each node
    
    setParameterId(tmp);
  }
  
  void Axial::derivs(Tissue &T, DataMatrix &cellData,
		     DataMatrix &wallData,
		     DataMatrix &vertexData,
		     DataMatrix &cellDerivs,
		     DataMatrix &wallDerivs,
		     DataMatrix &vertexDerivs) {
    // Do the update for each vertex .
    size_t numVertex = T.numVertex();
    for (size_t vertexIndex = 0; vertexIndex < numVertex; ++vertexIndex) {
      double Y0 = parameter(0);
      double a = parameter(1);
      double d = parameter(2);
      double F = parameter(3);
      
      DataMatrix position(1, vertexData[vertexIndex]);
      
      if (position[0][2] > Y0 + a && position[0][2] < Y0 + a + d) {
	// vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
	// vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
	vertexDerivs[vertexIndex][2] += F;
      }
      if (position[0][2] < Y0 - a && position[0][2] > Y0 - a - d) {
	// vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
	// vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
	vertexDerivs[vertexIndex][2] -= F;
      }
    }
  }
  
  // HJ: Old(?) version of this reaction (before renamed Force::Axial), more advanced given that it provides
  // axial and radial forces (using 8 parameters). It might be that one of these versions is used for the future PLoS Comp Biol
  // and one is for the (3D cell) hypocotyl work with Siobhan... 
  // VertexFromHypocotylGrowth::
  // VertexFromHypocotylGrowth(std::vector<double> &paraValue,
  // 		       std::vector< std::vector<size_t> >
  // 		       &indValue )
  // {
  //   //Do some checks on the parameters and variable indeces
  //   //
  //   if( paraValue.size()!=8 ) {
  //     std::cerr << "VertexFromHypocotylGrowth::"
  // 	      << "VertexFromHypocotylGrowth()"
  // 	      << "Applies forces by simple springs from a growing center in the
  // Hypocotyl"
  //               << "uses 7 parameters: Y0, a, d, delta, lambda, b, K, epsilon"
  //               << std::endl
  //               << std::endl;
  //     exit(0);
  //   }
  //   if( indValue.size() != 0 ) {
  //     std::cerr << "VertexFromHypocotylGrowth::"
  // 	      << "VertexFromHypocotylGrowth() "
  // 	      << "No variable indices used." << std::endl;
  //     exit(0);
  //   }
  //   //Set the variable values
  //   //
  //   setId("VertexFromHypocotylGrowth");
  //   setParameter(paraValue);
  //   setVariableIndex(indValue);
  
  //   //Set the parameter identities
  //   //
  //   std::vector<std::string> tmp( numParameter() );
  //   tmp[0] = "y0";      // the center of the growing region
  //   tmp[1] = "a";       // growing region interwall
  //   tmp[2] = "d";       // distance between the boundary of the growing region
  //   and epidermis tmp[3] = "delta";   // maximum displacement ( at the center
  //   of the growing region) tmp[4] = "lambda";  // Gaussian distribution of
  //   displacement of the growing center tmp[5] = "b";       // interwall on
  //   epidermis recieving the force tmp[6] = "Ks";       // spring constant
  //   tmp[7] = "epsilon"; // deltaX for evaluating the integral force
  
  //   setParameterId( tmp );
  // }
  
  // void VertexFromHypocotylGrowth::
  // derivs(Tissue &T,
  //        DataMatrix &cellData,
  //        DataMatrix &wallData,
  //        DataMatrix &vertexData,
  //        DataMatrix &cellDerivs,
  //        DataMatrix &wallDerivs,
  //        DataMatrix &vertexDerivs ) {
  
  //   //Do the update for each vertex .
  //   size_t numVertex = T.numVertex();
  //   for (size_t vertexIndex=0 ; vertexIndex<numVertex; ++vertexIndex) {
  //     double Y0=parameter(0);
  //     double a=parameter(1);
  //     double d=parameter(2);
  //     double delta=parameter(3);
  //     double lambda=parameter(4);
  //     double b=parameter(5);
  //     double Ks=parameter(6);
  //     double epsilon=parameter(7);
  
  //     DataMatrix position(1,vertexData[vertexIndex]);
  
  //     if( position[0][2] > Y0-b && position[0][2] < Y0+b ){
  
  //       double Faxial=0;
  //       double Fradial=0;
  //       double X=Y0-a;
  //       double Y=position[0][2];
  //       while (X < Y0+a){
  // 	int sgnx=0;
  //         if ( X-Y0 > 0 ) sgnx=1;
  // 	else sgnx=-1;
  
  // 	Faxial -=(epsilon/(2*a))*Ks
  // 	  *(    (  std::sqrt( (Y-X)*(Y-X)+d*d )
  // 		   /std::sqrt((Y-X-sgnx*delta*exp(-lambda*X*X))*(Y-X-sgnx*delta*exp(-lambda*X*X))+d*d
  // )
  // 		   )
  // 		-1 )
  //           *(Y-X-sgnx*delta*exp(-lambda*X*X));
  
  // 	Fradial =(epsilon/(2*a))*Ks
  // 	  *(    (  std::sqrt( (Y-X)*(Y-X)+d*d )
  // 		   /std::sqrt((Y-X-sgnx*delta*exp(-lambda*X*X))*(Y-X-sgnx*delta*exp(-lambda*X*X))+d*d
  // )
  // 		   )
  // 		-1 )
  //           *(d);
  // 	  X+=epsilon;
  //       }
  //       //vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
  //       //vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
  //       vertexDerivs[vertexIndex][2]+=Faxial;
  //     }
  //   }
  // }
  
} //end namespace Force
