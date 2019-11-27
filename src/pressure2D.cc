//
// Filename     : pressure2D.cc
// Description  : Classes describing forces acting on 2D surfaces
//              : (perpendicular to edges)
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : November 2019
// Revision     : $Id:$
//
#include "pressure2D.h"
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "tissue.h"

namespace Pressure2D {

AreaPotential::AreaPotential(std::vector<double> &paraValue,
                     std::vector<std::vector<size_t>> &indValue) {
  // Do some checks on the parameters and variable indices
  //
  if (paraValue.size() != 2 || (paraValue[1] != 0.0 && paraValue[1] != 1.0)) {
    std::cerr
        << "Pressure2D::AreaPotential::"
        << "AreaPotential() "
        << "Uses two parameters P_force and normalizeVolumeFlag (= 0 or 1)."
        << std::endl;
    exit(EXIT_FAILURE);
  }

  if (indValue.size() != 0) {
    std::cerr << "Pressure2D::AreaPotential::"
              << "AreaPotential() "
              << "No index should be given." << std::endl;
    exit(EXIT_FAILURE);
  }


  // Set the variable values
  //
  setId("Pressure2D::AreaPotential");
  setParameter(paraValue);
  setVariableIndex(indValue);

  // Set the parameter identities
  //
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "P_force";
  tmp[1] = "f_V_norm";
  setParameterId(tmp);
}

void AreaPotential::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
			   DataMatrix &vertexData, DataMatrix &cellDerivs,
			   DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
  // NOTE: Assuming cells and vertices are sorted, and that we are working in
  // 2 dimensions.

  // Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition();
  bool normalize = parameter(1) == 1;
  assert(dimension == 2);

  for (size_t cellI = 0; cellI < numCells; ++cellI) {
    Cell &tmpCell = T.cell(cellI);

    double factor = 0.5 * parameter(0);
    if (normalize) {
      double cellVolume = tmpCell.calculateVolume(vertexData);
      factor /= std::fabs(cellVolume);
    }

    for (size_t k = 0; k < tmpCell.numVertex(); ++k) {
      size_t v1I = tmpCell.vertex(k)->index();
      size_t v1PlusI = tmpCell.vertex((k + 1) % (tmpCell.numVertex()))->index();
      size_t v1MinusK = k > 0 ? k - 1 : tmpCell.numVertex() - 1;
      size_t v1MinusI = tmpCell.vertex(v1MinusK)->index();

      vertexDerivs[v1I][0] +=
          factor * (vertexData[v1PlusI][1] - vertexData[v1MinusI][1]);
      vertexDerivs[v1I][1] +=
          factor * (vertexData[v1MinusI][0] - vertexData[v1PlusI][0]);
    }
  }
}

AreaPotentialTri::AreaPotentialTri(std::vector<double> &paraValue,
				   std::vector<std::vector<size_t>> &indValue) {
  // Do some checks on the parameters and variable indeces
  //
  if (paraValue.size() != 2 && paraValue.size() != 3) {
    std::cerr << "Pressure2D::AreaPotentialTri::"
              << "AreaPotentialTri() "
              << "Uses two or three parameters p0=P_force, "
              << "p1=flag_Vnorm(=1/0), [p2=flag_internalCellsOnly(=1/0)]."
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1] != 0 && paraValue[1] != 1) {
    std::cerr << "Pressure2D::AreaPotentialTri::"
              << "AreaPotentialTri() "
              << "p1=flag_Vnorm(=1/0), is a flag and needs to be set to 0 or 1."
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue.size() > 2 && (paraValue[2] != 0 && paraValue[2] != 1)) {
    std::cerr << "Pressure2D::AreaPotentialTri::"
              << "AreaPotentialTri() "
              << "p2=flag_InternalCellsOnly(=1/0), is a flag and needs to be "
                 "set to 0 or 1."
              << std::endl
              << "This parameter is optional." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 0) {
    std::cerr << "Pressure2D::AreaPotentialTri::"
              << "AreaPotentialTri() "
              << "No index should be given." << std::endl;
    exit(EXIT_FAILURE);
  }
  // Set the variable values
  //
  setId("Pressure2D::AreaPotentialTri");
  setParameter(paraValue);
  setVariableIndex(indValue);

  // Set the parameter identities
  //
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "P_force";
  tmp[1] = "f_V_norm";
  if (numParameter() > 2) {
    tmp[2] = "f_internal";
  }
  setParameterId(tmp);
}

void AreaPotentialTri::derivs(Tissue &T, DataMatrix &cellData,
                           DataMatrix &wallData, DataMatrix &vertexData,
                           DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                           DataMatrix &vertexDerivs) {
  // Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition();

  // For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {
    Cell &tmpCell = T.cell(cellI);
    // If p2 flag is set, this will only be done for internal cells
    if (numParameter() < 3 || parameter(2) != 1 ||
        !(tmpCell.isNeighbor(T.background()))) {
      // Calculate cell position and size from vertices
      std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
      assert(xCenter.size() == dimension);
      double cellVolume = tmpCell.calculateVolume(vertexData);
      double fac = parameter(0) * 0.5;
      if (parameter(1) == 1) {
        fac /= cellVolume;
      }
      // Calculate derivative contributions to vertices from each wall
      for (size_t k = 0; k < tmpCell.numWall(); ++k) {
        Wall &tmpWall = tmpCell.wallRef(k);
        size_t v1I = tmpWall.vertex1()->index();
        size_t v2I = tmpWall.vertex2()->index();
        std::vector<double> n(dimension), dx(dimension), x0(dimension);
        double b = 0;
        for (size_t d = 0; d < dimension; ++d) {
          n[d] = vertexData[v2I][d] - vertexData[v1I][d];
          b += n[d] * n[d];
          x0[d] = 0.5 * (vertexData[v1I][d] + vertexData[v2I][d]);
          dx[d] = xCenter[d] - x0[d];
        }
        assert(b > 0.0);
        b = std::sqrt(b);
        for (size_t d = 0; d < dimension; ++d) n[d] /= b;
        double bInv = 1.0 / b;
        double h =
            dx[0] * dx[0] + dx[1] * dx[1] -
            (n[0] * dx[0] + n[1] * dx[1]) * (n[0] * dx[0] + n[1] * dx[1]);
        assert(h > 0.0);
        h = std::sqrt(h);
        double hInv = 1.0 / h;
        vertexDerivs[v1I][0] +=
            fac * (0.5 * b * hInv *
                       (-dx[0] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (-n[1] * n[1] * bInv * dx[0] - 0.5 * n[0] +
                                      n[0] * n[1] * bInv * dx[1])) -
                   h * n[0]);
        vertexDerivs[v2I][0] +=
            fac * (0.5 * b * hInv *
                       (-dx[0] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (n[1] * n[1] * bInv * dx[0] - 0.5 * n[0] -
                                      n[0] * n[1] * bInv * dx[1])) +
                   h * n[0]);
        vertexDerivs[v1I][1] +=
            fac * (0.5 * b * hInv *
                       (-dx[1] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (-n[0] * n[0] * bInv * dx[1] - 0.5 * n[1] +
                                      n[0] * n[1] * bInv * dx[0])) -
                   h * n[1]);
        vertexDerivs[v2I][1] +=
            fac * (0.5 * b * hInv *
                       (-dx[1] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (n[0] * n[0] * bInv * dx[1] - 0.5 * n[1] -
                                      n[0] * n[1] * bInv * dx[0])) +
                   h * n[1]);
      }
    }
  }
}

AreaPotentialTriSpatialThreshold::AreaPotentialTriSpatialThreshold(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t>> &indValue) {
  // Do some checks on the parameters and variable indeces
  //
  if (paraValue.size() != 2) {
    std::cerr << "Pressure2D::AreaPotentialTriSpatialThreshold::"
              << "AreaPotentialTriSpatialThreshold() "
              << "Uses two parameters P_force and X_th.\n";
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "AreaPotentialTriSpatialThreshold::"
              << "AreaPotentialTriSpatialThreshold() "
              << "One index given (direction for threshold).\n";
    exit(EXIT_FAILURE);
  }
  // Set the variable values
  //
  setId("Pressure2D::AreaPotentialTriSpatialThreshold");
  setParameter(paraValue);
  setVariableIndex(indValue);

  // Set the parameter identities
  //
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "K_force";
  tmp[1] = "X_th";
  setParameterId(tmp);
}

void AreaPotentialTriSpatialThreshold::derivs(Tissue &T, DataMatrix &cellData,
                                           DataMatrix &wallData,
                                           DataMatrix &vertexData,
                                           DataMatrix &cellDerivs,
                                           DataMatrix &wallDerivs,
                                           DataMatrix &vertexDerivs) {
  // Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition();
  size_t maxDim = variableIndex(0, 0);
  assert(maxDim < dimension);

  // Find max pos in given direction;
  double max = vertexData[0][maxDim];
  for (size_t i = 1; i < vertexData.size(); ++i) {
    if (vertexData[i][maxDim] > max) max = vertexData[i][maxDim];
  }
  // For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {
    Cell &tmpCell = T.cell(cellI);
    // Calculate cell position from vertices
    std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
    assert(xCenter.size() == dimension);
    // Only update if close to apex
    if (max - xCenter[maxDim] < parameter(1)) {
      // Calculate derivative contributions to vertices from each wall
      for (size_t k = 0; k < tmpCell.numWall(); ++k) {
        Wall &tmpWall = tmpCell.wallRef(k);
        size_t v1I = tmpWall.vertex1()->index();
        size_t v2I = tmpWall.vertex2()->index();
        std::vector<double> n(dimension), dx(dimension), x0(dimension);
        double b = 0;
        for (size_t d = 0; d < dimension; ++d) {
          n[d] = vertexData[v2I][d] - vertexData[v1I][d];
          b += n[d] * n[d];
          x0[d] = 0.5 * (vertexData[v1I][d] + vertexData[v2I][d]);
          dx[d] = xCenter[d] - x0[d];
        }
        assert(b > 0.0);
        b = std::sqrt(b);
        for (size_t d = 0; d < dimension; ++d) n[d] /= b;
        double bInv = 1.0 / b;
        double h =
            dx[0] * dx[0] + dx[1] * dx[1] -
            (n[0] * dx[0] + n[1] * dx[1]) * (n[0] * dx[0] + n[1] * dx[1]);
        assert(h > 0.0);
        h = std::sqrt(h);
        double hInv = 1.0 / h;
        double fac = parameter(0) * 0.5;

        vertexDerivs[v1I][0] +=
            fac * (0.5 * b * hInv *
                       (-dx[0] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (-n[1] * n[1] * bInv * dx[0] - 0.5 * n[0] +
                                      n[0] * n[1] * bInv * dx[1])) -
                   h * n[0]);
        vertexDerivs[v2I][0] +=
            fac * (0.5 * b * hInv *
                       (-dx[0] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (n[1] * n[1] * bInv * dx[0] - 0.5 * n[0] -
                                      n[0] * n[1] * bInv * dx[1])) +
                   h * n[0]);
        vertexDerivs[v1I][1] +=
            fac * (0.5 * b * hInv *
                       (-dx[1] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (-n[0] * n[0] * bInv * dx[1] - 0.5 * n[1] +
                                      n[0] * n[1] * bInv * dx[0])) -
                   h * n[1]);
        vertexDerivs[v2I][1] +=
            fac * (0.5 * b * hInv *
                       (-dx[1] - 2.0 * (n[0] * dx[0] + n[1] * dx[1]) *
                                     (n[0] * n[0] * bInv * dx[1] - 0.5 * n[1] -
                                      n[0] * n[1] * bInv * dx[0])) +
                   h * n[1]);

        // vertexDerivs[v1I][0] += fac*( b*hInv*(bInv*(n[0]*dx[0]+n[1]*dx[1])*
        //				    (n[1]*n[1]*dx[0]+b*n[0]-n[0]*n[1]*dx[1]) +
        //dx[0] ) - h*n[0] );  vertexDerivs[v1I][1] += fac*(
        // b*hInv*(bInv*(n[1]*dx[1]+n[0]*dx[0])*
        //				    (-n[0]*n[0]*dx[1]+b*n[1]+n[1]*n[0]*dx[0])
        //+ dx[1] ) - h*n[1] );  vertexDerivs[v2I][0] += fac*( n[0]*h +
        // n[1]*hInv*(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[1]-n[1]*dx[0]) );
        // vertexDerivs[v2I][1] += fac*( n[1]*h +
        // n[0]*hInv*(n[1]*dx[1]+n[0]*dx[0])*(n[1]*dx[0]-n[0]*dx[1]) );
      }
    }
  }
}

AreaPotentialTargetArea::AreaPotentialTargetArea(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t>> &indValue) {
  if (paraValue.size() != 2) {
    std::cerr << "AreaPotentialTargetArea::AreaPotentialTargetArea() "
              << "Uses two parameter: P no_contraction_flag" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "AreaPotentialTargetArea::AreaPotentialTargetArea() "
              << "No target volume index given.\n";
    exit(EXIT_FAILURE);
  }

  setId("Pressure2D::AreaPotentialTargetArea");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "P";
  tmp[1] = "f_no_contraction";
  setParameterId(tmp);
}

void AreaPotentialTargetArea::derivs(Tissue &T, DataMatrix &cellData,
                                     DataMatrix &wallData,
                                     DataMatrix &vertexData,
                                     DataMatrix &cellDerivs,
                                     DataMatrix &wallDerivs,
                                     DataMatrix &vertexDerivs) {
  // This function assumes two dimensions (x, y) and that
  // the vertices are sorted in a clock-wise order.

  assert(T.vertex(0).numPosition());

  double area;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    std::vector<std::pair<double, double>> vertices;

    for (size_t i = 0; i < cell.numVertex(); ++i) {
      std::pair<double, double> vertex;
      vertex.first  = vertexData[cell.vertex(i)->index()][0];
      vertex.second = vertexData[cell.vertex(i)->index()][1];
      vertices.push_back(vertex);
    }
    area = polygonArea(vertices);
    double sa = area < 0 ? -1 : +1;
    area = std::fabs(area);

    // Don't allow lower water volume than volume
    if (parameter(1) != 0.0 && area > cellData[n][variableIndex(0, 0)]) {
      continue;
    }

    for (size_t i = 1; i < (cell.numVertex() + 1); ++i) {
      Vertex *vertex  = cell.vertex(i % cell.numVertex());       // Current vertex in polygon.
      Vertex *pvertex = cell.vertex((i - 1) % cell.numVertex()); // Previous vertex in polygon.
      Vertex *nvertex = cell.vertex((i + 1) % cell.numVertex()); // Next vertex in polygon.

      double px = vertexData[pvertex->index()][0];
      double py = vertexData[pvertex->index()][1];

      double nx = vertexData[nvertex->index()][0];
      double ny = vertexData[nvertex->index()][1];

      double dAdx = sa * 0.5 * (-py + ny);
      double dAdy = sa * 0.5 * (px - nx);

      vertexDerivs[vertex->index()][0] +=
          parameter(0) * (1 - area / cellData[n][variableIndex(0, 0)]) * dAdx;
      vertexDerivs[vertex->index()][1] +=
          parameter(0) * (1 - area / cellData[n][variableIndex(0, 0)]) * dAdy;
    }
  }
}

double AreaPotentialTargetArea::polygonArea(
    std::vector<std::pair<double, double>> vertices) {
  double area = 0.0;
  size_t N = vertices.size();
  for (size_t n = 0; n < N; ++n) {
    area += vertices[n].first * vertices[(n + 1) % N].second;
    area -= vertices[(n + 1) % N].first * vertices[n].second;
  }
  area *= 0.5;
  return area;
}

}  // end namespace Pressure2D
