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
           std::vector<std::vector<size_t> >
               &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::Radial::"
                  << "Radial() "
                  << "Uses two parameters k_growth and r_pow (0,1)" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
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
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[0] = "r_pow";
    setParameterId(tmp);
}

void Radial::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    size_t dimension = vertexData[0].size();

    for (size_t i = 0; i < numVertices; ++i) {
        double fac = parameter(0);
        if (parameter(1) == 0.0) {
            double r = 0.0;
            for (size_t d = 0; d < dimension; ++d)
                r += vertexData[i][d] * vertexData[i][d];
            if (r > 0.0)
                r = std::sqrt(r);
            if (r > 0.0)
                fac /= r;
            else
                fac = 0.0;
        }
        for (size_t d = 0; d < dimension; ++d)
            vertexDerivs[i][d] += fac * vertexData[i][d];
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
                           DataMatrix &sdydtVertex) {
    return derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
}

namespace CenterTriangulation {
Radial::
    Radial(std::vector<double> &paraValue,
           std::vector<std::vector<size_t> >
               &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::CenterTriangulation::Radial::"
                  << "Radial() " << std::endl
                  << "Uses two parameters k_growth and r_pow (0,1)" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
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
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[0] = "r_pow";
    setParameterId(tmp);
}

void Radial::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    size_t numCells = T.numCell();
    size_t dimension = vertexData[0].size();

    // Move vertices / apply force to vertices
    for (size_t i = 0; i < numVertices; ++i) {
        double fac = parameter(0);
        if (parameter(1) == 0.0) {
            double r = 0.0;
            for (size_t d = 0; d < dimension; ++d)
                r += vertexData[i][d] * vertexData[i][d];
            if (r > 0.0)
                r = std::sqrt(r);
            if (r > 0.0)
                fac /= r;
            else
                fac = 0.0;
        }
        for (size_t d = 0; d < dimension; ++d)
            vertexDerivs[i][d] += fac * vertexData[i][d];
    }
    // Move / apply force to vertices defined in cell centers
    for (size_t i = 0; i < numCells; ++i) {
        double fac = parameter(0);
        if (parameter(1) == 0.0) {
            double r = 0.0;
            for (size_t d = variableIndex(0, 0); d < variableIndex(0, 0) + dimension; ++d)
                r += cellData[i][d] * cellData[i][d];
            if (r > 0.0)
                r = std::sqrt(r);
            if (r > 0.0)
                fac /= r;
            else
                fac = 0.0;
        }
        for (size_t d = variableIndex(0, 0); d < variableIndex(0, 0) + dimension; ++d)
            cellDerivs[i][d] += fac * cellData[i][d];
    }
}

ForceToCell::
    ForceToCell(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> >
                    &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::CenterTriangulation::ForceToCell::"
                  << "ForceToCell() " << std::endl
                  << "Uses two parameters k_growth and r_pow (0,1)" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 2) {
        std::cerr << "GrowthForce::CenterTriangulation::ForceToCell::"
                  << "ForceToCell() " << std::endl
                  << "CellData index of cell to point force towards and"
                  << "Start of additional Cell variable indices (center(x,y,z) "
                  << "L_1,...,L_n, n=num vertex) are given in first level. "
                  << "See Documentation for namespace CenterTriangulation."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("GrowthForce::CenterTriangulation::ForceToCell");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[0] = "r_pow";
    setParameterId(tmp);
}

void ForceToCell::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    size_t numCells = T.numCell();
    size_t dimension = vertexData[0].size();
    size_t forceCellData = variableIndex(0, 1);

    double fac = parameter(0);
    int v = 0;
    // Move vertices / apply force to vertices
    for (size_t i = 0; i < numVertices; ++i) {
        std::vector<double> force{0.0, 0.0, 0.0};
        std::vector<double> vertPos{vertexData[i][0], vertexData[i][1], vertexData[i][2]};
        for (size_t c = 0; c < T.vertex(i).numCell(); c++) {
            size_t cid = T.vertex(i).cell(c)->index();
            int abv_id = cellData[cid][forceCellData];
            if (abv_id > 0) {
                std::vector<double> abv_ctrd = T.cell(abv_id).positionFromVertex(vertexData);
                for (size_t d = 0; d < dimension; d++) {
                    force[d] += (abv_ctrd[d] - vertPos[d]);
                }
            }
        }
        double r = 0.0;
        for (size_t d = 0; d < dimension; d++)
            r += force[d] * force[d];
        if (r > 0.0)
            r = std::sqrt(r);
        if (r > 0.0)
            fac /= r;
        else
            fac = 0.0;

        for (size_t d = 0; d < dimension; ++d)
            vertexDerivs[i][d] += fac * force[d];
    }

    // Move / apply force to vertices defined in cell centers
    for (size_t i = 0; i < numCells; ++i) {
        double fac = parameter(0);
        std::vector<double> force{0.0, 0.0, 0.0};
        int abv_id = cellData[i][forceCellData];
        std::vector<double> vertPos{cellData[i][variableIndex(0, 0) + 0],
                                    cellData[i][variableIndex(0, 0) + 1],
                                    cellData[i][variableIndex(0, 0) + 2]};
        if (abv_id > 0) {
            std::vector<double> abv_ctrd = T.cell(abv_id).positionFromVertex(vertexData);
            for (size_t d = 0; d < dimension; d++) {
                force[d] += (abv_ctrd[d] - vertPos[d]);
            }
        }
        double r = 0.0;
        for (size_t d = variableIndex(0, 0); d < variableIndex(0, 0) + dimension; ++d)
            r += force[d] * force[d];
        if (r > 0.0)
            r = std::sqrt(r);
        if (r > 0.0)
            fac /= r;
        else
            fac = 0.0;

        for (size_t d = 0; d < dimension; ++d)
            cellDerivs[i][variableIndex(0, 0) + d] += fac * force[d];
    }
}

}  // end namespace CenterTriangulation

EpidermalRadial::
    EpidermalRadial(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> >
                        &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::EpidermalRadial::"
                  << "EpidermalRadial() "
                  << "Uses two parameters k_growth and r_pow (0,1)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr << "GrowthForce::EpidermalRadial::"
                  << "EpidermalRadial() "
                  << "No variable index is used."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("GrowthForce::EpidermalRadial");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[0] = "r_pow";
    setParameterId(tmp);
}

void EpidermalRadial::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    size_t dimension = vertexData[0].size();

    for (size_t i = 0; i < numVertices; ++i) {
        if (T.vertex(i).isBoundary(T.background())) {
            double fac = parameter(0);
            if (parameter(1) == 0.0) {
                double r = 0.0;
                for (size_t d = 0; d < dimension; ++d)
                    r += vertexData[i][d] * vertexData[i][d];
                if (r > 0.0)
                    r = std::sqrt(r);
                if (r > 0.0)
                    fac /= r;
                else
                    fac = 0.0;
            }
            for (size_t d = 0; d < dimension; ++d)
                vertexDerivs[i][d] += fac * vertexData[i][d];
        }
    }
}

X::
    X(std::vector<double> &paraValue,
      std::vector<std::vector<size_t> >
          &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::X::"
                  << "X() "
                  << "Uses two parameters k_growth and growth_mode (0,1)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr << "GrowthForce::X::"
                  << "X() "
                  << "No variable index is used." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("GrowthForce::X");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[1] = "growth_mode";
    setParameterId(tmp);
}

void X::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    size_t s_i = 0;  // spatial index
    double fac = parameter(0);
    // size_t dimension=vertexData[s_i].size();
    size_t growth_mode = parameter(1);

    for (size_t i = 0; i < numVertices; ++i) {
        if (growth_mode == 1) {
            vertexDerivs[i][s_i] += fac * vertexData[i][s_i];
        } else {
            if (vertexData[i][s_i] >= 0) {
                vertexDerivs[i][s_i] += fac;
            } else {
                vertexDerivs[i][s_i] -= fac;
            }
        }
    }
}

void X::
    derivsWithAbs(Tissue &T,
                  DataMatrix &cellData,
                  DataMatrix &wallData,
                  DataMatrix &vertexData,
                  DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs,
                  DataMatrix &vertexDerivs,
                  DataMatrix &sdydtCell,
                  DataMatrix &sdydtWall,
                  DataMatrix &sdydtVertex) {
    return derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
}

Y::
    Y(std::vector<double> &paraValue,
      std::vector<std::vector<size_t> >
          &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::Y::"
                  << "Y() "
                  << "Uses two parameters k_growth and growth_mode (=0/1)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr << "GrowthForce::Y::"
                  << "Y() "
                  << "No variable index is used." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("GrowthForce::Y");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[1] = "growth_mode";
    setParameterId(tmp);
}

void Y::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    size_t s_i = 1;  // spatial index
    // size_t dimension=vertexData[s_i].size();
    double fac = parameter(0);
    size_t growth_mode = parameter(1);
    // std::cout <<  "fac = " << fac << "\n";

    for (size_t i = 0; i < numVertices; ++i) {
        if (growth_mode == 1) {
            vertexDerivs[i][s_i] += fac * vertexData[i][s_i];
        } else {
            if (vertexData[i][s_i] >= 0) {
                vertexDerivs[i][s_i] += fac;
            } else {
                vertexDerivs[i][s_i] -= fac;
            }
        }
    }
}

void Y::
    derivsWithAbs(Tissue &T,
                  DataMatrix &cellData,
                  DataMatrix &wallData,
                  DataMatrix &vertexData,
                  DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs,
                  DataMatrix &vertexDerivs,
                  DataMatrix &sdydtCell,
                  DataMatrix &sdydtWall,
                  DataMatrix &sdydtVertex) {
    return derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
}

SphereCylinder::
    SphereCylinder(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> >
                       &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0 && paraValue[1] != 1)) {
        std::cerr << "GrowthForce::SphereCylinder::"
                  << "SphereCylinder() "
                  << "Uses two parameters k_growth and r_pow (0,1)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr << "GrowthForce::SphereCylinder::"
                  << "SphereCylinder() "
                  << "No variable index is used." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("GrowthForce::SphereCylinder");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "k_growth";
    tmp[0] = "r_pow";
    setParameterId(tmp);
}

void SphereCylinder::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    if (vertexData[0].size() != 3) {
        std::cerr << "GrowthForceSphereCylinder::derivs() "
                  << "Only works for 3 dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t xI = 0;
    size_t yI = 1;
    size_t zI = 2;

    for (size_t i = 0; i < numVertices; ++i) {
        if (vertexData[i][zI] < 0.0) {  // on cylinder
            if (parameter(1) == 0.0) {
                vertexDerivs[i][zI] -= parameter(0);
            } else {
                double r = std::sqrt(vertexData[i][xI] * vertexData[i][xI] +
                                     vertexData[i][yI] * vertexData[i][yI]);
                vertexDerivs[i][zI] -= parameter(0) * (3.14159265 * 0.5 * r - vertexData[i][zI]);
            }
        } else {  // on half sphere
            double r = std::sqrt(vertexData[i][xI] * vertexData[i][xI] +
                                 vertexData[i][yI] * vertexData[i][yI] +
                                 vertexData[i][zI] * vertexData[i][zI]);
            double rPrime = std::sqrt(vertexData[i][xI] * vertexData[i][xI] +
                                      vertexData[i][yI] * vertexData[i][yI]);
            double theta = std::asin(rPrime / r);

            double fac = parameter(0) * theta;
            if (parameter(0) == 1) {
                fac *= r;
            }
            vertexDerivs[i][xI] += fac * vertexData[i][xI] * vertexData[i][zI] / rPrime;
            vertexDerivs[i][yI] += fac * vertexData[i][yI] * vertexData[i][zI] / rPrime;
            vertexDerivs[i][zI] -= fac * rPrime;
        }
    }
}
}  // end namespace GrowthForce
