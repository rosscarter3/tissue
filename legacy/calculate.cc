//
// Filename     : calculate.cc
// Description  : Classes for calculating and updating variables
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : April 2018
// Revision     : $Id:$
//

#include "calculate.h"

#include <vector>

#include "baseReaction.h"
#include "tissue.h"

namespace Calculate {

AngleVectors::
    AngleVectors(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> >
            &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "Calculate::AngleVectors:: "
                  << "AngleVectors() "
                  << "Calculates abs(cos(...)) of angle between two 3d "
                  << "vectors(starting from given indices) in cellData vector "
                  << "and stores it in the given index in cellData vector, uses no "
                  << "parameter " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "Calculate::AngleVectors:: "
                  << "AngleVectors() "
                  << "1st level with 2 variable indices used: "
                  << "start index for 1st vector   start index for 2nd vector "
                  << "2nd level with 1 variable index used: "
                  << "store index for angle(deg) " << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Calculate::AngleVectors");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void AngleVectors::
    derivs(Tissue &T, DataMatrix &cellData,
           DataMatrix &wallData, DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
        double teta = 0;
        double tmp = std::sqrt(cellData[cellIndex][variableIndex(0, 0)] *
                                   cellData[cellIndex][variableIndex(0, 0)] +
                               cellData[cellIndex][variableIndex(0, 0) + 1] *
                                   cellData[cellIndex][variableIndex(0, 0) + 1] +
                               cellData[cellIndex][variableIndex(0, 0) + 2] *
                                   cellData[cellIndex][variableIndex(0, 0) + 2]);
        if (tmp > 0) {
            cellData[cellIndex][variableIndex(0, 0)] =
                cellData[cellIndex][variableIndex(0, 0)] / tmp;
            cellData[cellIndex][variableIndex(0, 0) + 1] =
                cellData[cellIndex][variableIndex(0, 0) + 1] / tmp;
            cellData[cellIndex][variableIndex(0, 0) + 2] =
                cellData[cellIndex][variableIndex(0, 0) + 2] / tmp;
        }

        tmp = std::sqrt(cellData[cellIndex][variableIndex(0, 1)] *
                            cellData[cellIndex][variableIndex(0, 1)] +
                        cellData[cellIndex][variableIndex(0, 1) + 1] *
                            cellData[cellIndex][variableIndex(0, 1) + 1] +
                        cellData[cellIndex][variableIndex(0, 1) + 2] *
                            cellData[cellIndex][variableIndex(0, 1) + 2]);
        if (tmp > 0) {
            cellData[cellIndex][variableIndex(0, 1)] =
                cellData[cellIndex][variableIndex(0, 1)] / tmp;
            cellData[cellIndex][variableIndex(0, 1) + 1] =
                cellData[cellIndex][variableIndex(0, 1) + 1] / tmp;
            cellData[cellIndex][variableIndex(0, 1) + 2] =
                cellData[cellIndex][variableIndex(0, 1) + 2] / tmp;
        }
        teta = std::abs(cellData[cellIndex][variableIndex(0, 0)] *
                            cellData[cellIndex][variableIndex(0, 1)] +
                        cellData[cellIndex][variableIndex(0, 0) + 1] *
                            cellData[cellIndex][variableIndex(0, 1) + 1] +
                        cellData[cellIndex][variableIndex(0, 0) + 2] *
                            cellData[cellIndex][variableIndex(0, 1) + 2]);

        cellData[cellIndex][variableIndex(1, 0)] = teta;
    }
}

AngleVectorXYplane::
    AngleVectorXYplane(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> >
            &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "Calculate::AngleVectorXYplane:: "
                  << "AngleVectorXYplane() "
                  << "Calculates abs(cos(...)) of angle between a vector and "
                  << "XY plane and stores it in the given index in cellData "
                  << "vector, uses no parameter." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 ||
        indValue[1].size() != 1) {
        std::cerr << "Calculate::AngleVectorXYplane:: "
                  << "AngleVectorXYplane() "
                  << "1st level with 1 variable index used: "
                  << "start index for the vector  "
                  << "2nd level with 1 variable index used: "
                  << "store index for cos(angle) " << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Calculate::AngleVectorXYplane");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void AngleVectorXYplane::
    derivs(Tissue &T, DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
        double teta = 0;
        double temp = std::sqrt(cellData[cellIndex][variableIndex(0, 0)] *
                                    cellData[cellIndex][variableIndex(0, 0)] +
                                cellData[cellIndex][variableIndex(0, 0) + 1] *
                                    cellData[cellIndex][variableIndex(0, 0) + 1]);
        double pi = 3.14159265;
        if (temp < 0.00000001)
            teta = pi / 2;
        else
            teta = std::atan(cellData[cellIndex][variableIndex(0, 0) + 2] / temp);

        cellData[cellIndex][variableIndex(1, 0)] = teta;
    }
}

AngleVector::
    AngleVector(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> >
            &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "Calculate::AngleVector:: "
                  << "AngleVector() "
                  << "Calculates  angle between a 3d vector(starting from given "
                  << "index) in cellData vector and "
                  << "a given axes(x,y,z) and stores it in the given index in "
                  << "cellData vector, uses one parameter "
                  << "for specifying the axes " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 ||
        indValue[1].size() != 1) {
        std::cerr << "Calculate::AngleVector:: "
                  << "AngleVector() "
                  << "1st level with 1 variable index used: "
                  << "start index for the vector  "
                  << "2nd level with 1 variable index used: "
                  << "store index for angle(deg) " << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Calculate::AngleVector");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "axis";  // axis specifier

    setParameterId(tmp);

    if (parameter(0) != 0 && parameter(0) != 1 && parameter(0) != 2) {
        std::cerr << "Calculate::AngleVector::AngleVector() "
                  << "Parameter 0 can only be 0(X), 1(Y) or 2(Z) " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (parameter(0) == 1 || parameter(0) == 2) {
        std::cerr << "Calculate::AngleVector::AngleVector() "
                  << "The code should be modified for 3d " << std::endl;
        exit(0);
    }
}

void AngleVector::
    derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
           DataMatrix &vertexData, DataMatrix &cellDerivs,
           DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
        double teta = 0;
        double tmp = std::sqrt(cellData[cellIndex][variableIndex(0, 0)] *
                                   cellData[cellIndex][variableIndex(0, 0)] +
                               cellData[cellIndex][variableIndex(0, 0) + 1] *
                                   cellData[cellIndex][variableIndex(0, 0) + 1] +
                               cellData[cellIndex][variableIndex(0, 0) + 2] *
                                   cellData[cellIndex][variableIndex(0, 0) + 2]);
        cellData[cellIndex][variableIndex(0, 0)] =
            cellData[cellIndex][variableIndex(0, 0)] / tmp;
        cellData[cellIndex][variableIndex(0, 0) + 1] =
            cellData[cellIndex][variableIndex(0, 0) + 1] / tmp;
        cellData[cellIndex][variableIndex(0, 0) + 2] =
            cellData[cellIndex][variableIndex(0, 0) + 2] / tmp;

        double pi = 3.141592;
        if (parameter(0) == 0) {
            teta = 180 * std::acos(cellData[cellIndex][variableIndex(0, 0)]) / pi;
            if (cellData[cellIndex][variableIndex(0, 0) + 1] < 0) {
                teta = 180 - teta;
            }
            cellData[cellIndex][variableIndex(1, 0)] = teta;
        }

        // should be modified for 3d
        if (parameter(0) == 1) {
            teta = 180 * std::acos(cellData[cellIndex][variableIndex(0, 0) + 1]) / pi;
            cellData[cellIndex][variableIndex(1, 0)] = teta;
        }

        if (parameter(0) == 2) {
            teta = 180 * std::acos(cellData[cellIndex][variableIndex(0, 0) + 2]) / pi;
            cellData[cellIndex][variableIndex(1, 0)] = teta;
        }
    }
}

VertexVelocity::
    VertexVelocity(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> >
            &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "Calculate::VertexVelocity:: "
                  << "VertexVelocity() "
                  << "Calculates  average velocity of vertices per face/cell in tissue "
                  << "for e.g. checking closeness to mechanical equilibrum."
                  << "Uses no parameter. " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Calculate::VertexVelocity:: "
                  << "VertexVelocity() "
                  << "1 level with 1 variable index used: "
                  << "store index for vertex velocity. " << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Calculate::VertexVelocity");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    std::vector<std::string> tmp(numParameter());
    setParameterId(tmp);
}

void VertexVelocity::
    derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
           DataMatrix &vertexData, DataMatrix &cellDerivs,
           DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t velocityIndex = variableIndex(0, 0);
    size_t dimension = vertexData[0].size();

    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
        size_t numCellVertex = T.cell(cellIndex).numVertex();
        double cellVelocity = 0.;
        for (size_t k = 0; k < numCellVertex; ++k) {
            size_t VIndex = T.cell(cellIndex).vertex(k)->index();
            double vertexVelocity = 0.;
            for (size_t d = 0; d < dimension; ++d) {
                vertexVelocity += vertexDerivs[VIndex][d] * vertexDerivs[VIndex][d];
            }
            cellVelocity += std::sqrt(vertexVelocity);
        }
        cellData[cellIndex][velocityIndex] = cellVelocity / numCellVertex;
    }
}

TissueVolumeChange::
    TissueVolumeChange(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> >
            &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "Calculate::TissueVolumeChange:: "
                  << "TissueVolumeChange() "
                  << "Calculates change in template volume and its time derivative "
                  << "and total Derivative and stores them in the given indices in "
                  << "cellData vector "
                  << "uses no parameter " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 6) {
        std::cerr << "Calculate::TissueVolumeChange:: "
                  << "TissueVolumeChange() "
                  << "1 level with 6 variable indices used: "
                  << "cell-index-VolumeChange       component-index-VolumeChange "
                  << "cell-index-deltaVolumeChange  component-index-deltaVolumeChange "
                  << "cell-index-totalDerivative    component-index-totalDerivative "
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Calculate::TissueVolumeChange");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void TissueVolumeChange::
    initiate(Tissue &T, DataMatrix &cellData,
             DataMatrix &wallData,
             DataMatrix &vertexData,
             DataMatrix &cellDerivs,
             DataMatrix &wallDerivs,
             DataMatrix &vertexDerivs) {
    size_t dimension = 3;  // Only implemented for 3D models
    assert(dimension == vertexData[0].size());
    vertexDataRest.resize(vertexData.size());
    for (size_t i = 0; i < vertexDataRest.size(); ++i)
        vertexDataRest[i].resize(vertexData[i].size());

    size_t numVertices = T.numVertex();

    for (size_t vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex)
        for (size_t d = 0; d < dimension; ++d) {
            vertexDataRest[vertexIndex][d] = vertexData[vertexIndex][d];
        }
}

void TissueVolumeChange::
    derivs(Tissue &T, DataMatrix &cellData,
           DataMatrix &wallData, DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numVertices = T.numVertex();
    VolumeChange = 0;
    for (size_t vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) {
        VolumeChange += std::sqrt(
            (vertexData[vertexIndex][0] -
             vertexDataRest[vertexIndex][0]) *
                (vertexData[vertexIndex][0] -
                 vertexDataRest[vertexIndex][0]) +
            (vertexData[vertexIndex][1] -
             vertexDataRest[vertexIndex][1]) *
                (vertexData[vertexIndex][1] -
                 vertexDataRest[vertexIndex][1]) +
            (vertexData[vertexIndex][2] -
             vertexDataRest[vertexIndex][2]) *
                (vertexData[vertexIndex][2] -
                 vertexDataRest[vertexIndex][2]));
    }
    deltaVolumeChange =
        VolumeChange - cellData[variableIndex(0, 0)][variableIndex(0, 1)];
    totalDerivative = 0;
    for (size_t vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) {
        if (vertexDerivs[0].size() > 2)
            totalDerivative += std::sqrt(
                vertexDerivs[vertexIndex][0] *
                    vertexDerivs[vertexIndex][0] +
                vertexDerivs[vertexIndex][1] *
                    vertexDerivs[vertexIndex][1] +
                vertexDerivs[vertexIndex][2] *
                    vertexDerivs[vertexIndex][2]);
    }

    cellData[variableIndex(0, 0)][variableIndex(0, 1)] = VolumeChange;
    cellData[variableIndex(0, 2)][variableIndex(0, 3)] = deltaVolumeChange;
    cellData[variableIndex(0, 4)][variableIndex(0, 5)] = totalDerivative;
    // std:: cerr << VolumeChange<<" "<<deltaVolumeChange<<"
    // "<<VolumeChangeDerT<<std::endl;
}

}  // end namespace Calculate
