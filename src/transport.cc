//
// Filename     : transport.cc
// Description  : Classes describing transport reactions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2013
// Revision     : $Id:$
//

#include "transport.h"

MembraneDiffusionSimple::
    MembraneDiffusionSimple(std::vector<double> &paraValue,
                            std::vector<std::vector<size_t> >
                                &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "MembraneDiffusionSimple::"
                  << "MembraneDiffusionSimple() "
                  << "One parameters (Diffusion constant) used." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "MembraneDiffusionSimple::"
                  << "MembraneDiffusionSimple() "
                  << "One variable (pair) used for diffusing wall variable." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("MembraneDiffusionSimple");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "D_mem";

    setParameterId(tmp);
}

void MembraneDiffusionSimple::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t pwI = variableIndex(0, 0);  // diffusing molecule (membrane/wall)

    assert(pwI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();
        for (size_t k = 0; k < numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();

            size_t pwIadd = 1;  // Two variables per wall, setting which to use
            if (T.cell(i).wall(k)->cell1()->index() == i) {
                pwIadd = 0;
            }
            double fac = parameter(0) * wallData[j][pwI + pwIadd];
            wallDerivs[j][pwI + pwIadd] -= 2. * fac;

            size_t kNext = (k + 1) % numWalls;
            size_t jNext = T.cell(i).wall(kNext)->index();
            if (T.cell(i).wall(kNext)->cell1()->index() == i) {
                wallDerivs[jNext][pwI] += fac;
            } else {
                wallDerivs[jNext][pwI + 1] += fac;
            }
            size_t kBefore = k > 0 ? k - 1 : numWalls - 1;
            size_t jBefore = T.cell(i).wall(kBefore)->index();
            if (T.cell(i).wall(kBefore)->cell1()->index() == i) {
                wallDerivs[jBefore][pwI] += fac;
            } else {
                wallDerivs[jBefore][pwI + 1] += fac;
            }
        }
    }
}

DiffusionSimple::
    DiffusionSimple(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> >
                        &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "DiffusionModelSimple::"
                  << "DiffusionModelSimple() "
                  << "One parameter (diffusion constant) used: p_0" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "DiffusionSimple::"
                  << "DiffusionSimple() "
                  << "One level of one variable indices used" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("DiffusionSimple");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "p_0";

    setParameterId(tmp);
}

void DiffusionSimple::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);
    assert(aI < cellData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }
                if (i < neighIndex) {  // Both directions at once
                    cellDerivs[i][aI] -= parameter(0) * (cellData[i][aI] - cellData[neighIndex][aI]);
                    cellDerivs[neighIndex][aI] += parameter(0) * (cellData[i][aI] - cellData[neighIndex][aI]);
                }
            }
        }
    }
}

void DiffusionSimple::
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
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);
    assert(aI < cellData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }
                if (i < neighIndex) {  // Both directions at once
                    cellDerivs[i][aI] -= parameter(0) * (cellData[i][aI] - cellData[neighIndex][aI]);
                    cellDerivs[neighIndex][aI] += parameter(0) * (cellData[i][aI] - cellData[neighIndex][aI]);
                    sdydtCell[i][aI] += 0;
                    sdydtCell[neighIndex][aI] += 0;
                }
            }
        }
    }
}

DiffusionSimpleOne::
    DiffusionSimpleOne(std::vector<double> &paraValue,
                       std::vector<std::vector<size_t> >
                           &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "DiffusionSimpleOne::"
                  << "DiffusionSimpleOne() "
                  << "One parameter (diffusion constant) used: p_0" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) {
        std::cerr << "DiffusionSimpleOne::"
                  << "DiffusionSimpleOne() "
                  << "Two levels of one variable index each is used" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Diffusion::SimpleOne");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "p_0";

    setParameterId(tmp);
}

void DiffusionSimpleOne::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);
    size_t hI = variableIndex(1, 0);
    assert(aI < cellData[0].size() && hI < cellData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }
                if (i < neighIndex) {  // Both directions at once
                    double fac = parameter(0) * (cellData[i][aI] * cellData[i][hI] -
                                                 cellData[neighIndex][aI] * cellData[neighIndex][hI]);
                    cellDerivs[i][aI] -= fac;
                    cellDerivs[neighIndex][aI] += fac;
                }
            }
        }
    }
}

void DiffusionSimpleOne::
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
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);
    size_t hI = variableIndex(1, 0);
    assert(aI < cellData[0].size() && hI < cellData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }
                if (i < neighIndex) {  // Both directions at once
                    double fac = parameter(0) * (cellData[i][aI] * cellData[i][hI] -
                                                 cellData[neighIndex][aI] * cellData[neighIndex][hI]);
                    cellDerivs[i][aI] -= fac;
                    cellDerivs[neighIndex][aI] += fac;
                    sdydtCell[i][aI] += 0.0;
                    sdydtCell[neighIndex][aI] += 0.0;
                }
            }
        }
    }
}

DiffusionConductiveSimple::
    DiffusionConductiveSimple(std::vector<double> &paraValue,
                              std::vector<std::vector<size_t> >
                                  &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 5) {
        std::cerr << "DiffusionConductiveSimple::"
                  << "DiffusionConductiveSimple() "
                  << "Five parameters used (see documentation)" << std::endl
                  << "p_0 - diffusion constant (multiplied with conductivity)" << std::endl
                  << "p_1 - conductivity update rate" << std::endl
                  << "p_2 - flux feedback power" << std::endl
                  << "p_3 - control parameter (gamma)" << std::endl
                  << "p_4 - conductivity \'degradation\' rate." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) {
        std::cerr << "DiffusionConductiveSimple::"
                  << "DiffusionConductiveSimple() "
                  << "Two levels of variable indices used, "
                  << "first for diffusive cell variable, "
                  << "second for wall conductivity variable" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("DiffusionConductiveSimple");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "p_0";
    tmp[1] = "p_1";
    tmp[2] = "p_2";
    tmp[3] = "p_3";
    tmp[4] = "p_4";

    setParameterId(tmp);
}

void DiffusionConductiveSimple::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    // Cell variable updates
    //
    size_t numCells = T.numCell();
    size_t cI = variableIndex(0, 0);
    size_t CI = variableIndex(1, 0);
    assert(cI < cellData[0].size() && CI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t wallI = T.cell(i).wall(n)->index();
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }
                if (i < neighIndex) {
                    double conductance = wallData[wallI][CI];
                    double flux = parameter(0) * conductance *
                                  (cellData[i][cI] - cellData[neighIndex][cI]);
                    cellDerivs[i][cI] -= flux;
                    cellDerivs[neighIndex][cI] += flux;
                    // Update wall variables (conductance)
                    if (conductance > 0.0) {
                        flux = std::fabs(flux);
                        wallDerivs[wallI][CI] += parameter(1) *
                                                 ((std::pow(flux, parameter(2)) / std::pow(conductance, parameter(3) + 1)) - parameter(4)) * conductance;
                        // std::cerr << wallI << " " << CI << " "
                        //	      << wallData[wallI][CI] << " "
                        //	      << wallDerivs[wallI][CI] << std::endl;
                    }
                }
            }
        }
    }
}

Diffusion2d::
    Diffusion2d(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> >
                    &indValue) {
    // Do some checks on the parameters and variable indices
    //
    if (paraValue.size() != 1) {
        std::cerr << "Diffusion2d::"
                  << "Diffusion2d() "
                  << "One parameter (diffusion constant) used: p_0" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Diffusion2d::"
                  << "Diffusion2d() "
                  << "One level of one variable indices used" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("Diffusion2d");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "p_0";

    setParameterId(tmp);
}

void Diffusion2d::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);
    size_t dimension = vertexData[0].size();
    assert(aI < cellData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }

                // calculate distance between compartments

                std::vector<double> neighpos = T.cell(neighIndex).positionFromVertex();
                std::vector<double> cellpos = T.cell(i).positionFromVertex();

                double distance = 0;

                for (size_t d = 0; d < dimension; ++d)
                    distance += (neighpos[d] - cellpos[d]) * (neighpos[d] - cellpos[d]);

                distance = std::sqrt(distance);

                // std::cerr << "distance = " << distance << "\n";

                size_t v1 = T.cell(i).wall(n)->vertex1()->index();
                size_t v2 = T.cell(i).wall(n)->vertex2()->index();
                double contactLength = 0;
                for (size_t d = 0; d < dimension; ++d)
                    contactLength += (vertexData[v1][d] - vertexData[v2][d]) * (vertexData[v1][d] - vertexData[v2][d]);
                contactLength = std::sqrt(contactLength);
                double cellVolume = T.cell(i).calculateVolume(vertexData);
                cellDerivs[i][aI] -=
                    parameter(0) * contactLength * (cellData[i][aI] - cellData[neighIndex][aI]) / (cellVolume * distance);
                cellVolume = T.cell(neighIndex).calculateVolume(vertexData);
                cellDerivs[neighIndex][aI] +=
                    parameter(0) * contactLength * (cellData[i][aI] - cellData[neighIndex][aI]) / (cellVolume * distance);
            }
        }
    }
}

void Diffusion2d::
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
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);
    size_t dimension = vertexData[0].size();
    assert(aI < cellData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        size_t numWalls = T.cell(i).numWall();

        for (size_t n = 0; n < numWalls; ++n) {
            if (T.cell(i).wall(n)->cell1() != T.background() &&
                T.cell(i).wall(n)->cell2() != T.background()) {
                size_t neighIndex;
                if (T.cell(i).wall(n)->cell1()->index() == i)
                    neighIndex = T.cell(i).wall(n)->cell2()->index();
                else {
                    neighIndex = T.cell(i).wall(n)->cell1()->index();
                }

                // calculate distance between compartments

                std::vector<double> neighpos = T.cell(neighIndex).positionFromVertex();
                std::vector<double> cellpos = T.cell(i).positionFromVertex();

                double distance = 0;

                for (size_t d = 0; d < dimension; ++d)
                    distance += (neighpos[d] - cellpos[d]) * (neighpos[d] - cellpos[d]);

                distance = std::sqrt(distance);

                // std::cerr << "distance = " << distance << "\n";

                size_t v1 = T.cell(i).wall(n)->vertex1()->index();
                size_t v2 = T.cell(i).wall(n)->vertex2()->index();
                double contactLength = 0;
                for (size_t d = 0; d < dimension; ++d)
                    contactLength += (vertexData[v1][d] - vertexData[v2][d]) * (vertexData[v1][d] - vertexData[v2][d]);
                contactLength = std::sqrt(contactLength);
                double cellVolume = T.cell(i).calculateVolume(vertexData);
                cellDerivs[i][aI] -=
                    parameter(0) * contactLength * (cellData[i][aI] - cellData[neighIndex][aI]) / (cellVolume * distance);
                cellVolume = T.cell(neighIndex).calculateVolume(vertexData);
                cellDerivs[neighIndex][aI] +=
                    parameter(0) * contactLength * (cellData[i][aI] - cellData[neighIndex][aI]) / (cellVolume * distance);
                sdydtCell[i][aI] += 0;
                sdydtCell[neighIndex][aI] += 0;
            }
        }
    }
}

ActiveTransportCellEfflux::
    ActiveTransportCellEfflux(std::vector<double> &paraValue,
                              std::vector<std::vector<size_t> >
                                  &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "ActiveTransportCellEfflux::"
                  << "ActiveTransportCellEfflux() "
                  << "1 parameters used (see transport.h)\n";
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) {
        std::cerr << "ActiveTransportCellEfflux::"
                  << "ActiveTransportCellEfflux() "
                  << "One cell variable indices (auxin) and one wall variable"
                  << " indices are used (PIN)." << std::endl;
        exit(0);
    }
    // Set the variable values
    //
    setId("ActiveTransportCellEfflux");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "T";
    setParameterId(tmp);
}

void ActiveTransportCellEfflux::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);   // auxin
    size_t pwI = variableIndex(1, 0);  // pin (membrane/wall)

    assert(aI < cellData[0].size() &&
           pwI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        // Auxin transport and protein cycling
        size_t numWalls = T.cell(i).numWall();
        for (size_t k = 0; k < numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();
            if (T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
                double fac = parameter(0) * cellData[i][aI] * wallData[j][pwI];  // p_2 A_i + p_4 A_i P_ij

                cellDerivs[i][aI] -= fac;
                cellDerivs[iNeighbor][aI] += fac;
            } else if (T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
                double fac = parameter(0) * cellData[i][aI] * wallData[j][pwI + 1];

                cellDerivs[i][aI] -= fac;
                cellDerivs[iNeighbor][aI] += fac;
            }
        }
    }
}

ActiveTransportCellEffluxMM::
    ActiveTransportCellEffluxMM(std::vector<double> &paraValue,
                                std::vector<std::vector<size_t> >
                                    &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
        std::cerr << "ActiveTransportCellEfflux::"
                  << "ActiveTransportCellEfflux() "
                  << "2 parameters used (see transport.h)\n";
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) {
        std::cerr << "ActiveTransportCellEffluxMM::"
                  << "ActiveTransportCellEffluxMM() "
                  << "One cell variable indices (auxin) and one wall variable"
                  << " indices are used (PIN)." << std::endl;
        exit(0);
    }
    // Set the variable values
    //
    setId("ActiveTransportCellEffluxMM");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "T";
    tmp[1] = "K_T";
    setParameterId(tmp);
}

void ActiveTransportCellEffluxMM::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);   // auxin
    size_t pwI = variableIndex(1, 0);  // pin (membrane/wall)

    assert(aI < cellData[0].size() &&
           pwI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        // Auxin transport and protein cycling
        size_t numWalls = T.cell(i).numWall();
        for (size_t k = 0; k < numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();
            if (T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
                double fac = parameter(0) * (cellData[i][aI] / (cellData[i][aI] + parameter(1))) * wallData[j][pwI];  // p_2 A_i + p_4 A_i P_ij

                cellDerivs[i][aI] -= fac;
                cellDerivs[iNeighbor][aI] += fac;
            } else if (T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
                double fac = parameter(0) * (cellData[i][aI] / (cellData[i][aI] + parameter(1))) * wallData[j][pwI + 1];

                cellDerivs[i][aI] -= fac;
                cellDerivs[iNeighbor][aI] += fac;
            }
        }
    }
}

DiffusionActiveTransportCell::
    DiffusionActiveTransportCell(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> >
                                     &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
        std::cerr << "DiffusionActiveTransportCell::"
                  << "DiffusionActiveTransportCell() "
                  << "2 parameters used (see Documentation or transport.h)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) &&
        (indValue.size() != 3 || indValue[0].size() != 1 || indValue[1].size() != 1 || indValue[2].size() != 1)) {
        std::cerr << "ActiveTransportCellEfflux::"
                  << "ActiveTransportCellEfflux() "
                  << "One cell variable index (auxin) at first level and one wall variable"
                  << " index (PIN) at second level are always used." << std::endl;
        std::cerr << "An extra wall index for saving flux can be given in third level." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("DiffusionActiveTransportCell");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "D";
    tmp[1] = "T";
    setParameterId(tmp);
}

void DiffusionActiveTransportCell::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);   // auxin
    size_t pwI = variableIndex(1, 0);  // pin (membrane/wall)

    assert(aI < cellData[0].size() &&
           pwI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        // Auxin transport and protein cycling
        size_t numWalls = T.cell(i).numWall();
        for (size_t k = 0; k < numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();
            if (T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
                if (i < iNeighbor) {
                    double fac = (parameter(0) + parameter(1) * wallData[j][pwI]) * cellData[i][aI];
                    fac -= (parameter(0) + parameter(1) * wallData[j][pwI + 1]) * cellData[iNeighbor][aI];
                    cellDerivs[i][aI] -= fac;
                    cellDerivs[iNeighbor][aI] += fac;
                    if (numVariableIndexLevel() == 3) {  // save flux
                        if (fac >= 0.0) {
                            wallData[j][variableIndex(2, 0)] = fac;
                            wallData[j][variableIndex(2, 0) + 1] = 0.0;
                        } else {
                            wallData[j][variableIndex(2, 0)] = 0.0;
                            wallData[j][variableIndex(2, 0) + 1] = -fac;
                        }
                    }
                }
            } else if (T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
                if (i < iNeighbor) {
                    double fac = (parameter(0) + parameter(1) * wallData[j][pwI + 1]) * cellData[i][aI];
                    fac -= (parameter(0) + parameter(1) * wallData[j][pwI]) * cellData[iNeighbor][aI];
                    cellDerivs[i][aI] -= fac;
                    cellDerivs[iNeighbor][aI] += fac;
                    if (numVariableIndexLevel() == 3) {  // save flux
                        if (fac >= 0.0) {
                            wallData[j][variableIndex(2, 0) + 1] = fac;
                            wallData[j][variableIndex(2, 0)] = 0.0;
                        } else {
                            wallData[j][variableIndex(2, 0) + 1] = 0.0;
                            wallData[j][variableIndex(2, 0)] = -fac;
                        }
                    }
                }
            }
        }
    }
}
ActiveTransportWall::
    ActiveTransportWall(std::vector<double> &paraValue,
                        std::vector<std::vector<size_t> >
                            &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 5) {
        std::cerr << "ActiveTransportWall::"
                  << "ActiveTransportWall() "
                  << "5 parameters used (see documentation or transport.h)\n";
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2) {
        std::cerr << "ActiveTransportWall::"
                  << "ActiveTransportWall() "
                  << "Two cell variable indices (auxin,,  AUX ; first row) and two wall variable"
                  << " indices (paired) are used (auxin, PIN)." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("ActiveTransportWall");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "p_IAAH(in)";
    tmp[1] = "p_IAA-(in,AUX)";
    tmp[2] = "p_IAAH(out)";
    tmp[3] = "p_IAA-(out,PIN)";
    tmp[4] = "D_IAA";

    setParameterId(tmp);
}

void ActiveTransportWall::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);    // auxin
    size_t auxI = variableIndex(0, 1);  // AUX
    size_t awI = variableIndex(1, 0);   // auxin (wall)
    size_t pwI = variableIndex(1, 1);   // PIN (membrane/wall)

    assert(aI < cellData[0].size() &&
           auxI < cellData[0].size() &&
           awI < wallData[0].size() &&
           pwI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        // Auxin transport and protein cycling
        size_t numWalls = T.cell(i).numWall();
        for (size_t k = 0; k < numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();
            if (T.cell(i).wall(k)->cell1()->index() == i &&
                T.cell(i).wall(k)->cell2() != T.background()) {
                // cell-wall transport
                double fac = (parameter(2) + parameter(3) * wallData[j][pwI]) * cellData[i][aI] - (parameter(0) + parameter(1) * cellData[i][auxI]) * wallData[j][awI];

                wallDerivs[j][awI] += fac;
                cellDerivs[i][aI] -= fac;

                // wall-wall diffusion
                fac = parameter(4) * wallData[j][awI];
                wallDerivs[j][awI] -= fac;
                wallDerivs[j][awI + 1] += fac;
            } else if (T.cell(i).wall(k)->cell2()->index() == i &&
                       T.cell(i).wall(k)->cell1() != T.background()) {
                // size_t cellNeigh = T.cell(i).wall(k)->cell1()->index();
                //  cell-wall transport
                double fac = (parameter(2) + parameter(3) * wallData[j][pwI + 1]) * cellData[i][aI] - (parameter(0) + parameter(1) * cellData[i][auxI]) * wallData[j][awI + 1];

                wallDerivs[j][awI + 1] += fac;
                cellDerivs[i][aI] -= fac;

                // wall-wall diffusion
                wallDerivs[j][awI + 1] -= parameter(4) * wallData[j][awI + 1];
                wallDerivs[j][awI] += parameter(4) * wallData[j][awI + 1];
            }
            // else {
            // std::cerr << "ActiveTransportWallEfflux::derivs() Cell-wall neighborhood wrong."
            //	  << std::endl;
            // exit(-1);
            // }
        }
    }
}

InfluxActiveTransportCell::
    InfluxActiveTransportCell(std::vector<double> &paraValue,
                              std::vector<std::vector<size_t> >
                                  &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
        std::cerr << "InfluxActiveTransportCell::"
                  << "InfluxActiveTransportCell() "
                  << "2 parameters used (see Documentation or transport.h)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) &&
        (indValue.size() != 3 || indValue[0].size() != 1 || indValue[1].size() != 1 || indValue[2].size() != 1)) {
        std::cerr << "InfluxActiveTransportCell::"
                  << "InfluxActiveTransportCell() "
                  << "One cell variable index (auxin) at first level and one wall variable"
                  << "index (AUX/LAX) at second level are always used." << std::endl;
        std::cerr << "An extra wall index for saving/updating flux can be given in third level." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("InfluxActiveTransportCell");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "T";
    tmp[1] = "k0";
    setParameterId(tmp);
}

void InfluxActiveTransportCell::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t numCells = T.numCell();
    size_t aI = variableIndex(0, 0);   // auxin (cell)
    size_t awI = variableIndex(1, 0);  // aux/lax (membrane/wall)

    assert(aI < cellData[0].size() &&
           awI < wallData[0].size());

    for (size_t i = 0; i < numCells; ++i) {
        // Auxin transport and protein cycling
        size_t numWalls = T.cell(i).numWall();
        for (size_t k = 0; k < numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();
            if (T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
                if (i < iNeighbor) {
                    double facnorm = parameter(1) + wallData[j][awI] + wallData[j][awI + 1];
                    double fac = parameter(0) * (wallData[j][awI] * cellData[iNeighbor][aI] / facnorm - wallData[j][awI + 1] * cellData[i][aI] / facnorm);
                    cellDerivs[i][aI] += fac;
                    cellDerivs[iNeighbor][aI] -= fac;
                    if (numVariableIndexLevel() == 3) {  // update flux
                        if (fac >= 0.0) {
                            wallData[j][variableIndex(2, 0)] += fac;
                            wallData[j][variableIndex(2, 0) + 1] += 0.0;
                        } else {
                            wallData[j][variableIndex(2, 0)] += 0.0;
                            wallData[j][variableIndex(2, 0) + 1] += -fac;
                        }
                    }
                }
            } else if (T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
                if (i < iNeighbor) {
                    double facnorm = parameter(1) + wallData[j][awI] + wallData[j][awI + 1];
                    double fac = parameter(0) * (wallData[j][awI + 1] * cellData[iNeighbor][aI] / facnorm - wallData[j][awI] * cellData[i][aI] / facnorm);
                    cellDerivs[i][aI] += fac;
                    cellDerivs[iNeighbor][aI] -= fac;
                    if (numVariableIndexLevel() == 3) {  // update flux
                        if (fac >= 0.0) {
                            wallData[j][variableIndex(2, 0) + 1] += fac;
                            wallData[j][variableIndex(2, 0)] += 0.0;
                        } else {
                            wallData[j][variableIndex(2, 0) + 1] += 0.0;
                            wallData[j][variableIndex(2, 0)] += -fac;
                        }
                    }
                }
            }
        }
    }
}
