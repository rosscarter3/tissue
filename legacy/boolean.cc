//
// Filename     : boolean.cc
// Description  : Classes describing some boolean reactions applied as updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2020
// Revision     : $Id:$
//

#include "boolean.h"

#include "baseReaction.h"
#include "tissue.h"

namespace Boolean {

AndGate::AndGate(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "AndGate::AndGate()  parameter used "
                  << "gatetype" << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "AndGate::AndGate() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are two indices "
                  << ", and one for the output variables" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "gatetype";
    setParameterId(tmp);
}

void AndGate::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                     DataMatrix &vertexData, DataMatrix &cellDerivs,
                     DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndGate::derivsWithAbs(Tissue &T, DataMatrix &cellData,
                            DataMatrix &wallData, DataMatrix &vertexData,
                            DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                            DataMatrix &vertexDerivs, DataMatrix &sdydtCell,
                            DataMatrix &sdydtWall, DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndGate::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                     DataMatrix &vertexData, double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        size_t gate_condition = 0;

        if (cellData[cellI][cIndex_input1] == 1 &&
            cellData[cellI][cIndex_input2] == 1) {
            gate_condition = 1;
        }

        if (gate_condition == 1) {
            cellData[cellI][cIndex_output] = 1;
        } else if (gate_condition == 0 && parameter(0) == 0) {
            cellData[cellI][cIndex_output] = 0;
        }
    }
}

AndNotGate::AndNotGate(std::vector<double> &paraValue,
                       std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "AndNotGate::AndNotGate()  parameter used "
                  << "gatetype" << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "AndNotGate::AndNotGate() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are two indices "
                  << ", and one for the output variables" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "gatetype";
    setParameterId(tmp);
}

void AndNotGate::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                        DataMatrix &vertexData, DataMatrix &cellDerivs,
                        DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndNotGate::derivsWithAbs(Tissue &T, DataMatrix &cellData,
                               DataMatrix &wallData, DataMatrix &vertexData,
                               DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                               DataMatrix &vertexDerivs, DataMatrix &sdydtCell,
                               DataMatrix &sdydtWall, DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndNotGate::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                        DataMatrix &vertexData, double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        size_t gate_condition = 0;

        if (cellData[cellI][cIndex_input1] == 1 &&
            cellData[cellI][cIndex_input2] == 0) {
            gate_condition = 1;
        }

        if (gate_condition == 1) {
            cellData[cellI][cIndex_output] = 1;
        } else if (gate_condition == 0 && parameter(0) == 0) {
            cellData[cellI][cIndex_output] = 0;
        }
    }
}

AndSpecialGate::AndSpecialGate(std::vector<double> &paraValue,
                               std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr
            << "AndSpecialGate::AndSpecialGate()  does not use any parameters."
            << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 3 ||
        indValue[1].size() != 1) {
        std::cerr << "AndSpecialGate::AndSpecialGate() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are three indices "
                  << ", and one for the output variables" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
    //
}

void AndSpecialGate::derivs(Tissue &T, DataMatrix &cellData,
                            DataMatrix &wallData, DataMatrix &vertexData,
                            DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                            DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndSpecialGate::derivsWithAbs(Tissue &T, DataMatrix &cellData,
                                   DataMatrix &wallData, DataMatrix &vertexData,
                                   DataMatrix &cellDerivs,
                                   DataMatrix &wallDerivs,
                                   DataMatrix &vertexDerivs,
                                   DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                                   DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndSpecialGate::update(Tissue &T, DataMatrix &cellData,
                            DataMatrix &wallData, DataMatrix &vertexData,
                            double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_input3 = variableIndex(0, 2);
    size_t cIndex_output = variableIndex(1, 0);

    // size_t cond=0;
    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        cellData[cellI][cIndex_output] = 0;
        // if (cellData[cellI][cIndex_input2]==1 ||
        // cellData[cellI][cIndex_input2]==0){cond=1;} if
        // (cellData[cellI][cIndex_input2]==0){cond=1;}
        //    if (cellData[cellI][cIndex_input1]==1 && cond==1 &&
        //    cellData[cellI][cIndex_input3]==0)

        if (cellData[cellI][cIndex_input1] == 1 &&
            cellData[cellI][cIndex_input2] == 0 &&
            cellData[cellI][cIndex_input3] == 0) {
            cellData[cellI][cIndex_output] = 1;
        }
    }
}

AndSpecialGate2::AndSpecialGate2(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr
            << "AndSpecialGate2::AndSpecialGate2()  does not use any parameters."
            << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 3 ||
        indValue[1].size() != 1) {
        std::cerr << "AndSpecialGate2::AndSpecialGate2() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are three indices "
                  << ", and one for the output variables" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
    //
}

void AndSpecialGate2::derivs(Tissue &T, DataMatrix &cellData,
                             DataMatrix &wallData, DataMatrix &vertexData,
                             DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                             DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndSpecialGate2::derivsWithAbs(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs, DataMatrix &sdydtCell, DataMatrix &sdydtWall,
    DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndSpecialGate2::update(Tissue &T, DataMatrix &cellData,
                             DataMatrix &wallData, DataMatrix &vertexData,
                             double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_input3 = variableIndex(0, 2);
    size_t cIndex_output = variableIndex(1, 0);

    // size_t cond=0;
    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        if (cellData[cellI][cIndex_input1] == 1 &&
            cellData[cellI][cIndex_input2] == 1 &&
            cellData[cellI][cIndex_input3] == 0) {
            cellData[cellI][cIndex_output] = 1;
        }
    }
}

AndSpecialGate3::AndSpecialGate3(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "AndSpecialGate3::AndSpecialGate3()  needs one parameter."
                  << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 3 ||
        indValue[1].size() != 1) {
        std::cerr << "AndSpecialGate3::AndSpecialGate3() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are three indices "
                  << ", and one for the output variables" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
    //
}

void AndSpecialGate3::derivs(Tissue &T, DataMatrix &cellData,
                             DataMatrix &wallData, DataMatrix &vertexData,
                             DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                             DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndSpecialGate3::derivsWithAbs(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs, DataMatrix &sdydtCell, DataMatrix &sdydtWall,
    DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndSpecialGate3::update(Tissue &T, DataMatrix &cellData,
                             DataMatrix &wallData, DataMatrix &vertexData,
                             double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_input3 = variableIndex(0, 2);
    size_t cIndex_output = variableIndex(1, 0);

    // size_t cond=0;
    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        // if (cellData[cellI][cIndex_input2]==1 ||
        // cellData[cellI][cIndex_input2]==0){cond=1;} if
        // (cellData[cellI][cIndex_input2]==0){cond=1;}
        //    if (cellData[cellI][cIndex_input1]==1 && cond==1 &&
        //    cellData[cellI][cIndex_input3]==0)

        if (cellData[cellI][cIndex_input1] > parameter(0) &&
            cellData[cellI][cIndex_input2] == 1 &&
            cellData[cellI][cIndex_input3] == 0) {
            cellData[cellI][cIndex_output] = 1;
        }
        //   else
        //      {cellData[cellI][cIndex_output]=0;}

        // if (cellData[cellI][cIndex_input1]==1 &&
        // cellData[cellI][cIndex_input2]==1) {if(cellData[cellI][cIndex_input3]==0)
        // {cellData[cellI][cIndex_output]=1;}
        // else
        // {cellData[cellI][cIndex_output]=0;}
        //}
        // else
        //  {cellData[cellI][cIndex_output]=0;}
    }
}

AndGateCount::AndGateCount(std::vector<double> &paraValue,
                           std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "No parameters are used in AndGateCount" << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "AndGateCount::AndGateCount() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are two indices "
                  << ", and one for the output variable" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void AndGateCount::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                          DataMatrix &vertexData, DataMatrix &cellDerivs,
                          DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndGateCount::derivsWithAbs(Tissue &T, DataMatrix &cellData,
                                 DataMatrix &wallData, DataMatrix &vertexData,
                                 DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                                 DataMatrix &vertexDerivs,
                                 DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                                 DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndGateCount::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                          DataMatrix &vertexData, double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        if (cellData[cellI][cIndex_input1] == 1 &&
            cellData[cellI][cIndex_input2] == 1) {
            cellData[cellI][cIndex_output] += 1;
        }
    }
}

OrGateCount::OrGateCount(std::vector<double> &paraValue,
                         std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "No parameters are used in OrGateCount" << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "OrGateCount::OrGateCount() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are two indices "
                  << ", and one for the output variable" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void OrGateCount::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                         DataMatrix &vertexData, DataMatrix &cellDerivs,
                         DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void OrGateCount::derivsWithAbs(Tissue &T, DataMatrix &cellData,
                                DataMatrix &wallData, DataMatrix &vertexData,
                                DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                                DataMatrix &vertexDerivs, DataMatrix &sdydtCell,
                                DataMatrix &sdydtWall,
                                DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void OrGateCount::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                         DataMatrix &vertexData, double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        if (cellData[cellI][cIndex_input1] == 1 ||
            cellData[cellI][cIndex_input2] == 1) {
            cellData[cellI][cIndex_output] += 1;
        }
    }
}

OrSpecialGateCount::OrSpecialGateCount(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "No parameters are used in OrSpecialGateCount" << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "OrSpecialGateCount::OrSpecialGateCount() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are two indices "
                  << ", and one for the output variable" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void OrSpecialGateCount::derivs(Tissue &T, DataMatrix &cellData,
                                DataMatrix &wallData, DataMatrix &vertexData,
                                DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                                DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void OrSpecialGateCount::derivsWithAbs(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs, DataMatrix &sdydtCell, DataMatrix &sdydtWall,
    DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void OrSpecialGateCount::update(Tissue &T, DataMatrix &cellData,
                                DataMatrix &wallData, DataMatrix &vertexData,
                                double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        if (cellData[cellI][cIndex_input1] == 1 ||
            cellData[cellI][cIndex_input2] > 0) {
            cellData[cellI][cIndex_output] += 1;
        }
    }
}

AndThresholdsGate::AndThresholdsGate(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2) {
        std::cerr << "AndThresholdsGate::AndThresholdsGate() 2 parameter used "
                  << "which are the two thresholds." << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 2 ||
        indValue[1].size() != 1) {
        std::cerr << "AndThresholdsGate::AndThresholdsGate() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variables, which are two indices "
                  << ", and one for the output variables" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "gatetype";
    setParameterId(tmp);
}

void AndThresholdsGate::derivs(Tissue &T, DataMatrix &cellData,
                               DataMatrix &wallData, DataMatrix &vertexData,
                               DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                               DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void AndThresholdsGate::derivsWithAbs(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs, DataMatrix &sdydtCell, DataMatrix &sdydtWall,
    DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void AndThresholdsGate::update(Tissue &T, DataMatrix &cellData,
                               DataMatrix &wallData, DataMatrix &vertexData,
                               double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input1 = variableIndex(0, 0);
    size_t cIndex_input2 = variableIndex(0, 1);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        if (cellData[cellI][cIndex_input1] > parameter(0) &&
            cellData[cellI][cIndex_input2] > parameter(1)) {
            cellData[cellI][cIndex_output] = 1;
        }
    }
}

Count::Count(std::vector<double> &paraValue,
             std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "No parameters are used in Count" << std::endl;
        exit(0);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Count::Count() "
                  << "One level of variable index is used " << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void Count::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                   DataMatrix &vertexData, DataMatrix &cellDerivs,
                   DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void Count::derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                          DataMatrix &vertexData, DataMatrix &cellDerivs,
                          DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                          DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                          DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void Count::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                   DataMatrix &vertexData, double h) {
    size_t numCells = T.numCell();

    size_t cIndex_output = variableIndex(0, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        cellData[cellI][cIndex_output] += 1.0;
    }
}

FlagCount::FlagCount(std::vector<double> &paraValue,
                     std::vector<std::vector<size_t> > &indValue) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 0) {
        std::cerr << "No parameters are used in FlagCount" << std::endl;
        exit(0);
    }
    if (indValue.size() != 2 || indValue[0].size() != 1 ||
        indValue[1].size() != 1) {
        std::cerr << "FlagCount::FlagCount() "
                  << "Two levels of variable indices are used, "
                  << "One for the input variable, "
                  << ", and one for the output variable" << std::endl;
        exit(0);
    }
    //
    // Set the variable values
    //
    setId("add");
    setParameter(paraValue);
    setVariableIndex(indValue);
}

void FlagCount::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
    // Nothing to be done for the derivative function.
}

void FlagCount::derivsWithAbs(Tissue &T, DataMatrix &cellData,
                              DataMatrix &wallData, DataMatrix &vertexData,
                              DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                              DataMatrix &vertexDerivs, DataMatrix &sdydtCell,
                              DataMatrix &sdydtWall, DataMatrix &sdydtVertex) {
    // Nothing to be done for the derivative function.
}

void FlagCount::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, double h) {
    size_t numCells = T.numCell();

    size_t cIndex_input = variableIndex(0, 0);
    size_t cIndex_output = variableIndex(1, 0);

    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        if (cellData[cellI][cIndex_input] == 1) {
            cellData[cellI][cIndex_output] += 1.0;
        }
    }
}

}  // end namespace Boolean
