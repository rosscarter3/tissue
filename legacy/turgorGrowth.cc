//
// Filename     : turgorGrowth.cc
// Description  : Classes describing growth updates based on active updates of turgor pressure
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//
#include "turgorGrowth.h"

#include <cmath>
#include <vector>

#include "baseReaction.h"

namespace TurgorGrowth {

WaterVolume::
    WaterVolume(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 5) {
        std::cerr << "TurgorGrowth::WaterVolume::WaterVolume() "
                  << "Uses five parameters: k_p, P_max, k_pp and "
                  << "denyShrink_flag allowNegTurgor_flag." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() < 1 || indValue.size() > 2 || indValue[0].size() != 1 || (indValue.size() == 2 && indValue[1].size() != 1)) {
        std::cerr << "TurgorGrowth::WaterVolume::WaterVolume() "
                  << "Water volume index must be given in "
                  << "first level.\n"
                  << "Optionally index for saving the turgor pressure can be"
                  << " given at second level." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("TurgorGrowth::WaterVolume");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_p";
    tmp[1] = "P_max";
    tmp[2] = "k_pp";
    tmp[3] = "denyShrink_flag";
    tmp[4] = "allowNegTurgor_flag";
    setParameterId(tmp);
}

void WaterVolume::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);

        double P = 0.0;
        double totalLength = 0.0;
        for (size_t i = 0; i < cell.numWall(); ++i) {
            size_t vertex1Index = cell.wall(i)->vertex1()->index();
            size_t vertex2Index = cell.wall(i)->vertex2()->index();
            size_t dimensions = vertexData[vertex1Index].size();

            double distance = 0.0;
            for (size_t d = 0; d < dimensions; ++d) {
                distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]) * (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
            }
            distance = std::sqrt(distance);
            totalLength += distance;

            // Old turgor measure from wall extensions
            // for (size_t j = 0; j < numVariableIndex(1); ++j)
            // P += wallData[cell.wall(i)->index()][variableIndex(1, j)]/distance;
        }

        // Calculate turgor measure from volume and 'water volume'
        // P ~ p_2(V_w-V)/V
        double cellVolume = cell.calculateVolume(vertexData);
        P = (cellData[cell.index()][variableIndex(0, 0)] - cellVolume) / cellVolume;
        if (P < 0.0 && !parameter(4))
            P = 0.0;

        P *= parameter(2);

        if (numVariableIndexLevel() == 2)
            cellData[n][variableIndex(1, 0)] = P;

        if (!parameter(3) || parameter(1) - P > 0.0)
            cellDerivs[cell.index()][variableIndex(0, 0)] +=
                parameter(0) * (parameter(1) - P) * totalLength;
    }
}

}  // namespace TurgorGrowth
