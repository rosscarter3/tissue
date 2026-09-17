//
// Filename     : mechanical.cc
// Description  : Classes describing updates due to mechanical interaction
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include "mechanical.h"

#include <utility>
#include <vector>

#include "baseReaction.h"
#include "tissue.h"

namespace CenterTriangulation {
VertexFromCellPressure::VertexFromCellPressure(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 2 || (paraValue[1] != 0.0 && paraValue[1] != 1.0)) {
        std::cerr
            << "CenterTriangulation::VertexFRomCellPressure::"
            << "VertexFromCellPressure() " << std::endl
            << "Uses two parameters K_force and normalizeVolumeFlag (= 0 or 1).\n";
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 2) {
        std::cerr << "CenterTriangulation::VertexFRomCellPressure::"
                  << "VertexFromCellPressure() " << std::endl
                  << "Start of additional Cell variable indices (center(x,y,z) "
                  << "L_1,...,L_n, n=num vertex) is given in first level."
                  << "concentration index is given in second level. " << std::endl;
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("CenterTriangulation::VertexFromCellPressure");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    tmp[1] = "f_V_norm";
    setParameterId(tmp);
}

void VertexFromCellPressure::derivs(Tissue &T, DataMatrix &cellData,
                                    DataMatrix &wallData,
                                    DataMatrix &vertexData,
                                    DataMatrix &cellDerivs,
                                    DataMatrix &wallDerivs,
                                    DataMatrix &vertexDerivs) {
    // Do the update for each vertex via each wall in each cell
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    std::vector<double> cellCenter(dimension);

    // Assumming vertices and walls are sorted in 3 dimensions
    if (dimension != 3) {
        std::cerr << "CellTriangulationVertexFromCellPressure::derivs() "
                  << "assumes vertices in 3 dimensions." << std::endl;
        exit(-1);
    }
    //
    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        Cell &tmpCell = T.cell(cellI);

        for (size_t d = 0; d < dimension; ++d) {
            cellCenter[d] = cellData[cellI][variableIndex(0, 0) + d];
        }

        for (size_t k = 0; k < tmpCell.numVertex(); ++k) {
            size_t v1I = tmpCell.vertex(k)->index();
            size_t v2I = tmpCell.vertex((k + 1) % (tmpCell.numVertex()))->index();

            // //normal to the triangle plane
            // std::vector<double> nCell(dimension), nWall(dimension);
            // nCell[0]=(vertexData[v1I][1]-cellCenter[1])*(vertexData[v2I][2]-cellCenter[2])
            //   -(vertexData[v1I][2]-cellCenter[2])*(vertexData[v2I][1]-cellCenter[1]);
            // nCell[1]=(vertexData[v1I][2]-cellCenter[2])*(vertexData[v2I][0]-cellCenter[0])
            //   -(vertexData[v1I][0]-cellCenter[0])*(vertexData[v2I][2]-cellCenter[2]);
            // nCell[2]=(vertexData[v1I][0]-cellCenter[0])*(vertexData[v2I][1]-cellCenter[1])
            //   -(vertexData[v1I][1]-cellCenter[1])*(vertexData[v2I][0]-cellCenter[0]);
            // //NOTE normal to the wall outward counter clockwise/clockwise dependent
            // nWall[0]=(vertexData[v2I][1]-vertexData[v1I][1])*nCell[2]
            //   -(vertexData[v2I][2]-vertexData[v1I][2])*nCell[1];
            // nWall[1]=(vertexData[v2I][2]-vertexData[v1I][2])*nCell[0]
            //   -(vertexData[v2I][0]-vertexData[v1I][0])*nCell[2];
            // nWall[2]=(vertexData[v2I][0]-vertexData[v1I][0])*nCell[1]
            //   -(vertexData[v2I][1]-vertexData[v1I][1])*nCell[0];

            // double dxNorm = 0.0;
            // dxNorm =nWall[0]*nWall[0]+nWall[1]*nWall[1]+nWall[2]*nWall[2];
            // if (dxNorm>0.0) {
            //   dxNorm = std::sqrt(dxNorm);
            //   for( size_t d=0 ; d<dimension ; ++d ) {
            //     nWall[d] /=dxNorm;
            //   }
            // }
            // else {
            //   std::cerr << "CellTriangulationVertexFromCellPressure::derivs() "
            // 	    << "strange wall length or direction." << std::endl;
            //   exit(-1);
            // }

            std::vector<double> x0(dimension), dx(dimension);
            double dxNorm = 0.0;
            for (size_t d = 0; d < dimension; ++d) {
                x0[d] = 0.5 * (vertexData[v1I][d] + vertexData[v2I][d]);
                // Caveat: This is NOT  normal to the wall in the triangle plane
                dx[d] = x0[d] - cellCenter[d];
            }

            double wallLength = 0.0;
            double wallFactor = 0.0;  // dx.wallVector/wallVector^2 for making dx
                                      // perpendicular to the wall
            std::vector<double> wallVector(dimension);

            for (size_t d = 0; d < dimension; ++d) {
                wallVector[d] = vertexData[v1I][d] - vertexData[v2I][d];
                wallLength += wallVector[d] * wallVector[d];
            }
            if (wallLength > 0.0) {
                wallLength = std::sqrt(wallLength);
            } else {
                std::cerr << "CenterTriangulation::VertexFromCellPressure::derivs() "
                          << "Strange wall length." << std::endl;
                exit(-1);
            }
            wallFactor = (dx[0] * wallVector[0] + dx[1] * wallVector[1] +
                          dx[2] * wallVector[2]) /
                         (wallLength * wallLength);
            for (size_t d = 0; d < dimension; ++d) {
                dx[d] = dx[d] - wallFactor * wallVector[d];
                dxNorm += dx[d] * dx[d];
            }
            if (dxNorm > 0.0) {
                dxNorm = std::sqrt(dxNorm);
                double dxInv = 1.0 / dxNorm;
                for (size_t d = 0; d < dimension; ++d) {
                    dx[d] = dx[d] * dxInv;
                }
            } else {
                std::cerr << "CellTriangulationVertexFromCellPressure::derivs() "
                          << "Force direction undetermined." << std::endl;
                exit(-1);
            }

            double factor = 0.5 * parameter(0);

            if (parameter(1) == 1) {
                // NOTE maybe this one should be calculated using the central mesh
                // vertex?
                // Behruz: it is ok now
                double cellVolume = tmpCell.calculateVolume(vertexData);
                factor /= std::fabs(cellVolume);
            }

            if (variableIndex(0, 1) != 0) {
                factor *= cellData[cellI][variableIndex(0, 1)];
            }

            factor *= wallLength;
            for (size_t d = 0; d < dimension; ++d) {
                vertexDerivs[v1I][d] += factor * dx[d];
                vertexDerivs[v2I][d] += factor * dx[d];
            }
        }
    }
}

VertexFromCellPressureLinear::VertexFromCellPressureLinear(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 3 || (paraValue[1] != 0.0 && paraValue[1] != 1.0)) {
        std::cerr
            << "CenterTriangulation::VertexFromCellPressureLinear::"
            << "VertexFromCellPressureLinear() " << std::endl
            << "Uses three parameters K_force and normalizeVolumeFlag (= 0 or 1)"
            << "and deltaT that sets the time the linear increase is applied.\n";
        exit(0);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "CenterTriangulation::VertexFromCellPressureLinear::"
                  << "VertexFromCellPressureLinear() " << std::endl
                  << "Start of additional Cell variable indices (center(x,y,z) "
                  << "L_1,...,L_n, n=num vertex) is given in first level."
                  << std::endl;
        exit(0);
    }
    // Set the variable values
    //
    setId("CenterTriangulation::VertexFromCellPressureLinear");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    tmp[1] = "f_V_norm";
    tmp[2] = "deltaT";
    timeFactor_ = 0.0;
    setParameterId(tmp);
}

void VertexFromCellPressureLinear::derivs(Tissue &T, DataMatrix &cellData,
                                          DataMatrix &wallData,
                                          DataMatrix &vertexData,
                                          DataMatrix &cellDerivs,
                                          DataMatrix &wallDerivs,
                                          DataMatrix &vertexDerivs) {
    // Do the update for each vertex via each wall in each cell
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    std::vector<double> cellCenter(dimension);

    // Assumming vertices and walls are sorted
    //
    // For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
        Cell &tmpCell = T.cell(cellI);

        for (size_t d = 0; d < dimension; ++d) {
            cellCenter[d] = cellData[cellI][variableIndex(0, 0) + d];
        }

        for (size_t k = 0; k < tmpCell.numVertex(); ++k) {
            size_t v1I = tmpCell.vertex(k)->index();
            size_t v2I = tmpCell.vertex((k + 1) % (tmpCell.numVertex()))->index();
            std::vector<double> x0(dimension), dx(dimension);
            double dxNorm = 0.0;
            for (size_t d = 0; d < dimension; ++d) {
                x0[d] = 0.5 * (vertexData[v1I][d] + vertexData[v2I][d]);
                // Caveat: This is  not normal to the wall in the triangle plane
                dx[d] = x0[d] - cellCenter[d];
            }
            double wallLength = 0.0;
            double wallFactor = 0.0;  // dx.wallVector/wallVector^2 for making dx
                                      // perpendicular to the wall
            std::vector<double> wallVector(dimension);

            for (size_t d = 0; d < dimension; ++d) {
                wallVector[d] = vertexData[v1I][d] - vertexData[v2I][d];
                wallLength += wallVector[d] * wallVector[d];
            }
            if (wallLength > 0.0) {
                wallLength = std::sqrt(wallLength);
            } else {
                std::cerr
                    << "CenterTriangulation::VertexFromCellPressureLinear::derivs() "
                    << "Strange wall length." << std::endl;
                std::cerr << wallLength << " " << cellI << std::endl;

                // for(size_t gg=0; gg<cellData.size(); ++gg){
                //   for(size_t g=39; g<cellData[gg].size(); ++g)
                //     std::cerr<<cellData[gg][g]<<"  ";
                //   std::cerr<<std::endl;
                //}
                exit(0);
            }
            wallFactor = (dx[0] * wallVector[0] + dx[1] * wallVector[1] +
                          dx[2] * wallVector[2]) /
                         (wallLength * wallLength);
            for (size_t d = 0; d < dimension; ++d) {
                dx[d] =
                    dx[d] - wallFactor * wallVector[d];  // now it is normal to the wall
                dxNorm += dx[d] * dx[d];
            }
            if (dxNorm > 0.0) {
                dxNorm = std::sqrt(dxNorm);
                double dxInv = 1.0 / dxNorm;
                for (size_t d = 0; d < dimension; ++d) {
                    dx[d] = dx[d] * dxInv;
                }
            } else {
                std::cerr
                    << "CenterTriangulation::VertexFromCellPressureLinear::derivs() "
                    << "Force direction undetermined." << std::endl;
                exit(0);
            }
            double factor = 0.5 * timeFactor_ * parameter(0);
            if (parameter(1) == 1) {
                // NOTE maybe this one should be calculated using the central mesh
                // vertex?
                double cellVolume = tmpCell.calculateVolume(vertexData);
                factor /= std::fabs(cellVolume);
            }
            factor *= wallLength;

            // factor *= 1.5 - cellData[cellI][11];

            // if(totaltime<100) for residual stress calculation
            for (size_t d = 0; d < dimension; ++d) {
                vertexDerivs[v1I][d] += factor * dx[d];
                vertexDerivs[v2I][d] += factor * dx[d];
            }
        }
    }
}

void VertexFromCellPressureLinear::update(Tissue &T, DataMatrix &cellData,
                                          DataMatrix &wallData,
                                          DataMatrix &vertexData, double h) {
    static double tt = 0;
    tt += h;
    totaltime = tt;
    // AGAIN, WHY THIS???
    // if (totaltime < 200)
    if (true) {
        if (timeFactor_ < 1.0) {
            timeFactor_ += h / parameter(numParameter() - 1);
        }
        if (timeFactor_ >= 1.0) timeFactor_ = 1.0;
    } else {
        timeFactor_ = 0;
    }
    // cellData[0][12]=timeFactor_*parameter(0);
}

}  // end namespace CenterTriangulation

TargetAreaFromPressure::TargetAreaFromPressure(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 4) {
        std::cerr << "TargetAreaFromPressure::TargetAreaFromPressure() "
                  << "Uses four parameters: k_p, P_max, k_pp and "
                  << "allowShrink_flag." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() < 2 || indValue.size() > 3 || indValue[0].size() != 2 ||
        (indValue.size() == 3 && indValue[2].size() != 1)) {
        std::cerr << "TargetAreaFromPressure::TargetAreaFromPressure() "
                  << "Wall length index and cell volume index must be given in "
                  << "first level.\n"
                  << "Force indices must be given in second level. "
                  << "Optionally index for saving the pressure can be"
                  << " given at third level." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("TargetAreaFromPressure");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_p";
    tmp[1] = "P_max";
    tmp[2] = "k_pp";
    tmp[3] = "f_allowShrink";

    setParameterId(tmp);
}

void TargetAreaFromPressure::derivs(Tissue &T, DataMatrix &cellData,
                                    DataMatrix &wallData,
                                    DataMatrix &vertexData,
                                    DataMatrix &cellDerivs,
                                    DataMatrix &wallDerivs,
                                    DataMatrix &vertexDerivs) {
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);

        double P = 0.0;
        double sum = 0.0;

        // Go through all the cell walls and calculate the pressure
        for (size_t i = 0; i < cell.numWall(); ++i) {
            size_t vertex1Index = cell.wall(i)->vertex1()->index();
            size_t vertex2Index = cell.wall(i)->vertex2()->index();
            size_t dimensions = vertexData[vertex1Index].size();

            // Calculate the length of the cell wall, and sum up the cell wall
            // lengths.
            double distance = 0.0;
            for (size_t d = 0; d < dimensions; ++d) {
                distance +=
                    (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]) *
                    (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
            }
            distance = std::sqrt(distance);
            sum += distance;

            // Retrieve the different wall forces and divide them by the distance
            // between the vertices. That is: sum up the forces / wall length applied
            // on a single wall.
            for (size_t j = 0; j < numVariableIndex(1); ++j) {
                P += wallData[cell.wall(i)->index()][variableIndex(1, j)] / distance;
            }
        }
        P *= parameter(2);  // Multiply by k_pp (to 'normalise/rescale' the Forces)

        // If we have set to store the pressure in a variable, do this here.
        if (numVariableIndexLevel() == 3) {
            cellData[n][variableIndex(2, 0)] = P;
        }

        // If the pressure is smaller than the max pressure (param(1)), or if
        // we allow for shrinkage, update the volume.
        if (parameter(3) || parameter(1) - P > 0.0) {
            cellDerivs[cell.index()][variableIndex(0, 1)] +=
                parameter(0) * (parameter(1) - P) * sum;
        }
    }
}

VertexFromCellPowerdiagram::VertexFromCellPowerdiagram(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    // Do some checks on the parameters and variable indeces
    //
    if (paraValue.size() != 1) {
        std::cerr << "VertexFromCellPowerdiagram::"
                  << "VertexFromCellPowerdiagram() "
                  << "Uses one parameter K_force.\n";
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "VertexFromCellPowerdiagram::"
                  << "VertexFromCellPowerdiagram() "
                  << "Cell radius index given.\n";
        exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("VertexFromCellPowerdiagram");
    setParameter(paraValue);
    setVariableIndex(indValue);

    // Set the parameter identities
    //
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "K_force";
    setParameterId(tmp);
}

void VertexFromCellPowerdiagram::derivs(Tissue &T, DataMatrix &cellData,
                                        DataMatrix &wallData,
                                        DataMatrix &vertexData,
                                        DataMatrix &cellDerivs,
                                        DataMatrix &wallDerivs,
                                        DataMatrix &vertexDerivs) {
    // Do the update for each vertex
    // size_t numCells = T.numCell();
    size_t numVertices = T.numVertex();
    size_t cellRIndex = variableIndex(0, 0);
    size_t dimension = T.vertex(0).numPosition();

    for (size_t i = 0; i < numVertices; ++i) {
        // Only calculate if vertex connected to three cells
        size_t numCellForVertex = T.vertex(i).numCell();
        if (numCellForVertex == 3) {
            // Calculate position for all three cells
            DataMatrix cellPos(numCellForVertex);
            std::vector<double> cellR(numCellForVertex);
            for (size_t cellI = 0; cellI < numCellForVertex; ++cellI) {
                cellR[cellI] = cellData[T.vertex(i).cell(cellI)->index()][cellRIndex];
                // cellR[cellI] = T.cell( T.vertex(i).cell(cellI)->index()
                // ).variable(0);
                cellPos[cellI] = T.cell(T.vertex(i).cell(cellI)->index())
                                     .positionFromVertex(vertexData);
            }
            // std::cerr << "* " << cellPos[0][0] << " " << cellPos[0][1] << "  "
            //	<< cellPos[1][0] << " " << cellPos[1][1] << "  "
            //	<< cellPos[2][0] << " " << cellPos[2][1] << std::endl;

            // Calculate optimal position from cell positions and raddii
            // according to power diagram formulation
            if (dimension != 2) {
                std::cerr << "VertexFromCellPowerdiagram::derivs() Only defined "
                          << "for two dimensions so far..." << std::endl;
                exit(-1);
            }
            double Kji =
                cellPos[1][0] * cellPos[1][0] + cellPos[1][1] * cellPos[1][1] -
                cellPos[0][0] * cellPos[0][0] - cellPos[0][1] * cellPos[0][1] -
                (cellR[1] * cellR[1] - cellR[0] * cellR[0]);
            double Kki =
                cellPos[2][0] * cellPos[2][0] + cellPos[2][1] * cellPos[2][1] -
                cellPos[0][0] * cellPos[0][0] - cellPos[0][1] * cellPos[0][1] -
                (cellR[2] * cellR[2] - cellR[0] * cellR[0]);
            std::vector<double> powPos(dimension);
            powPos[1] =
                0.5 *
                (Kji * (cellPos[2][0] - cellPos[0][0]) -
                 Kki * (cellPos[1][0] - cellPos[0][0])) /
                ((cellPos[1][1] - cellPos[0][1]) * (cellPos[2][0] - cellPos[0][0]) -
                 (cellPos[2][1] - cellPos[0][1]) * (cellPos[1][0] - cellPos[0][0]));
            powPos[0] = 0.5 * Kji / (cellPos[1][0] - cellPos[0][0]) -
                        powPos[1] * (cellPos[1][1] - cellPos[0][1]) /
                            (cellPos[1][0] - cellPos[0][0]);
            // std::cerr << i << " " << vertexData[i][0] << " " << vertexData[i][1] <<
            // "\t"
            //	<< powPos[0] << " " << powPos[1] << "\t"
            //	<< cellPos[0][0] << " " << cellPos[0][1] << "  "
            //	<< cellPos[1][0] << " " << cellPos[1][1] << "  "
            //	<< cellPos[2][0] << " " << cellPos[2][1] << "\t"
            //	<< cellR[0] << " " << cellR[1] << " " << cellR[2] << std::endl;

            // Update vertex for each dimension
            for (size_t d = 0; d < dimension; d++)
                vertexDerivs[i][d] -= parameter(0) * (vertexData[i][d] - powPos[d]);
        }
    }
}

PerpendicularWallPressure::PerpendicularWallPressure(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 1) {
        std::cerr << "PerpendicularWallPressure::PerpendicularWallPressure() "
                  << "Uses one parameter: k_force" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "PerpendicularWallPressure::PerpendicularWallPressure() "
                  << "First level gives pressure index in cell." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("PerpendicularWallPressure");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";

    setParameterId(tmp);
}

void PerpendicularWallPressure::derivs(Tissue &T, DataMatrix &cellData,
                                       DataMatrix &wallData,
                                       DataMatrix &vertexData,
                                       DataMatrix &cellDerivs,
                                       DataMatrix &wallDerivs,
                                       DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension == 2) {
        for (size_t n = 0; n < T.numWall(); ++n) {
            Wall wall = T.wall(n);
            std::vector<Cell *> cells(2);
            cells[0] = wall.cell1();
            cells[1] = wall.cell2();
            std::vector<Vertex *> vertices(2);
            vertices[0] = wall.vertex1();
            vertices[1] = wall.vertex2();

            double nx = vertexData[vertices[1]->index()][0] -
                        vertexData[vertices[0]->index()][0];
            double ny = vertexData[vertices[1]->index()][1] -
                        vertexData[vertices[0]->index()][1];
            double A = std::sqrt(nx * nx + ny * ny);
            nx = nx / A;
            ny = ny / A;

            for (size_t i = 0; i < cells.size(); ++i) {
                Cell *cell = cells[i];

                if (cell == T.background()) {
                    continue;
                }

                double force = parameter(0) *
                               cellData[cell->index()][variableIndex(0, 0)] *
                               wall.lengthFromVertexPosition(vertexData);

                double x = 0.0;
                double y = 0.0;
                for (size_t j = 0; j < cell->numVertex(); ++j) {
                    x += vertexData[cell->vertex(j)->index()][0];
                    y += vertexData[cell->vertex(j)->index()][1];
                }
                x = (x / cell->numVertex()) - vertexData[vertices[0]->index()][0];
                y = (y / cell->numVertex()) - vertexData[vertices[0]->index()][1];

                int sign = (nx * y - ny * x) >= 0 ? 1 : -1;

                for (size_t k = 0; k < vertices.size(); ++k) {
                    Vertex *vertex = vertices[k];
                    vertexDerivs[vertex->index()][0] += ny * sign * force;
                    vertexDerivs[vertex->index()][1] += -nx * sign * force;
                }
            }
        }
    } else if (dimension == 3) {
        size_t numCell = T.numCell();
        for (size_t i = 0; i < numCell; ++i) {
            Cell *cell1 = &(T.cell(i));
            size_t c1I = cell1->index();
            std::vector<double> n1 = cell1->getNormalToPCAPlane();
            int n1Sign = cell1->vectorSignFromSort(n1, vertexData);
            std::vector<double> n = n1;
            std::vector<double> x1 = cell1->positionFromVertex(vertexData);

            size_t numNeigh = cell1->numWall();
            for (size_t k = 0; k < numNeigh; ++k) {
                Cell *cell2 = cell1->cellNeighbor(k);
                size_t c2I = cell2->index();
                // Get normal to second cell (if not background) and calculate combined
                // normal (n)
                if (cell2 != T.background() && c2I < c1I) {
                    continue;
                } else if (cell2 != T.background()) {
                    std::vector<double> n2(dimension);
                    n2 = cell2->getNormalToPCAPlane();
                    int n2Sign = cell2->vectorSignFromSort(n2, vertexData);
                    for (size_t d = 0; d < dimension; ++d)
                        n[d] = n1Sign * n1[d] + n2Sign * n2[d];
                }
                // Get wall direction
                size_t v1I = cell1->wall(k)->vertex1()->index();
                size_t v2I = cell1->wall(k)->vertex2()->index();
                std::vector<double> nw(dimension), xw(dimension);
                for (size_t d = 0; d < dimension; ++d) {
                    nw[d] = vertexData[v2I][d] - vertexData[v1I][d];
                    xw[d] = 0.5 * (vertexData[v2I][d] + vertexData[v1I][d]);
                }
                // Normalize n and nw
                double norm = 0.0, normW = 0.0;
                for (size_t d = 0; d < dimension; ++d) {
                    norm += n[d] * n[d];
                    normW += nw[d] * nw[d];
                }
                assert(norm > 0.0 && normW > 0.0);
                norm = std::sqrt(norm);
                normW = std::sqrt(normW);
                double normFac = 1.0 / norm, normWFac = 1.0 / normW;
                for (size_t d = 0; d < dimension; ++d) {
                    n[d] *= normFac;
                    nw[d] *= normWFac;
                }

                // Calculate force direction perpendicular to n and nw, plus sign
                std::vector<double> nF(dimension);
                nF[0] = n[1] * nw[2] - n[2] * nw[1];
                nF[1] = n[2] * nw[0] - n[0] * nw[2];
                nF[2] = n[0] * nw[1] - n[1] * nw[0];
                norm = 0.0;
                for (size_t d = 0; d < dimension; ++d) {
                    norm += nF[d] * nF[d];
                }
                assert(norm > 0.0);
                norm = std::sqrt(norm);
                normFac = 1.0 / norm;
                for (size_t d = 0; d < dimension; ++d) {
                    nF[d] *= normFac;
                }

                // sign from scalar product?
                double scalar = 0.0;
                for (size_t d = 0; d < dimension; ++d)
                    scalar += (xw[d] - x1[d]) * nF[d];
                int sign = 1;
                if (scalar < 0.0) sign = -1;
                double forceFactor =
                    parameter(0) * cellData[c1I][variableIndex(0, 0)] * normW;
                if (cell2 != T.background())
                    forceFactor -=
                        parameter(0) * cellData[c2I][variableIndex(0, 0)] * normW;

                // Print result (wall and nF) for debugging
                // double fac=0.5*sign;
                // std::cerr << c1I << " " << vertexData[v1I][0] << " "
                //				<< vertexData[v1I][1] << " "
                //				<< vertexData[v1I][2] << " "
                //				<< xw[0] << " " << xw[1] << " " << xw[2] <<
                // std::endl;  std::cerr << c1I << " " << vertexData[v2I][0] << " "
                //				<< vertexData[v2I][1] << " "
                //				<< vertexData[v2I][2] << " "
                //				<< xw[0]+fac*nF[0] << " " << xw[1]+fac*nF[1]
                //<< " " << xw[2]+fac*nF[2]
                //				<< std::endl << std::endl << std::endl;

                for (size_t d = 0; d < dimension; ++d) {
                    vertexDerivs[v1I][d] += forceFactor * sign * nF[d];
                    vertexDerivs[v2I][d] += forceFactor * sign * nF[d];
                }
            }
        }
    } else {
        std::cerr << "PerpendicularWallPressure::derivs() Only applicable for two "
                  << "or three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
}
namespace Pressure3D {
Constant::
    Constant(std::vector<double> &paraValue,
             std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 2) {
        std::cerr << "Pressure3D::Constant::Constant() "
                  << "Uses two parameters: k_force and areaFlag" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 0.0 && paraValue[1] != 1.0) {
        std::cerr << "Pressure3D::Constant::Constant() " << std::endl
                  << "areaFlag (p_1) must be zero (no area included) or one (area included)."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr << "Pressure3D::Constant::Constant() " << std::endl
                  << "No variable index used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::Constant");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";
    tmp[1] = "areaFlag";

    setParameterId(tmp);
}

void Constant::
    derivs(Tissue &T, DataMatrix &cellData,
           DataMatrix &wallData, DataMatrix &vertexData,
           DataMatrix &cellDerivs, DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::Constant::derivs() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    unsigned int numFlipNormal = 0;
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        unsigned int flipFlag = 0;
        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            double normFac = 1.0 / norm;
            for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
        }

        std::vector<int> scalarProdSign(cell.numVertex());
        std::vector<double> scalarProdVal(cell.numVertex());
        double scalarProdSum = 0.0;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            size_t k2 = (k + 1) % cell.numVertex();
            size_t k3 = (k + 2) % cell.numVertex();
            // Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
            assert(cell.numVertex() > 2);
            std::vector<double> nw1(dimension), nw2(dimension);
            for (size_t d = 0; d < dimension; ++d) {
                nw1[d] = vertexData[cell.vertex(k2)->index()][d] -
                         vertexData[cell.vertex(k)->index()][d];
                nw2[d] = vertexData[cell.vertex(k3)->index()][d] -
                         vertexData[cell.vertex(k2)->index()][d];
            }
            // cross product
            double scalarProd = 0.0;
            for (size_t d1 = 0; d1 < dimension; ++d1) {
                size_t d2 = (d1 + 1) % dimension;
                size_t d3 = (d1 + 2) % dimension;
                scalarProd += (nw1[d1] * nw2[d2] - nw1[d2] * nw2[d1]) * normal[d3];
            }
            scalarProdVal[k] = scalarProd;
            scalarProdSum += scalarProd;
            if (scalarProd > 0.0)
                scalarProdSign[k] = 1;
            else
                scalarProdSign[k] = -1;
        }
        // for (size_t k=1; k<cell.numVertex(); ++k)
        //   if (scalarProdSign[k]!=scalarProdSign[0]) {
        //     std::cerr << "Cell " << n << " has diverging signs on scalar
        //     product." << std::endl; break;
        //   }
        int scalarProdSignSum = 0;
        for (size_t k = 0; k < scalarProdSign.size(); ++k)
            scalarProdSignSum += scalarProdSign[k];

        if (scalarProdSignSum < 0) {
            numFlipNormal++;
            flipFlag = 1;
            for (size_t d = 0; d < dimension; ++d) {
                normal[d] = -normal[d];
            }
        } else if (scalarProdSignSum == 0) {
            // std::cerr << "Pressure3D::Constant::derivs() Warning! "
            //	  << "cell " << n
            //	  << " has no majority sign in right hand rule expression."
            //	  << std::endl;
            if (std::fabs(scalarProdSum) > 0.01) {
                if (scalarProdSum < 0.0) {
                    numFlipNormal++;
                    flipFlag = 1;
                    for (size_t d = 0; d < dimension; ++d) {
                        normal[d] = -normal[d];
                    }
                }
            } else {
                // Print all walls and exit
                std::cerr << "#Pressure3D::Constant::derivs() failed finding normal direction, "
                          << "printing the wall information." << std::endl;
                std::vector<double> center = cell.positionFromVertex(vertexData);
                for (size_t k = 0; k < cell.numWall(); ++k) {
                    std::cerr << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
                    std::cerr << std::endl
                              << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
                    std::cerr << std::endl
                              << std::endl;
                }
                std::cerr << "1 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "2 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << center[d] << " ";
                std::cerr << std::endl
                          << "2 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << center[d] + normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "3 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "4 ";
                for (size_t i = 0; i < scalarProdVal.size(); ++i)
                    std::cerr << scalarProdVal[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "5 ";
                for (size_t i = 0; i < scalarProdSign.size(); ++i)
                    std::cerr << scalarProdSign[i] << " ";
                std::cerr << std::endl
                          << std::endl;

                exit(EXIT_FAILURE);
            }
        }
        // Get the cell size
        double A = 1.0;
        if (parameter(1) == 1.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();

        double coeff = parameter(0) * A;
        // update the vertex derivatives
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            double vCoeff = coeff;
            // For boundary (only force from a single cell)
            // if (cell.vertex(k)->isBoundary(T.background())) vCoeff *= 1.5;
            for (size_t d = 0; d < dimension; ++d) {
                vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
            }
        }
        // For saving normals in direction used for test plotting
        // for (size_t d=0; d<dimension; ++d)
        // cellData[cell.index()][d] = normal[d];
    }
    // std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has
    // flipped normal." << std::endl;
}

Linear::
    Linear(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 3) {
        std::cerr << "Pressure3D::Linear::Linear() "
                  << "Uses three parameters: k_force and areaFlag and time span."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 0.0 && paraValue[1] != 1.0) {
        std::cerr << "Pressure3D::Linear::Linear() "
                  << "areaFlag must be zero (no area included) or one (area included)."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[2] <= 0.0) {
        std::cerr << "Pressure3D::Linear::Linear() "
                  << "time span (for ramping up forces) must be positive." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr << "Pressure3D::Linear::Linear() "
                  << std::endl
                  << "No variable index used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::Linear");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";
    tmp[1] = "areaFlag";
    tmp[2] = "deltaT";

    timeFactor_ = 0.0;
    setParameterId(tmp);
}

void Linear::
    derivs(Tissue &T, DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::Linear::Linear() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    unsigned int numFlipNormal = 0;
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        unsigned int flipFlag = 0;
        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            double normFac = 1.0 / norm;
            for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
        }

        std::vector<int> scalarProdSign(cell.numVertex());
        std::vector<double> scalarProdVal(cell.numVertex());
        double scalarProdSum = 0.0;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            size_t k2 = (k + 1) % cell.numVertex();
            size_t k3 = (k + 2) % cell.numVertex();
            // Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
            assert(cell.numVertex() > 2);
            std::vector<double> nw1(dimension), nw2(dimension);
            for (size_t d = 0; d < dimension; ++d) {
                nw1[d] = vertexData[cell.vertex(k2)->index()][d] -
                         vertexData[cell.vertex(k)->index()][d];
                nw2[d] = vertexData[cell.vertex(k3)->index()][d] -
                         vertexData[cell.vertex(k2)->index()][d];
            }
            // cross product
            double scalarProd = 0.0;
            for (size_t d1 = 0; d1 < dimension; ++d1) {
                size_t d2 = (d1 + 1) % dimension;
                size_t d3 = (d1 + 2) % dimension;
                scalarProd += (nw1[d1] * nw2[d2] - nw1[d2] * nw2[d1]) * normal[d3];
            }
            scalarProdVal[k] = scalarProd;
            scalarProdSum += scalarProd;
            if (scalarProd > 0.0)
                scalarProdSign[k] = 1;
            else
                scalarProdSign[k] = -1;
        }
        // for (size_t k=1; k<cell.numVertex(); ++k)
        //   if (scalarProdSign[k]!=scalarProdSign[0]) {
        //     std::cerr << "Cell " << n << " has diverging signs on scalar
        //     product." << std::endl; break;
        //   }
        int scalarProdSignSum = 0;
        for (size_t k = 0; k < scalarProdSign.size(); ++k)
            scalarProdSignSum += scalarProdSign[k];

        if (scalarProdSignSum < 0) {
            numFlipNormal++;
            flipFlag = 1;
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        } else if (scalarProdSignSum == 0) {
            // std::cerr << "Pressure3D::Linear::derivs() Warning! "
            //	  << "cell " << n
            //	  << " has no majority sign in right hand rule expression."
            //	  << std::endl;
            if (std::fabs(scalarProdSum) > 0.01) {
                if (scalarProdSum < 0.0) {
                    numFlipNormal++;
                    flipFlag = 1;
                    for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
                }
            } else {
                // Print all walls and exit.
                std::cerr << "#Pressure3D::Linear::derivs() failed finding normal direction, "
                          << "prinitng the wall information." << std::endl;

                std::vector<double> center = cell.positionFromVertex(vertexData);
                for (size_t k = 0; k < cell.numWall(); ++k) {
                    std::cerr << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
                    std::cerr << std::endl
                              << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
                    std::cerr << std::endl
                              << std::endl;
                }
                std::cerr << "1 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "2 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << center[d] << " ";
                std::cerr << std::endl
                          << "2 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << center[d] + normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "3 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "4 ";
                for (size_t i = 0; i < scalarProdVal.size(); ++i)
                    std::cerr << scalarProdVal[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "5 ";
                for (size_t i = 0; i < scalarProdSign.size(); ++i)
                    std::cerr << scalarProdSign[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Get the cell size
        double A = 1.0;
        if (parameter(1) == 1.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();

        double coeff = timeFactor_ * parameter(0) * A;
        // update the vertex derivatives
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            double vCoeff = coeff;
            if (cell.vertex(k)->isBoundary(T.background())) vCoeff *= 1.5;
            for (size_t d = 0; d < dimension; ++d) {
                vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
            }
            // std::cerr << timeFactor_ * parameter(0)  << std::endl;
        }
        // For saving normals in direction used for test plotting
        // for (size_t d=0; d<dimension; ++d)
        // cellData[cell.index()][d] = normal[d];
    }
    // std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has
    // flipped normal."
    //	      << std::endl;
}

void Linear::update(Tissue &T, DataMatrix &cellData,
                    DataMatrix &wallData,
                    DataMatrix &vertexData, double h) {
    if (timeFactor_ < 1.0) {
        timeFactor_ += h / parameter(numParameter() - 1);
    }
    if (timeFactor_ > 1.0) timeFactor_ = 1.0;
    // cellData[0][12]=timeFactor_*parameter(0);
}

Spatial::Spatial(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 5) {
        std::cerr << "Pressure3D::Spatial::Spatial() "
                  << "Uses five parameters: k_min, k_max, K_spatial, n_spatial and "
                  << "areaFlag." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[4] != 0.0 && paraValue[4] != 1.0) {
        std::cerr
            << "Pressure3D::Spatial::Spatial() "
            << "areaFlag must be zero (no area included) or one (area included)."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Pressure3D::Spatial::Spatial() "
                  << std::endl
                  << "Variable index for spatial max index is used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::Spatial");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_min";
    tmp[1] = "k_max";
    tmp[2] = "K_spatial";
    tmp[3] = "n_spatial";
    tmp[4] = "areaFlag";

    setParameterId(tmp);
    Kpow_ = std::pow(paraValue[2], paraValue[3]);
}

void Spatial::derivs(Tissue &T, DataMatrix &cellData,
                     DataMatrix &wallData,
                     DataMatrix &vertexData,
                     DataMatrix &cellDerivs,
                     DataMatrix &wallDerivs,
                     DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::Spatial::derivs() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t sI = variableIndex(0, 0);
    assert(sI < vertexData[0].size());
    double max = vertexData[0][sI];
    size_t maxI = 0;
    size_t numVertices = vertexData.size();
    for (size_t i = 1; i < numVertices; ++i)
        if (vertexData[i][sI] > max) {
            max = vertexData[i][sI];
            maxI = i;
        }
    std::vector<double> maxPos(dimension);
    for (size_t d = 0; d < dimension; ++d) maxPos[d] = vertexData[maxI][d];

    unsigned int numFlipNormal = 0;
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        unsigned int flipFlag = 0;
        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            double normFac = 1.0 / norm;
            for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
        }

        std::vector<int> scalarProdSign(cell.numVertex());
        std::vector<double> scalarProdVal(cell.numVertex());
        double scalarProdSum = 0.0;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            size_t k2 = (k + 1) % cell.numVertex();
            size_t k3 = (k + 2) % cell.numVertex();
            // Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
            assert(cell.numVertex() > 2);
            std::vector<double> nw1(dimension), nw2(dimension);
            for (size_t d = 0; d < dimension; ++d) {
                nw1[d] = vertexData[cell.vertex(k2)->index()][d] -
                         vertexData[cell.vertex(k)->index()][d];
                nw2[d] = vertexData[cell.vertex(k3)->index()][d] -
                         vertexData[cell.vertex(k2)->index()][d];
            }
            // cross product
            double scalarProd = 0.0;
            for (size_t d1 = 0; d1 < dimension; ++d1) {
                size_t d2 = (d1 + 1) % dimension;
                size_t d3 = (d1 + 2) % dimension;
                scalarProd += (nw1[d1] * nw2[d2] - nw1[d2] * nw2[d1]) * normal[d3];
            }
            scalarProdVal[k] = scalarProd;
            scalarProdSum += scalarProd;
            if (scalarProd > 0.0)
                scalarProdSign[k] = 1;
            else
                scalarProdSign[k] = -1;
        }
        //  		for (size_t k=1; k<cell.numVertex(); ++k)
        //  			if (scalarProdSign[k]!=scalarProdSign[0]) {
        //  				std::cerr << "Cell " << n << " has
        //  diverging signs on scalar product." << std::endl;
        //  				break;
        //  			}
        int scalarProdSignSum = 0;
        for (size_t k = 0; k < scalarProdSign.size(); ++k)
            scalarProdSignSum += scalarProdSign[k];

        if (scalarProdSignSum < 0) {
            numFlipNormal++;
            flipFlag = 1;
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        } else if (scalarProdSignSum == 0) {
            // std::cerr << "Pressure3D::Spatial::derivs(): Warning, Cell " << n
            //	  << " has no majority sign in right hand rule expression."
            //	  << std::endl;
            if (std::fabs(scalarProdSum) > 0.01) {
                if (scalarProdSum < 0.0) {
                    numFlipNormal++;
                    flipFlag = 1;
                    for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
                }
            } else {
                std::vector<double> center = cell.positionFromVertex(vertexData);
                std::cerr << "Pressure3D::Spatial::derivs() ERROR Cell " << n
                          << " force direction can not be established. Printing all walls:"
                          << std::endl;
                for (size_t k = 0; k < cell.numWall(); ++k) {
                    std::cerr << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
                    std::cerr << std::endl
                              << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
                    std::cerr << std::endl
                              << std::endl;
                }
                std::cerr << "1 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "2 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << center[d] << " ";
                std::cerr << std::endl
                          << "2 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << center[d] + normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "3 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "4 ";
                for (size_t i = 0; i < scalarProdVal.size(); ++i)
                    std::cerr << scalarProdVal[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "5 ";
                for (size_t i = 0; i < scalarProdSign.size(); ++i)
                    std::cerr << scalarProdSign[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Get the cell size
        double A = 1.0;
        if (parameter(4) == 1.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();
        // Calculate the spatial factor

        double maxDistance = 0.0;
        std::vector<double> cellPos = cell.positionFromVertex(vertexData);
        for (size_t d = 0; d < dimension; ++d)
            maxDistance += (maxPos[d] - cellPos[d]) * (maxPos[d] - cellPos[d]);
        maxDistance = std::sqrt(maxDistance);
        double sFactor = std::pow(maxDistance, parameter(3));
        sFactor = parameter(1) * Kpow_ / (Kpow_ + sFactor);

        double coeff = (parameter(0) + sFactor) * A;

        // update the vertex derivatives
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            double vCoeff = coeff;
            if (cell.vertex(k)->isBoundary(T.background())) vCoeff *= 1.5;
            for (size_t d = 0; d < dimension; ++d) {
                vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
            }
        }
    }
    // std::cerr << numFlipNormal << " cells out of " << T.numCell()
    // << " has flipped normal." << std::endl;
}

ConcentrationHill::ConcentrationHill(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 5) {
        std::cerr << "Pressure3D::ConcentrationHill::"
                  << "ConcentrationHill() "
                  << "Uses five parameters: k_min, k_add, K_H, n_H and areaFlag."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[4] != 0.0 && paraValue[4] != 1.0) {
        std::cerr
            << "Pressure3D::ConcentrationHill::"
            << "ConcentrationHill() "
            << "areaFlag must be zero (no area included) or one (area included)."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Pressure3D::ConcentrationHill::ConcentrationHill() "
                  << std::endl
                  << "Variable index for concentration needed." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::ConcentrationHill");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_min";
    tmp[1] = "k_max";
    tmp[2] = "K_H";
    tmp[3] = "n_H";
    tmp[4] = "areaFlag";

    setParameterId(tmp);
    Kpow_ = std::pow(paraValue[2], paraValue[3]);
}

void ConcentrationHill::derivs(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::ConcentrationHill::"
                  << "derivs() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t cI = variableIndex(0, 0);
    assert(cI < cellData[0].size());

    unsigned int numFlipNormal = 0;
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        unsigned int flipFlag = 0;
        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            double normFac = 1.0 / norm;
            for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
        }

        std::vector<int> scalarProdSign(cell.numVertex());
        std::vector<double> scalarProdVal(cell.numVertex());
        double scalarProdSum = 0.0;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            size_t k2 = (k + 1) % cell.numVertex();
            size_t k3 = (k + 2) % cell.numVertex();
            // Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
            assert(cell.numVertex() > 2);
            std::vector<double> nw1(dimension), nw2(dimension);
            for (size_t d = 0; d < dimension; ++d) {
                nw1[d] = vertexData[cell.vertex(k2)->index()][d] -
                         vertexData[cell.vertex(k)->index()][d];
                nw2[d] = vertexData[cell.vertex(k3)->index()][d] -
                         vertexData[cell.vertex(k2)->index()][d];
            }
            // cross product
            double scalarProd = 0.0;
            for (size_t d1 = 0; d1 < dimension; ++d1) {
                size_t d2 = (d1 + 1) % dimension;
                size_t d3 = (d1 + 2) % dimension;
                scalarProd += (nw1[d1] * nw2[d2] - nw1[d2] * nw2[d1]) * normal[d3];
            }
            scalarProdVal[k] = scalarProd;
            scalarProdSum += scalarProd;
            if (scalarProd > 0.0)
                scalarProdSign[k] = 1;
            else
                scalarProdSign[k] = -1;
        }
        // for (size_t k=1; k<cell.numVertex(); ++k)
        //   if (scalarProdSign[k]!=scalarProdSign[0]) {
        //     std::cerr << "Cell " << n << " has diverging signs on scalar
        //     product." << std::endl; break;
        //   }
        int scalarProdSignSum = 0;
        for (size_t k = 0; k < scalarProdSign.size(); ++k)
            scalarProdSignSum += scalarProdSign[k];

        if (scalarProdSignSum < 0) {
            numFlipNormal++;
            flipFlag = 1;
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        } else if (scalarProdSignSum == 0) {
            // std::cerr << "Pressure3D::ConcentrationHill::derivs() Warning, Cell " << n
            //	  << " has no majority sign in right hand rule expression."
            //	  << std::endl;
            if (std::fabs(scalarProdSum) > 0.01) {
                if (scalarProdSum < 0.0) {
                    numFlipNormal++;
                    flipFlag = 1;
                    for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
                }
            } else {
                std::vector<double> center = cell.positionFromVertex(vertexData);
                std::cerr << "Pressure3D::ConcentrationHill::derivs() ERROR Cell " << n
                          << " force direction can not be established. Printing all walls:"
                          << std::endl;
                for (size_t k = 0; k < cell.numWall(); ++k) {
                    std::cerr << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
                    std::cerr << std::endl
                              << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
                    std::cerr << std::endl
                              << std::endl;
                }
                std::cerr << "1 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "2 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << center[d] << " ";
                std::cerr << std::endl
                          << "2 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << center[d] + normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "3 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "4 ";
                for (size_t i = 0; i < scalarProdVal.size(); ++i)
                    std::cerr << scalarProdVal[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "5 ";
                for (size_t i = 0; i < scalarProdSign.size(); ++i)
                    std::cerr << scalarProdSign[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Get the cell size
        double A = 1.0;
        if (parameter(4) == 1.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();
        // Calculate the concentration dependent factor
        double cFactor = std::pow(cellData[n][cI], parameter(3));
        cFactor = parameter(1) * cFactor / (Kpow_ + cFactor);
        double coeff = (parameter(0) + cFactor) * A;

        // update the vertex derivatives
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            double vCoeff = coeff;
            if (cell.vertex(k)->isBoundary(T.background())) vCoeff *= 1.5;
            for (size_t d = 0; d < dimension; ++d) {
                vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
            }
        }
    }
    // std::cerr << numFlipNormal << " cells out of " << T.numCell()
    // << " has flipped normal." << std::endl;
}

Normalized::Normalized(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 2) {
        std::cerr
            << "Pressure3D::Normalized::Normalized() "
            << "Uses two parameters: k_force and areaFlag" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 0.0 && paraValue[1] != 1.0) {
        std::cerr
            << "Pressure3D::Normalized::Normalized() "
            << "areaFlag (p1) needs to be zero or one (1=area factor included)."
            << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() != 0) {
        std::cerr
            << "Pressure3D::Normalized::Normalized() "
            << std::endl
            << "No variable index used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::Normalized");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";
    tmp[1] = "areaFlag";

    setParameterId(tmp);
}

void Normalized::derivs(Tissue &T, DataMatrix &cellData,
                        DataMatrix &wallData,
                        DataMatrix &vertexData,
                        DataMatrix &cellDerivs,
                        DataMatrix &wallDerivs,
                        DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr
            << "Pressure3D::Normalized::derivs() "
            << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<double> tmpD(dimension);
    DataMatrix vertexDerivsTmp(vertexDerivs.size(), tmpD);
    unsigned int numFlipNormal = 0;
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        // This calculation should now be done in reaction CalculatePCAPlane
        // cell.calculatePCAPlane(vertexData);
        unsigned int flipFlag = 0;

        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            double normFac = 1.0 / norm;
            for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
        }

        std::vector<int> scalarProdSign(cell.numVertex());
        std::vector<double> scalarProdVal(cell.numVertex());
        double scalarProdSum = 0.0;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            size_t k2 = (k + 1) % cell.numVertex();
            size_t k3 = (k + 2) % cell.numVertex();
            // Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
            assert(cell.numVertex() > 2);
            std::vector<double> nw1(dimension), nw2(dimension);
            for (size_t d = 0; d < dimension; ++d) {
                nw1[d] = vertexData[cell.vertex(k2)->index()][d] -
                         vertexData[cell.vertex(k)->index()][d];
                nw2[d] = vertexData[cell.vertex(k3)->index()][d] -
                         vertexData[cell.vertex(k2)->index()][d];
            }
            // cross product
            double scalarProd = 0.0;
            for (size_t d1 = 0; d1 < dimension; ++d1) {
                size_t d2 = (d1 + 1) % dimension;
                size_t d3 = (d1 + 2) % dimension;
                scalarProd += (nw1[d1] * nw2[d2] - nw1[d2] * nw2[d1]) * normal[d3];
            }
            scalarProdVal[k] = scalarProd;
            scalarProdSum += scalarProd;
            if (scalarProd > 0.0)
                scalarProdSign[k] = 1;
            else
                scalarProdSign[k] = -1;
        }
        // for (size_t k=1; k<cell.numVertex(); ++k)
        //   if (scalarProdSign[k]!=scalarProdSign[0]) {
        //     std::cerr << "Cell " << n << " has diverging signs on scalar
        //     product." << std::endl; break;
        //   }
        int scalarProdSignSum = 0;
        for (size_t k = 0; k < scalarProdSign.size(); ++k)
            scalarProdSignSum += scalarProdSign[k];

        if (scalarProdSignSum < 0) {
            numFlipNormal++;
            flipFlag = 1;
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        } else if (scalarProdSignSum == 0) {
            // std::cerr << "Pressure3D::Normalized::derivs() Warning, Cell " << n
            //	  << " has no majority sign in right hand rule expression."
            //	  << std::endl;
            if (std::fabs(scalarProdSum) > 0.01) {
                if (scalarProdSum < 0.0) {
                    numFlipNormal++;
                    flipFlag = 1;
                    for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
                }
            } else {
                std::vector<double> center = cell.positionFromVertex(vertexData);
                std::cerr << "Pressure3D::Normalized::derivs() ERROR Cell " << n
                          << " force direction can not be established. Printing all walls:"
                          << std::endl;
                for (size_t k = 0; k < cell.numWall(); ++k) {
                    std::cerr << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
                    std::cerr << std::endl
                              << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
                    std::cerr << std::endl
                              << std::endl;
                }
                std::cerr << "1 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "2 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << center[d] << " ";
                std::cerr << std::endl
                          << "2 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << center[d] + normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "3 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "4 ";
                for (size_t i = 0; i < scalarProdVal.size(); ++i)
                    std::cerr << scalarProdVal[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "5 ";
                for (size_t i = 0; i < scalarProdSign.size(); ++i)
                    std::cerr << scalarProdSign[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Get the cell size
        double A = 1.0;
        if (parameter(1) != 0.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();

        // update the vertex derivatives
        for (size_t d = 0; d < dimension; ++d) {
            // escellData[n][d]=normal[d];
            for (size_t k = 0; k < cell.numVertex(); ++k)
                vertexDerivsTmp[cell.vertex(k)->index()][d] += A * normal[d];
        }
    }
    // Normalize and update all vertices
    size_t nVertex = vertexDerivs.size();
    for (size_t i = 0; i < nVertex; ++i) {
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d)
            norm += vertexDerivsTmp[i][d] * vertexDerivsTmp[i][d];
        double normFac = 1.0 / std::sqrt(norm);
        for (size_t d = 0; d < dimension; ++d) {
            vertexDerivsTmp[i][d] *= normFac;
            vertexDerivs[i][d] += parameter(0) * vertexDerivsTmp[i][d];
        }
    }
    // std::cerr << numFlipNormal << " cells out of " << T.numCell()
    // << " has flipped normal." << std::endl;
}

NormalizedSpatial::NormalizedSpatial(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 5) {
        std::cerr
            << "Pressure3D::NormalizedSpatial::NormalizedSpatial() "
            << "Uses four parameters: F_min F_max K_spatial n_spatial areaFlag"
            << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[4] != 0.0 && paraValue[4] != 1.0) {
        std::cerr << "Pressure3D::NormalizedSpatial::NormalizedSpatial() "
                  << "areaFlag needs to be zero or one (1=area factor included)."
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Pressure3D::NormalizedSpatial::NormalizedSpatial() "
                  << std::endl
                  << "Variable index for spatial coordinate used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::NormalizedSpatial");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "F_min";
    tmp[1] = "F_max";
    tmp[2] = "K_spatial";
    tmp[3] = "n_spatial";
    tmp[4] = "areaFlag";

    setParameterId(tmp);
    Kpow_ = std::pow(paraValue[2], paraValue[3]);
}

void NormalizedSpatial::
    derivs(
        Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
        DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::NormalizedSpatial::derivs() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    // Initiate spatial factor
    size_t sI = variableIndex(0, 0);
    size_t numVertices = vertexData.size();
    assert(sI < dimension);
    double sMax = vertexData[0][sI];
    double maxI = 0;
    for (size_t i = 1; i < numVertices; ++i)
        if (vertexData[i][sI] > sMax) {
            sMax = vertexData[i][sI];
            maxI = i;
        }

    std::vector<double> maxPos(dimension);
    for (size_t d = 0; d < dimension; ++d) maxPos[d] = vertexData[maxI][d];

    std::vector<double> tmpD(dimension);
    DataMatrix vertexDerivsTmp(vertexDerivs.size(), tmpD);
    unsigned int numFlipNormal = 0;
    size_t numCells = T.numCell();
    for (size_t n = 0; n < numCells; ++n) {
        Cell cell = T.cell(n);
        // This calculation should now be done in reaction CalculatePCAPlane
        // cell.calculatePCAPlane(vertexData);
        unsigned int flipFlag = 0;

        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            double normFac = 1.0 / norm;
            for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
        }

        std::vector<int> scalarProdSign(cell.numVertex());
        std::vector<double> scalarProdVal(cell.numVertex());
        double scalarProdSum = 0.0;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            size_t k2 = (k + 1) % cell.numVertex();
            size_t k3 = (k + 2) % cell.numVertex();
            // Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
            assert(cell.numVertex() > 2);
            std::vector<double> nw1(dimension), nw2(dimension);
            for (size_t d = 0; d < dimension; ++d) {
                nw1[d] = vertexData[cell.vertex(k2)->index()][d] -
                         vertexData[cell.vertex(k)->index()][d];
                nw2[d] = vertexData[cell.vertex(k3)->index()][d] -
                         vertexData[cell.vertex(k2)->index()][d];
            }
            // cross product
            double scalarProd = 0.0;
            for (size_t d1 = 0; d1 < dimension; ++d1) {
                size_t d2 = (d1 + 1) % dimension;
                size_t d3 = (d1 + 2) % dimension;
                scalarProd += (nw1[d1] * nw2[d2] - nw1[d2] * nw2[d1]) * normal[d3];
            }
            scalarProdVal[k] = scalarProd;
            scalarProdSum += scalarProd;
            if (scalarProd > 0.0)
                scalarProdSign[k] = 1;
            else
                scalarProdSign[k] = -1;
        }
        //  		for (size_t k=1; k<cell.numVertex(); ++k)
        //  			if (scalarProdSign[k]!=scalarProdSign[0]) {
        //  				std::cerr << "Cell " << n << " has
        //  diverging signs on scalar product." << std::endl;
        //  				break;
        //  			}
        int scalarProdSignSum = 0;
        for (size_t k = 0; k < scalarProdSign.size(); ++k)
            scalarProdSignSum += scalarProdSign[k];

        if (scalarProdSignSum < 0) {
            numFlipNormal++;
            flipFlag = 1;
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        } else if (scalarProdSignSum == 0) {
            // std::cerr << "Pressure3D::NormalizedSpatial::derivs() Warning, Cell " << n
            // << " has no majority sign in right hand rule expression." << std::endl;
            if (std::fabs(scalarProdSum) > 0.01) {
                if (scalarProdSum < 0.0) {
                    numFlipNormal++;
                    flipFlag = 1;
                    for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
                }
            } else {
                std::vector<double> center = cell.positionFromVertex(vertexData);
                std::cerr << "Pressure3D::NormalizedSpatial::derivs() ERROR Cell " << n
                          << " force direction can not be established. Printing all walls:"
                          << std::endl;
                for (size_t k = 0; k < cell.numWall(); ++k) {
                    std::cerr << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
                    std::cerr << std::endl
                              << "0 ";
                    for (size_t d = 0; d < dimension; ++d)
                        std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
                    std::cerr << std::endl
                              << std::endl;
                }
                std::cerr << "1 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "2 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << center[d] << " ";
                std::cerr << std::endl
                          << "2 ";
                for (size_t d = 0; d < dimension; ++d)
                    std::cerr << center[d] + normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "3 ";
                for (size_t d = 0; d < dimension; ++d) std::cerr << normal[d] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "4 ";
                for (size_t i = 0; i < scalarProdVal.size(); ++i)
                    std::cerr << scalarProdVal[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                std::cerr << "5 ";
                for (size_t i = 0; i < scalarProdSign.size(); ++i)
                    std::cerr << scalarProdSign[i] << " ";
                std::cerr << std::endl
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Get the cell size
        double A = 1.0;
        if (parameter(4) != 0.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();

        // update the vertex derivatives
        for (size_t d = 0; d < dimension; ++d) {
            // escellData[n][d]=normal[d];
            for (size_t k = 0; k < cell.numVertex(); ++k)
                vertexDerivsTmp[cell.vertex(k)->index()][d] += A * normal[d];
        }
    }
    // Normalize and update all vertices
    size_t nVertex = vertexDerivs.size();
    for (size_t i = 0; i < nVertex; ++i) {
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d)
            norm += vertexDerivsTmp[i][d] * vertexDerivsTmp[i][d];
        double normFac = 1.0 / std::sqrt(norm);
        // Calculate spatial factor
        double maxDistance = 0.0;
        for (size_t d = 0; d < dimension; ++d)
            maxDistance +=
                (maxPos[d] - vertexData[i][d]) * (maxPos[d] - vertexData[i][d]);
        maxDistance = std::sqrt(maxDistance);
        double sFactor = std::pow(maxDistance, parameter(3));
        sFactor = parameter(1) * Kpow_ / (Kpow_ + sFactor);
        for (size_t d = 0; d < dimension; ++d) {
            vertexDerivsTmp[i][d] *= normFac;
            vertexDerivs[i][d] += (parameter(0) + sFactor) * vertexDerivsTmp[i][d];
        }
    }
    // std::cerr << numFlipNormal << " cells out of " << T.numCell()
    // << " has flipped normal." << std::endl;
}

SphereCylinder::SphereCylinder(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 2) {
        std::cerr << "Pressure3D::SphereCylinder::"
                  << "SphereCylinder() "
                  << "Uses two parameters: k_force and areaFlag" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 0.0 && paraValue[1] != 1.0) {
        std::cerr << "Pressure3D::SphereCylinder::"
                  << "SphereCylinder() "
                  << "areaFlag (p1) needs to be zero or one (1=area factor included)."
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() != 0) {
        std::cerr << "Pressure3D::SphereCylinder::"
                  << "SphereCylinder() "
                  << std::endl
                  << "No variable index used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::SphereCylinder");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";
    tmp[1] = "areaFlag";

    setParameterId(tmp);
}

void SphereCylinder::derivs(Tissue &T, DataMatrix &cellData,
                            DataMatrix &wallData,
                            DataMatrix &vertexData,
                            DataMatrix &cellDerivs,
                            DataMatrix &wallDerivs,
                            DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::SphereCylinder::"
                  << "derivs() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        // This calculation should now be done in CalculatePCAPlane reaction.
        // cell.calculatePCAPlane(vertexData);
        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            for (size_t d = 0; d < dimension; ++d) normal[d] /= norm;
        }
        std::vector<double> cellPos = cell.positionFromVertex(vertexData);

        std::vector<double> scNorm(dimension);
        for (size_t d = 0; d < dimension; ++d) scNorm[d] = cellPos[d];
        if (cellPos[2] < 0.0) scNorm[2] = 0.0;
        double scalarProd = 0.0;
        for (size_t d = 0; d < dimension; ++d) scalarProd += scNorm[d] * normal[d];
        if (scalarProd < 0.0) {
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        }
        // Introduce cell size
        double A = 1.0;
        if (parameter(1) != 0.0)
            A = cell.calculateVolume(vertexData) / cell.numVertex();

        double coeff = parameter(0) * A;
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            double vCoeff = coeff;
            if (cell.vertex(k)->isBoundary(T.background())) vCoeff *= 1.5;
            for (size_t d = 0; d < dimension; ++d)
                vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
        }
    }
}

SphereCylinderConcentrationHill::
    SphereCylinderConcentrationHill(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 5) {
        std::cerr << "Pressure3D::SphereCylinderConcentrationHill::"
                  << "SphereCylinderConcentrationHill() "
                  << "Uses five parameters: k_forceConst, k_forceHill, K_Hill, "
                  << "n_Hill areaFlag"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[4] != 0.0 && paraValue[4] != 1.0) {
        std::cerr
            << "Pressure3D::SphereCylinderConcentrationHill::"
            << "SphereCylinderConcentrationHill() "
            << "areaFlag (p4) needs to be zero or one (1=area factor included)."
            << std::endl;
        exit(EXIT_FAILURE);
    }

    if (indValue.size() != 1 || indValue[0].size() != 1) {
        std::cerr << "Pressure3D::SphereCylinderConcentrationHill::"
                  << "SphereCylinderConcentrationHill() "
                  << std::endl
                  << "One variable index for concentration index is used."
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::SphereCylinderConcentrationHill");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_forceConst";
    tmp[1] = "k_forceHill";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "areaFlag";

    setParameterId(tmp);
}

void SphereCylinderConcentrationHill::
    derivs(
        Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
        DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::SphereCylinderConcentrationHill::"
                  << "derivs() "
                  << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }
    double Kpow = std::pow(parameter(2), parameter(3));
    for (size_t n = 0; n < T.numCell(); ++n) {
        Cell cell = T.cell(n);
        // This calculation should now be done in reaction CalculatePCAPlane
        // cell.calculatePCAPlane(vertexData);
        double concpow =
            std::pow(cellData[cell.index()][variableIndex(0, 0)], parameter(3));
        std::vector<double> normal = cell.getNormalToPCAPlane();
        double norm = 0.0;
        for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
        if (norm != 1.0) {
            norm = std::sqrt(norm);
            assert(norm > 0.0);
            for (size_t d = 0; d < dimension; ++d) normal[d] /= norm;
        }
        std::vector<double> cellPos = cell.positionFromVertex(vertexData);

        std::vector<double> scNorm(dimension);
        for (size_t d = 0; d < dimension; ++d) scNorm[d] = cellPos[d];
        if (cellPos[2] < 0.0) scNorm[2] = 0.0;
        double scalarProd = 0.0;
        for (size_t d = 0; d < dimension; ++d) scalarProd += scNorm[d] * normal[d];
        if (scalarProd < 0.0) {
            for (size_t d = 0; d < dimension; ++d) normal[d] = -normal[d];
        }
        double coeff = parameter(0) + parameter(1) * concpow / (concpow + Kpow);
        double A = 1.0;
        if (parameter(4) != 0.0) {
            A = cell.calculateVolume(vertexData) / cell.numVertex();
            coeff *= A;
        }
        for (size_t k = 0; k < cell.numVertex(); ++k) {
            double vCoeff = coeff;
            if (cell.vertex(k)->isBoundary(T.background())) vCoeff *= 1.5;
            for (size_t d = 0; d < dimension; ++d)
                vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
        }
    }
}

Triangular::Triangular(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
    if (paraValue.size() != 3 && paraValue.size() != 6) {
        std::cerr
            << "Pressure3D::Triangular::Triangular()"
            << " Uses three parameters: k_force, areaFlag(0/1 without/with "
            << "area-constant pressure) and Zplane" << std::endl
            << " or " << std::endl
            << "six parameters: k_force, areaFlag (2/3 without/with area-pressure"
            << " increase with area decrease), Zplane, equilibrium volume,volume"
            << " decrease factor and pressure increase factor." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 0.0 && paraValue[1] != 1.0 && paraValue[1] != 2.0 &&
        paraValue[1] != 3.0) {
        std::cerr
            << "Pressure3D::Triangular::Triangular() "
            << "areaFlag must be zero/two (no area included) or one/three (area "
            << "included)."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 0) {
        std::cerr
            << "Pressure3D::Triangular::Triangular() "
            << std::endl
            << "No variable index used." << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::Triangular");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";
    tmp[1] = "areaFlag";
    tmp[2] = "Zplane";

    if (paraValue.size() == 6) {
        tmp[3] = "V0";
        tmp[4] = "Vfactor";
        tmp[5] = "Pfactor";
    }
    setParameterId(tmp);
}

void Triangular::derivs(Tissue &T, DataMatrix &cellData,
                        DataMatrix &wallData,
                        DataMatrix &vertexData,
                        DataMatrix &cellDerivs,
                        DataMatrix &wallDerivs,
                        DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr
            << "Pressure3D::Triangular::derivs() "
            << "Only implemented for three dimensions." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (parameter(1) == 2 || parameter(1) == 3) {
        double totalVolume = 0;
        size_t numCell = cellData.size();
        assert(numCell == T.numCell());

        for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
            Cell cell = T.cell(cellIndex);
            if (cell.numVertex() != 3) {
                std::cerr
                    << "Pressure3D::Triangular::Triangular() "
                    << "Only implemented for triangular cells." << std::endl;
                exit(EXIT_FAILURE);
            }
            size_t v1 = T.cell(cellIndex).vertex(0)->index();
            size_t v2 = T.cell(cellIndex).vertex(1)->index();
            size_t v3 = T.cell(cellIndex).vertex(2)->index();
            DataMatrix position(3, vertexData[v1]);
            position[1] = vertexData[v2];
            position[2] = vertexData[v3];

            std::vector<double> normal = cell.getNormalTriangular(vertexData);
            double norm = 0.0;
            for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
            if (norm != 1.0) {
                norm = std::sqrt(norm);
                assert(norm > 0.0);
                double normFac = 1.0 / norm;
                for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
            }
            // double A= cell.calculateVolume(vertexData)*normal[2];
            // totalVolume +=
            // A*(((position[0][2]+position[1][2]+position[2][2])/3)-parameter(2));
        }

        for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
            Cell cell = T.cell(cellIndex);

            std::vector<double> normal = cell.getNormalTriangular(vertexData);
            double norm = 0.0;
            for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
            if (norm != 1.0) {
                norm = std::sqrt(norm);
                assert(norm > 0.0);
                double normFac = 1.0 / norm;
                for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
            }

            // Get the cell size
            double A = 1.0;
            if (parameter(1) == 3.0)
                A = cell.calculateVolume(vertexData) / cell.numVertex();

            double coeff =
                ((parameter(5) - 1) * parameter(0) * totalVolume /
                     ((parameter(4) - 1) * parameter(3)) +
                 (parameter(4) - parameter(5)) * parameter(0) / (parameter(4) - 1)) *
                A;
            // update the vertex derivatives
            for (size_t k = 0; k < cell.numVertex(); ++k) {
                double vCoeff = coeff;
                // if (cell.vertex(k)->isBoundary(T.background()))
                // vCoeff *= 1.5;
                for (size_t d = 0; d < dimension; ++d) {
                    vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
                }
            }
        }
    }

    if (parameter(1) == 0 || parameter(1) == 1) {
        // double totalVolume=0;
        size_t numCell = cellData.size();
        assert(numCell == T.numCell());
        for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
            Cell cell = T.cell(cellIndex);
            if (cell.numVertex() != 3) {
                std::cerr
                    << "Pressure3D::Triangular::Triangular() "
                    << "Only implemented for triangular cells." << std::endl;
                exit(EXIT_FAILURE);
            }
            size_t v1 = T.cell(cellIndex).vertex(0)->index();
            size_t v2 = T.cell(cellIndex).vertex(1)->index();
            size_t v3 = T.cell(cellIndex).vertex(2)->index();
            DataMatrix position(3, vertexData[v1]);
            position[1] = vertexData[v2];
            position[2] = vertexData[v3];

            std::vector<double> normal = cell.getNormalTriangular(vertexData);
            double norm = 0.0;
            for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
            if (norm != 1.0) {
                norm = std::sqrt(norm);
                assert(norm > 0.0);
                double normFac = 1.0 / norm;
                for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
            }
            // double A= cell.calculateVolume(vertexData)*normal[2];
            // totalVolume +=
            // A*(((position[0][2]+position[1][2]+position[2][2])/3)-parameter(2));
        }

        //  std:: cerr<<"total Volume:  "<<totalVolume<<std::endl;

        for (size_t n = 0; n < T.numCell(); ++n) {
            Cell cell = T.cell(n);
            if (cell.numVertex() != 3) {
                std::cerr
                    << "Pressure3D::Triangular::Triangular() "
                    << "Only implemented for triangular cells." << std::endl;
                exit(EXIT_FAILURE);
            }
            std::vector<double> normal = cell.getNormalTriangular(vertexData);
            double norm = 0.0;
            for (size_t d = 0; d < dimension; ++d) norm += normal[d] * normal[d];
            if (norm != 1.0) {
                norm = std::sqrt(norm);
                assert(norm > 0.0);
                double normFac = 1.0 / norm;
                for (size_t d = 0; d < dimension; ++d) normal[d] *= normFac;
            }

            // Get the cell size
            double A = 1.0;
            if (parameter(1) == 1.0)
                A = cell.calculateVolume(vertexData) / cell.numVertex();

            double coeff = parameter(0) * A;
            // update the vertex derivatives
            for (size_t k = 0; k < cell.numVertex(); ++k) {
                double vCoeff = coeff;
                // if (cell.vertex(k)->isBoundary(T.background()))
                // vCoeff *= 1.5;
                for (size_t d = 0; d < 3; ++d) {
                    vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
                }
            }
            // For saving normals in direction used for test plotting
            // for (size_t d=0; d<dimension; ++d)
            // cellData[cell.index()][d] = normal[d];
        }
    }
}

namespace CenterTriangulation {

Linear::
    Linear(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> > &indValue) {
    int n = paraValue.size();
    if ((n < 3) || (n > 4)) {
        std::cerr << "Pressure3D::CenterTriangulation::Linear::"
                  << "Linear() "
                  << "Uses three parameters: k_force and areaFlag and time span "
                  << "or four parameters: k_force, areaFlag, time span and starting pressure value"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[1] != 0.0 && paraValue[1] != 1.0 && paraValue[1] != 2.0) {
        std::cerr
            << "Pressure3D::CenterTriangulation::Linear::"
            << "Linear() "
            << "areaFlag must be zero (no area included) or one (area included) "
            << "or two(area included pressure in z direction only)." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (paraValue[2] < 0.0) {
        std::cerr << "Pressure3D::CenterTriangulation::Linear::"
                  << "Linear() "
                  << "time span must be positive." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (indValue.size() != 1) {
        std::cerr << "Pressure3D::CenterTriangulation::Linear::"
                  << "Linear() "
                  << std::endl
                  << "Start of additional Cell variable indices (center(x,y,z) "
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    setId("Pressure3D::CenterTriangulation::Linear");
    setParameter(paraValue);
    setVariableIndex(indValue);

    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_force";
    tmp[1] = "areaFlag";
    tmp[2] = "deltaT";

    timeFactor1 = 0.0;
    timeFactor2 = 0.0;
    setParameterId(tmp);
}

void Linear::derivs(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, DataMatrix &cellDerivs, DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs) {
    size_t dimension = vertexData[0].size();
    if (dimension != 3) {
        std::cerr << "Pressure3D::CenterTriangulation::Linear::"
                  << "derivs() "
                  << " Only implemented for three dimensions ." << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t numCells = T.numCell();
    size_t comIndex = variableIndex(0, 0);
    // size_t lengthInternalIndex = comIndex+dimension;

    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex)
        // if(cellData[cellIndex][27]!=-10)  // everywhere except l3 bottom
        // if(cellData[cellIndex][27]==-2) // only epidermis (L1 top)
        // if(cellData[cellIndex][28]==1 && cellData[cellIndex][29]==2) // only L2
        // top)  if(cellData[cellIndex][38]==-1) // for hypocotyl

        if (true) {
            size_t numWalls = T.cell(cellIndex).numWall();
            // Cell cell = T.cell(n);
            // double normal[3]={0,0,0};
            // normal[0]=cellData[n][normalVectorIndex  ];
            // normal[1]=cellData[n][normalVectorIndex+1];
            // normal[2]=cellData[n][normalVectorIndex+2];

            if (T.cell(cellIndex).numVertex() != numWalls) {
                std::cerr << "Pressure3D::CenterTriangulation::Linear::"
                          << "derivs() same number of vertices and walls."
                          << " Not for cells with " << T.cell(cellIndex).numWall()
                          << " walls and " << T.cell(cellIndex).numVertex()
                          << " vertices!" << std::endl;
                exit(EXIT_FAILURE);
            }
            // One triangle per 'vertex' in cyclic order
            for (size_t wallindex = 0; wallindex < numWalls; ++wallindex) {
                size_t kPlusOneMod = (wallindex + 1) % numWalls;
                // size_t v1 = com;
                size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
                size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();

                // size_t w1 = internal wallindex
                size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
                // size_t w3 = internal wallindex+1
                // Position matrix holds in rows positions for com, vertex(wallindex),
                // vertex(wallindex+1)
                DataMatrix position(3, vertexData[v2]);
                for (size_t d = 0; d < dimension; ++d)
                    position[0][d] = cellData[cellIndex][comIndex + d];  // com position
                // position[1] = vertexData[v2]; // given by initiation
                position[2] = vertexData[v3];
                // position[0][2] z for vertex 1 of the current element

                // Lengths are from com-vertex(wallindex),
                // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)),
                // com-vertex(wallindex+1)
                std::vector<double> length(numWalls);
                length[0] = std::sqrt((position[0][0] - position[1][0]) *
                                          (position[0][0] - position[1][0]) +
                                      (position[0][1] - position[1][1]) *
                                          (position[0][1] - position[1][1]) +
                                      (position[0][2] - position[1][2]) *
                                          (position[0][2] - position[1][2]));

                length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);

                length[2] = std::sqrt((position[0][0] - position[2][0]) *
                                          (position[0][0] - position[2][0]) +
                                      (position[0][1] - position[2][1]) *
                                          (position[0][1] - position[2][1]) +
                                      (position[0][2] - position[2][2]) *
                                          (position[0][2] - position[2][2]));

                // current Area of the element (using Heron's formula)
                double Area = std::sqrt((length[0] + length[1] + length[2]) *
                                        (-length[0] + length[1] + length[2]) *
                                        (length[0] - length[1] + length[2]) *
                                        (length[0] + length[1] - length[2])) *
                              0.25;

                double tempA = std::sqrt((position[2][0] - position[1][0]) *
                                             (position[2][0] - position[1][0]) +
                                         (position[2][1] - position[1][1]) *
                                             (position[2][1] - position[1][1]) +
                                         (position[2][2] - position[1][2]) *
                                             (position[2][2] - position[1][2]));

                double tempB = std::sqrt((position[0][0] - position[1][0]) *
                                             (position[0][0] - position[1][0]) +
                                         (position[0][1] - position[1][1]) *
                                             (position[0][1] - position[1][1]) +
                                         (position[0][2] - position[1][2]) *
                                             (position[0][2] - position[1][2]));

                double Xcurrent[3];
                Xcurrent[0] = (position[2][0] - position[1][0]) / tempA;
                Xcurrent[1] = (position[2][1] - position[1][1]) / tempA;
                Xcurrent[2] = (position[2][2] - position[1][2]) / tempA;

                double Bcurrent[3];
                Bcurrent[0] = (position[0][0] - position[1][0]) / tempB;
                Bcurrent[1] = (position[0][1] - position[1][1]) / tempB;
                Bcurrent[2] = (position[0][2] - position[1][2]) / tempB;

                double normal[3];
                normal[0] = Xcurrent[1] * Bcurrent[2] - Xcurrent[2] * Bcurrent[1];
                normal[1] = Xcurrent[2] * Bcurrent[0] - Xcurrent[0] * Bcurrent[2];
                normal[2] = Xcurrent[0] * Bcurrent[1] - Xcurrent[1] * Bcurrent[0];

                tempA = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                                  normal[2] * normal[2]);
                normal[0] = normal[0] / tempA;
                normal[1] = normal[1] / tempA;
                normal[2] = normal[2] / tempA;

                // normal[0]=cellData[cellIndex][32];
                // normal[1]=cellData[cellIndex][33];
                // normal[2]=cellData[cellIndex][34];

                // Get the cell size
                double A = 1.0 / 3;
                double coeff = 0.0;
                if (parameter(1) == 1.0 || parameter(1) == 2.0) A = Area / 2;

                // update the vertex derivatives
                if (parameter(1) == 0.0 || parameter(1) == 1.0) {
                    if (numParameter() == 3) {
                        coeff = timeFactor1 * parameter(0) * A;
                    } else if (numParameter() == 4) {
                        size_t pIndex = variableIndex(0, 1);
                        coeff = (parameter(3) + (timeFactor1 * (parameter(0) - parameter(3)))) * A;
                        cellData[cellIndex][pIndex] = parameter(3) + (timeFactor1 * (parameter(0) - parameter(3)));
                    }
                    cellDerivs[cellIndex][comIndex] += coeff * normal[0];
                    cellDerivs[cellIndex][comIndex + 1] += coeff * normal[1];
                    cellDerivs[cellIndex][comIndex + 2] += coeff * normal[2];

                    vertexDerivs[v2][0] += coeff * normal[0];
                    vertexDerivs[v2][1] += coeff * normal[1];
                    vertexDerivs[v2][2] += coeff * normal[2];

                    vertexDerivs[v3][0] += coeff * normal[0];
                    vertexDerivs[v3][1] += coeff * normal[1];
                    vertexDerivs[v3][2] += coeff * normal[2];
                }

                if (parameter(1) == 2.0) {
                    double coeff = timeFactor2 * parameter(0) * A;
                    cellDerivs[cellIndex][comIndex + 2] += coeff * normal[2];
                    vertexDerivs[v2][2] += coeff * normal[2];
                    vertexDerivs[v3][2] += coeff * normal[2];
                }

            }  // over walls
        }      // over cells
}

void Linear::update(
    Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
    DataMatrix &vertexData, double h) {
    static double deltat = 0;
    deltat += h;

    if (parameter(1) == 1 || parameter(1) == 0) {
        if (timeFactor1 < 1.0) timeFactor1 += h / parameter(2);

        if (timeFactor1 > 1.0) timeFactor1 = 1.0;
    }

    // if(parameter(1)==2 && deltat>800)
    if (parameter(1) == 2) {
        if (timeFactor2 < 1.0) timeFactor2 += h / parameter(2);

        if (timeFactor2 > 1.0) timeFactor2 = 1.0;
    }
}
}  // end namespace CenterTriangulation
}  // end namespace Pressure3D
