//
// Filename     : adhocReaction.h
// Description  : Classes describing some ad hoc updates
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : September 2007
// Revision     : $Id:$
//
#ifndef ADHOCREACTION_H
#define ADHOCREACTION_H

#include <cmath>

#include "baseReaction.h"
#include "myRandom.h"
#include "tissue.h"

/// @brief Simulate stresses in an inflation deflation experiment
///
/// @details description reaction is designed for the patchtest style simulations. it
/// applies a constant stress (s0 for t0) in x/y then a linear increse (to s1 over t1-t0)
/// for a time
/// then another constant stress (at s1 for t2 - t1) and then a linear decrease
/// (to s2 over t3 - t2)
/// s_y = (1 - s_aniso) * s_x
///
/// assumes vertices are labelled as:
/// @verbatim
///            sy
///            ^
///            |
///        tl-----tr
///        |       |
/// sx <-- |       | --> sx
///        |       |
///        bl-----br
///            |
///            v
///            sy
/// @endverbatim
/// In a model file the reaction is added as
/// @verbatim
/// InflationDeflationStresses 8 1 4
/// s0 s1 s2
/// t0 t1 t2 t3
/// s_aniso
/// bl br tl tr
/// @endverbatim
class InflationDeflationStresses : public BaseReaction {
   public:
    InflationDeflationStresses(
        std::vector<double> &paravalue,
        std::vector<std::vector<size_t> > &indvalue);

    void derivs(Tissue &T, DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                DataMatrix &cellDerivs,
                DataMatrix &wallDerivs,
                DataMatrix &vertexDerivs);

    void update(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &walldata,
                DataMatrix &vertexData,
                double h);

    double _time = 0.0;
};

/// @brief Sets positional derivatives to zero for vertices in specified region
///
/// @details A threshold is specified in a specific dimension and vertices on
/// either side of this is not updated.
///
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromPosition 2 1 1
/// threshold directionFlag
/// vertexPositionIndex
/// @endverbatim
/// where threshold sets the threshold in the dimension (x,y,z) specified with
/// vertexPositionIndex and directionFlag is set to 1 if vertices above the
/// threshold are to be kept fixed and -1 for vertices below the threshold.
///
/// @note This function sets the derivatives to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
class VertexNoUpdateFromPosition : public BaseReaction {
   public:
    VertexNoUpdateFromPosition(std::vector<double> &paraValue,
                               std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives to zero for vertices with listed indices
///
/// @details A list of vertex indices are specified for which vertex positions
/// are not. In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromIndex 0 1 N
/// vertexIndex1 [vertexIndex2...vertexIndexN]
/// @endverbatim
/// where the list of indices are the vertices not to be updated.
///
/// @note This function sets the derivatives to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
class VertexNoUpdateFromIndex : public BaseReaction {
   public:
    VertexNoUpdateFromIndex(std::vector<double> &paraValue,
                            std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives in the X direction to zero for vertices with listed indices
///
/// @details A list of vertex indices are specified for which vertex positions
/// are not In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromIndexHoldX 0 1 N
/// vertexIndex1 [vertexIndex2...vertexIndexN]
/// @endverbatim
/// where the list if indices are the vertices not to be updated in X direction.
///
/// @note This function sets the derivatives in X to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
class VertexNoUpdateFromIndexHoldX : public BaseReaction {
   public:
    VertexNoUpdateFromIndexHoldX(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives in the Y direction to zero for vertices with listed indices
///
/// @details A list of vertex indices are specified for which vertex positions
/// are not In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromIndexHoldY 0 1 N
/// vertexIndex1 [vertexIndex2...vertexIndexN]
/// @endverbatim
/// where the list if indices are the vertices not to be updated in Y direction.
///
/// @note This function sets the derivatives in Y to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
///
class VertexNoUpdateFromIndexHoldY : public BaseReaction {
   public:
    VertexNoUpdateFromIndexHoldY(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives in the Z direction to zero for vertices with listed indices
///
/// @details A list of vertex indices are specified for which vertex positions
/// are not In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromIndexHoldZ 0 1 N
/// vertexIndex1 [vertexIndex2...vertexIndexN]
/// @endverbatim
/// where the list if indices are the vertices not to be updated in Z direction.
///
/// @note This function sets the derivatives in Z to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
///
class VertexNoUpdateFromIndexHoldZ : public BaseReaction {
   public:
    VertexNoUpdateFromIndexHoldZ(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives to zero for vertices with listed indices
///
/// @details Freezes all of the vertices but the leading ones(at the tip).
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromList 0 0
/// @endverbatim
/// where the list if indices are the vertices not to be updated.
///
/// @note This function sets the derivatives to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
///
class VertexNoUpdateFromList : public BaseReaction {
   public:
    VertexNoUpdateFromList(std::vector<double> &paraValue,
                           std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);

    std::vector<size_t> updateVertices;
};

/// @brief randomizes the growth direction of tip cells(for two dim only)
///
/// @details Freezes all of the vertices but the leading ones(at the tip).
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexRandTip 1 1 2
/// max_randomization_angle
/// tip signal index
/// mt vector starting index(2d)
/// @endverbatim
/// where the list if indices are the vertices not to be updated.
///
/// @note This function sets the derivatives to zero, which means it has to be
/// provided after reactions that update the vertex derivatives.
///
class VertexRandTip : public BaseReaction {
   public:
    VertexRandTip(std::vector<double> &paraValue,
                  std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Sets positional derivatives to zero for vertices at boundary so that
/// boundary vertices would be restricted from moving in x and/or y ...
/// direction(s)
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundary 0 0
///
/// or for holding them only in x and z directions
/// VertexNoUpdateBoundary 0 1 2
/// 0 2
///
/// @endverbatim
///
class VertexNoUpdateBoundary : public BaseReaction {
   public:
    VertexNoUpdateBoundary(std::vector<double> &paraValue,
                           std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives to zero for vertices at boundary so that
/// boundary vertices would be restricted from moving parallel to the template
/// edges. The boundary condition gets updated in each step.(not completely
/// debugged)
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundaryPtemplate 0 1 1
/// comIndex-cellVector
///
/// @endverbatim
///
class VertexNoUpdateBoundaryPtemplate : public BaseReaction {  // BB

   public:
    VertexNoUpdateBoundaryPtemplate(std::vector<double> &paraValue,
                                    std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Sets positional derivatives to zero for vertices at boundary so that
/// boundary vertices would be restricted from moving parallel to the initial
/// template edges the boundary is static.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundaryPtemplateStatic 0 1 1
/// comIndex-cellVector
///
/// @endverbatim
///
class VertexNoUpdateBoundaryPtemplateStatic : public BaseReaction {  // BB

   public:
    VertexNoUpdateBoundaryPtemplateStatic(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> > &indValue);
    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

   private:
    std::vector<size_t> boundaryVertices;
    std::vector<std::vector<double> > boundaryNormal;
    size_t numBoundaryVertices;
};

/// @brief Sets positional derivatives to zero for vertices at boundary so that
/// boundary vertices would be restricted from moving parallel to the initial
/// template edges the boundary is static.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundaryPtemplateStatic3D 0 1 1
/// comIndex-cellVector
///
/// @endverbatim
///
class VertexNoUpdateBoundaryPtemplateStatic3D : public BaseReaction {  // BB

   public:
    VertexNoUpdateBoundaryPtemplateStatic3D(
        std::vector<double> &paraValue,
        std::vector<std::vector<size_t> > &indValue);
    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

   private:
    //  std::vector<size_t> bottomVertices;
    //  std::vector<size_t> sideVertices;

    std::vector<std::vector<double> > bottomNormals;
    std::vector<std::vector<double> > sideNormals;
    size_t numBottomCells;
    size_t numSideCells;
};

/// @brief Sets positional derivatives to zero for vertices at boundary so that
/// boundary vertices would be restricted from moving parallel to the initial
/// template edges the boundary is static.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundary3D 0 1 2
/// comIndex-cellVector
/// neighborIndex
///
/// @endverbatim
///
class VertexNoUpdateBoundary3D : public BaseReaction {  // BB

   public:
    VertexNoUpdateBoundary3D(std::vector<double> &paraValue,
                             std::vector<std::vector<size_t> > &indValue);
    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

   private:
    std::vector<double> bottomVertices, sideVertices;
};

/// @brief Sets positional derivatives to zero for vertices at boundary so that
/// boundary vertices would be restricted from moving parallel to the initial
/// template edges the boundary is static.
///
/// @details In the model file, the reaction is specified as
/// (for holding the boundary vertices in all directions):
/// @verbatim
/// VertexFromConstStressBoundary 8 0
/// stress x
/// stress y
/// dx
/// dy
/// boundary right  x
/// boundary left   x
/// boundary top    y
/// boundary bottom y
/// @endverbatim
///
class VertexFromConstStressBoundary : public BaseReaction {  // BB

   public:
    VertexFromConstStressBoundary(std::vector<double> &paraValue,
                                  std::vector<std::vector<size_t> > &indValue);
    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);

   private:
    std::vector<std::vector<double> > rightVertices, leftVertices, topVertices,
        bottomVertices;
    size_t numOldVertices;
    double totaltime;
};

/// @brief For ad-hoc manipulation of the templates
///
/// @details In the model file, the reaction is specified as
/// (for holding the boundary vertices in all directions):
/// @verbatim
/// manipulate 0 0
/// @endverbatim
/// @note Do not use unless you check the adHoc code.
class manipulate : public BaseReaction {  // BB

   public:
    manipulate(std::vector<double> &paraValue,
               std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Calculates cell polarity vector based on a vector and a measure for
/// each cell wall and places it in the cell center. The vector field could be
/// printed via a seperate vtk file using its own printState function.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// cellPolarity3D 0 2 2 2
///
/// CELL_IDENTITY_INDEX
/// comIndex
///
/// polarity_vector_index
/// polarity_measure_INDEX
///
/// @endverbatim
///
class cellPolarity3D : public BaseReaction {  // BB

   public:
    cellPolarity3D(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);
    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
    void printState(Tissue *T, DataMatrix &cellData, DataMatrix &wallData,
                    DataMatrix &vertexData, std::ostream &os);

   private:
    std::vector<std::vector<double> > cellFaces;
    std::vector<std::vector<double> > cellCentPol;
};

/// @brief Diffusion in 3 dimensional templates assuming the same volume of the
/// cells (for now ad-hoc)
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// diffusion3D 1 1 3
/// diff-constant
///
/// conc_index
/// neighborWall_index
/// 3Dcell_index
/// @endverbatim
///
class diffusion3D : public BaseReaction {  // BB

   public:
    diffusion3D(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> > &indValue);
    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

   private:
    std::vector<std::vector<double> >
        Cells3d;                                 // holds the wall_indices and neighbohrhood info
    std::vector<std::vector<double> > faceArea;  // holds the area of the faces
};

/// @brief Moves the complete tissue such that the maximal value in specified
/// direction (x/y/z) is constant (defined by maxPos).
///
/// @details The translation is done at each update (i.e. after each ODE integration
/// step). In a model file the reaction is defined as: 
/// @verbatim
/// VertexTranslateToMax 1 1 1
/// p_0 maxPos
/// vI_00 positional index
/// @endverbatim
class VertexTranslateToMax : public BaseReaction {
   public:
    VertexTranslateToMax(std::vector<double> &paraValue,
                         std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Centers the tissue such that the center of mass is in origo
///
/// @details The translation is done at each update (i.e. after each ODE
/// integration step), mainly for plotting by keeping the tissue centered.
///
/// In a model file the reaction is given by:
/// @verbatim
/// CenterCOM 0 0
/// @endverbatim
///
class CenterCOM : public BaseReaction {
   public:
    CenterCOM(std::vector<double> &paraValue,
              std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Centers the tissue such that the center of mass of specific cell is in origo
///
/// @details The translation is done at each update (i.e. after each ODE
/// integration step), mainly for plotting by keeping the tissue centered around a specific cell.
/// Handy for following e.g. an apical cell.
///
/// In a model file the reaction is given by:
/// @verbatim
/// CenterCellCOM 0 0
/// @endverbatim
/// where it assumed to center around cell 0, or
/// @verbatim
/// CenterCellCOM 1 0
/// cellIndex
/// @endverbatim
/// where the centering will be around cell given as parameter.
///
class CenterCellCOM : public BaseReaction {
   public:
    CenterCellCOM(std::vector<double> &paraValue,
              std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Centers the tissue such that the center of mass is in origo including
/// for central mesh points.
///
/// @details The translation is done at each update (i.e. after each ODE
/// integration step). This function should be used in connection with
/// mechanical TRBScenterTriangulation since those reactions adds an extra
/// vertex in the cell data vector, which also have to be moved. In a model file
/// the reaction is defined as
/// @verbatim
/// CenterCOMcenterTriangulation 0 1 1
/// cellDataPositionIndex
/// @endverbatim
/// where the index have to be matched with what is given in the mechanical TRBS
/// reaction.
///
class CenterCOMcenterTriangulation : public BaseReaction {
   public:
    CenterCOMcenterTriangulation(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Helper reaction that only calculates the PCA plane for every cell.
///
/// Since calculating the PCA plane is computationally expensive, this
/// reaction is calculating it once for all, such that other reaction in need
/// of the PCA plane does not have to do the calculation.
///
/// Note that this reaction has to be before the reactions using the PCA
/// plane.
///
class CalculatePCAPlane : public BaseReaction {
   public:
    CalculatePCAPlane(std::vector<double> &paraValue,
                      std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Initiate the wall length variables to a factor times the distance
/// between the two vertices
///
/// @details Sets the wall length variables to be the distance between the two
/// vertices times a factor. In the model file the reaction is defined as:
/// @verbatim
/// InitiateWallLength 1 0
/// factor
/// @endverbatim
///
class InitiateWallLength : public BaseReaction {
   public:
    InitiateWallLength(std::vector<double> &paraValue,
                       std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Initiate the target area variable to a factor times the current cell
/// volumes.
///
/// @details Sets the target area variables to a factor times the current
/// corresponding cell volume. Note that this method is independent of the
/// number of dimensions, and can thus be equivalently used to set a target
/// volume. In the model file the reaction is defined as:
/// @verbatim
/// InitiateTargetArea 1 0
/// factor
/// @endverbatim
///
class InitiateTargetArea : public BaseReaction {
   public:
    InitiateTargetArea(std::vector<double> &paraValue,
                       std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Initiate a wall variable to a constant value for all walls
///
/// @details Initiate a wall variable to a constant value for all walls
/// In the model file the reaction is defined as:
/// @verbatim
/// InitiateWallVariableConstant 1 1 1
/// factor
/// c_wallindex
/// @endverbatim
///
class InitiateWallVariableConstant : public BaseReaction {
   public:
    InitiateWallVariableConstant(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Adds additional vertices to all walls
///
/// @details Wall 'meshing' is applied by inserting additional vertices for all
/// walls. This creates new vertices connected to one of the previous walls that
/// are now divided into several walls. Only a single parameter is given which
/// sets the number of new vertices per wall. The 'remeshing' is done in the
/// initiation step and there is no contribution to the derivative for this
/// reaction.
///
/// In the model file the reaction is defined as:
/// @verbatim
/// InitiateWallMesh 1 0
/// numVertex
/// @endverbatim
///
class InitiateWallMesh : public BaseReaction {
   public:
    InitiateWallMesh(std::vector<double> &paraValue,
                     std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

/// @brief Different well-defined derivatives used to check the strain direction
/// update
///
class StrainTest : public BaseReaction {
   public:
    StrainTest(std::vector<double> &paraValue,
               std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

class CalculateVertexStressDirection : public BaseReaction {
   public:
    CalculateVertexStressDirection(std::vector<double> &paraValue,
                                   std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);

   private:
    std::vector<size_t> wallForceIndexes_;
};

class MoveVerticesRandomlyCapCylinder : public BaseReaction {
   public:
    MoveVerticesRandomlyCapCylinder(std::vector<double> &paraValue,
                                    std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief scales the template by a factor via Initiate
/// position of vertices and wall length variables will be scaled.
///
/// @details In the model file the reaction is defined as:
/// @verbatim
/// scaleTemplate 1 0
/// factor
/// @endverbatim
///
class scaleTemplate : public BaseReaction {
   public:
    scaleTemplate(std::vector<double> &paraValue,
                  std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Copies vectors from one index to another in the cell vector
/// at initiation (4 component after the indices will be copied).
///
/// @details This reaction copies a vector (4 components) from one place to
/// another at the initiation of the simulation. The start index of the vector
/// to be copied and where it will be copied to needs to be provided, and then
/// these four values are copied for all cells.
///
/// In the model file the reaction is defined as:
/// @verbatim
/// copyCellVector 0 1 2
/// copy_from_index copy_to_index
/// @endverbatim
///
class copyCellVector : public BaseReaction {
   public:
    copyCellVector(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief scales the template by a factor via Initiate
/// randomizes the MT direction of the cells within the cell plane
///
/// @details In the model file the reaction is defined as:
/// @verbatim
/// randomizeMT 0 1 1
/// MT_index
/// @endverbatim
///
class randomizeMT : public BaseReaction {
   public:
    randomizeMT(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

class restrictVertexRadially : public BaseReaction {
   public:
    restrictVertexRadially(std::vector<double> &paraValue,
                           std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief A molecule is produced/created with constant rate in the cells
/// located around CZ to creat phyllotactic pattern.
///
/// @f[ \frac{dc}{dt} = k_c/Area @f] if cell index is in the given list
/// Area will be  one if constant concentration(instead of amount) is needed to
/// be produced
///
/// In a model file the reaction is defined as
/// @verbatim
/// creationPrimordiaTime 5 1 2
/// k_c # creation rate
/// tp # time delay
/// teta # delta teta
/// z # vertical distance from max
/// constant-amount-flag
///
/// c_index
/// com_index
/// @endverbatim
///
class CreationPrimordiaTime : public BaseReaction {
   private:
    std::vector<size_t> proCells;

   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    CreationPrimordiaTime(std::vector<double> &paraValue,
                          std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Updates list of vertices with a given force applied where the force
/// is linearly increased from zero across a given time span (deltaT).
///
/// @details In a model file the reaction is defined as
/// @verbatim
/// VertexFromRotationalForceLinear 1/2/3(dimension+1) 1 (no of vertices)
/// Force component(s) deltaT
/// 1st vertex index
/// 2nd vertex index
/// ...
/// @endverbatim
///
class VertexFromRotationalForceLinear : public BaseReaction {
   private:
    double timeFactor_;
    std::vector<std::vector<double> > boundVerticesUp, boundVerticesDn;

   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    VertexFromRotationalForceLinear(std::vector<double> &paraValue,
                                    std::vector<std::vector<size_t> > &indValue);

    void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief Switches an output variable from 0 to 1 if an input variable is
/// above a threshold.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdSwitch 2 2 1 1   # number of parameters is set to two (threshold
/// and switch_type)
/// threshold		      # threshold above which a variable is reset to zero.
/// switch_type		      # the switch_type parameter takes the values 0 and 1
///                           # for defining the reversible and irreversible
///                           # switch, respectively.
/// index_var   	      # index of the index variable upstream the switch.
/// index_var_out  	      # list of updated indices - for the moment it can
///                           # be just one index - where the output of
///                           # the switch is written. @endverbatim
/// @note This function makes a downstream species reversibly or irreversibly
/// switch from 0 to 1, upon being above a certain threshold of an upstream
/// variable.
///
class ThresholdSwitch : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    ThresholdSwitch(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species reversibly or irreversibly
/// switch from 1 to 0, upon being above a certain threshold of an upstream variable.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdReset 2 2 1 1 # number of parameters is set to two
/// threshold		   # threshold above which a variable is reset to zero.
/// switch_type		   # the switchtype parameter takes the values 0 and 1
///                        # for defining the reversible and irreversible switch, respectively.
/// index_var   	   # index of the index variable upstream the switch.
/// index_var_out  	   # index of updated variable @endverbatim
///
class ThresholdReset : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    ThresholdReset(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species irreversibly switch from 1 to 0 (reset)
/// with a bit of noise, upon being above a certain threshold of an upstream variable.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdNoisyReset 2 2 1 1 # number of parameters is set to two
/// threshold		 	# threshold above which a variable is reset to zero.
/// noise_amplitude             # noise amplitude for the resetting
/// index_var   	 	# index of the index variable upstream the switch.
/// index_var_out  	        # index of updated variable @endverbatim
///
class ThresholdNoisyReset : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    ThresholdNoisyReset(std::vector<double> &paraValue,
                        std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species reversibly or irreversibly switch from 1
/// to 0, upon being above a certain threshold of an upstream variable.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdResetAndCount 2 2 1 1 # number of parameters is set to two
/// threshold		           # threshold above which a variable is reset to zero.
/// switch_type                    # the switchtype parameter takes the values 0 and 1 for
///                                # defining the reversible and irreversible switch, respectively.
/// index_var                      # index of the index variable upstream the switch.
/// index_var_out                  # index of updated variable @endverbatim
///
class ThresholdResetAndCount : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    ThresholdResetAndCount(std::vector<double> &paraValue,
                           std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species reversibly or irreversibly
/// switch from 1 to 0, upon being above a certain threshold of an upstream variable.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// FlagNoisyReset 2 2 1 1 # number of parameters is set to two
/// flag_value		   # value of the flag that will make the output variable resetting
///                        # to zero with noise.
/// switch_type            # the switchtype parameter takes the values 0 and 1 for defining the
///                        # reversible and irreversible switch, respectively.
/// index_var              # index of the index variable upstream the switch.
/// index_var_out          # index of updated variable @endverbatim
///
class FlagNoisyReset : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    FlagNoisyReset(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species be reset to 0 with noise, when being above
/// a certain threshold of an upstream variable and another flag variable has a certain flag_value
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdAndFlagNoisyReset 3 2 1 1 # number of parameters is set to three
/// threshold		 	       # threshold above which a variable is reset to zero.
/// flag_value                         # value of the flag that will make the output variable
///                                    # is resetting to zero with noise.
/// noise_amplitude		       # noise amplitude for the resetting
/// index_var_in1   	 	       # index threshold variable upstream the switch.
/// index_var_in2   	 	       # index flag variable upstream the switch.
/// index_var_out  		       # index of updated variable @endverbatim
///
class ThresholdAndFlagNoisyReset : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    ThresholdAndFlagNoisyReset(std::vector<double> &paraValue,
                               std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species add an add_value if a
/// certain flag variable is 1.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// FlagAddValue 1 2 1 1 # number of parameters is set to one
/// add_value  		 # number that will be added
/// index_var_in  	 # index of the flag variable.
/// index_var_out  	 # index where the output of the gate is written. @endverbatim
///
class FlagAddValue : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    FlagAddValue(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function  copies one input variable into the output variable
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// CopyVariable 0 2 1 1 # number of parameters is set to 0
/// index_var   	 # input index variable
/// index_var_out  	 # output index variable
/// @endverbatim
///
class CopyVariable : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    CopyVariable(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

//
// @brief This reaction should not be used (restricted use) unless you are a developer.
//
class DebugReaction : public BaseReaction {
   public:
    DebugReaction(std::vector<double> &paraValue,
                  std::vector<std::vector<size_t> > &indValue);

    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

#endif  // ADHOCREACTION_H
