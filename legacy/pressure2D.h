//
// Filename     : pressure2D.h
// Description  : Classes describing forces acting on 2D surfaces
//              : (perpendicular to edges)
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : November 2019
// Revision     : $Id:$
//
#ifndef PRESSURE2D_H
#define PRESSURE2D_H

#include <cmath>

#include "baseReaction.h"
#include "tissue.h"

///
/// @brief Namespace for reactions updating vertices based on 2D cell pressure
/// forces, i.e. perpendicular to edges, always updating in 2D.
///
/// These functions are for generating forces from internal cell pressures, and
/// includes reactions describing perpendicular to edge forces and
/// forces from a potential based on area increase.
///
/// @note Might be interesting for some of the old AreaPotential versions to use the
/// more straightforward area calculation now implemented.
/// @note These reactions assumes a 2D tissue, sometimes allow 3D tissue but do the
/// update in two dimensions.
///
namespace Pressure2D {
///
/// @brief Updates vertices from a cell pressure potential implemented as an
/// area-based potential in 2D
///
/// @details This reaction is a 2D version of a pressure force calculated from a
/// potential given as an area expansion
/// @f[ U(v_{i}) = - \frac{1}{2} p_{0} A(v_{i}) @f]
/// where @f$v_{i}@f$ are the vertex positions for the cell, A is the area and
/// @f$p_{0}@f$ is a constant pressure. The time derivative of vertex positions
/// are then calculated as positional derivative of the potential
/// @f[ \frac{dv_{ix}}{dt} = - \frac{dU}{dv_{ix}}@f]
/// @f[ \frac{dv_{iy}}{dt} = - \frac{dU}{dv_{iy}}@f]
/// The area is calculated by the formula
/// @f[ 2A = abs(\sum_{i} v_{ix}v_{(i+1)y} - v_{(i+1)x}v_{iy}) @f]
/// where the vertices is covered in a circular fashion (@f$v_{N}=v_{0}@f$).
/// The forces on vertex @f$v@f$ from a cell is then given by
///
/// @f[\frac{dx_v}{dt} = 0.5*p_0 (y_{v_r} - y_{v_l}) @f]
/// @f[\frac{dy_v}{dt} = 0.5*p_0 (x_{v_l} - x_{v_r}) @f]
///
/// where @f$v_r,v_l@f$ are right and left vertices in the sorted order.
/// @f$p_0@f$ represents the pressure, and if @f$p_{1}=1@f$, the pressure will
/// be divided by the cell volume. The Force is in the 'outward' normal
/// direction for the two edges connected to the vertices and proportional to
/// the length of the edge.
///
/// In a model file, the reaction is given by:
/// @verbatim
/// Pressure2D::AreaPotential 2 0
/// P V_normflag(=0/1)
/// @endverbatim
/// @note Requires two dimensions with vertices sorted.
///
class AreaPotential : public BaseReaction {
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
    AreaPotential(std::vector<double> &paraValue,
                  std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Updates vertices from a cell pressure potential described as an area
/// expansion in 2D
///
/// @details This reaction is a 2D version of a pressure force calculated from a
/// potential given as an area expansion
/// @f[ U(v_{i}) = - \frac{1}{2} p_{0} A(v_{i}) @f]
/// where @f$v_{i}@f$ are the vertex positions for the cell, A is the area and
/// @f$p_{0}@f$ is a constant pressure. The time derivative of vertex positions
/// are then calculated as positional derivative of the potential
/// @f[ \frac{dv_{ix}}{dt} = - \frac{dU}{dv_{ix}}@f]
/// @f[ \frac{dv_{iy}}{dt} = - \frac{dU}{dv_{iy}}@f]
/// The area (expansion) is calculated for each triangle described
/// by an edge (two vertex positions, @f$v_{1},v_{2}@f$) and the center of mass
/// @f$x_{c}@f$, of the cell. The Area is given as @f$\frac{1}{2} h n@f$, where
/// n is the length of the edge (base of triangle) and h is the height. The
/// height 'vector' is extracted by first finding The vector from the center of
/// the edge, @f$x_{0}@f$, and the center of mass
/// @f[ dx = x_{c} - x_{0} @f]
/// then project this down onto the edge vector @f$n@f$ followd by extracting
/// the height vector @f$h@f$ perpendicular to the edge for the area calculation
/// @f[ h = dx + dx \frac{n}{|n|} @f]
/// In a model file, the reaction is given by:
/// @verbatim
/// Pressure2D::AreaPotentialTri 2[/3] 0
/// P V_normFlag(=0/1) [InternalCellsOnlyFlag(=0/1)]
/// @endverbatim
/// where @f$p_{0}=P@f$ is the pressure magnitude, @f$p_{1}@f$ is a flag set to
/// 1 if the calculated forces should be normalised with the cell area, and
/// @f$p_{2}@f$ is an optional flag set to 1 if the pressure should be applied
/// only to internal cells, i.e. non-epidermal cells not connected to the
/// background.
///
/// @note The secondary effect to the area of the movement of the center of mass
/// when a vertex is moved is not taken into account.
/// @note Maybe convert to the area description used in Pressure2D::AreaPotential?
/// @note Requires two dimensions.
///
class AreaPotentialTri : public BaseReaction {
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
    AreaPotentialTri(std::vector<double> &paraValue,
                     std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Updates vertices from a cell pressure potential described as an area
/// expansion if the cell is close enough to the maximal position
///
/// @details This reaction uses for each individual cell the same update as
/// Pressure2D::AreaPotentialTri and the only difference is an extra parameter
/// setting a threshold in space where the update is only done if the cell is
/// closer to the maximal position than this given threshold variable in the
/// direction given as (only) variable index. The idea is for example to define
/// a growth zone close to the apex of a tissue that follows the tip as the
/// tissue is growing.
/// In a model file, the reaction is given by:
/// @verbatim
/// Pressure2D::AreaPotentialTriSpatialThreshold 2 1 1
/// P Threshold
/// direction_index
/// @endverbatim
/// where P=@f$p_{0}@f$ is the constant pressure, Threshold=@f$p_{1}@f$ is the
/// spatial threshold for how close to the maximal position a cell needs to be
/// to be updated, and direction_index is the spatial direction the max and
/// threshold are calculated in.
///
/// @see Pressure2D::AreaPotentialTri
/// @note Currently, this reaction has not implemented a flag for area
/// normalized forces.
/// @note Maybe convert to the area description used in Pressure2D::AreaPotential?
/// @note Requires two dimensions.
///
class AreaPotentialTriSpatialThreshold : public BaseReaction {
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
    AreaPotentialTriSpatialThreshold(std::vector<double> &paraValue,
                                     std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Updates vertices with forces from a cell pressure potential described
/// by an area expansion and with an additional aim of a target area in 2D
///
/// @details This reaction is a 2D version of a pressure force calculated from a
/// potential given as an area expansion
/// @f[ U(v_{i}) = - \frac{1}{2} p_{0} A(v_{i}) @f]
/// where @f$v_{i}@f$ are the vertex positions for the cell, A is the area and
/// @f$p_{0}@f$ is a constant pressure. The time derivative of vertex positions
/// are then calculated as positional derivative of the potential
/// @f[ \frac{dv_{ix}}{dt} = - \frac{dU}{dv_{ix}}@f]
/// @f[ \frac{dv_{iy}}{dt} = - \frac{dU}{dv_{iy}}@f]
/// The area is calculated by the formula
/// @f[ 2A = abs(\sum_{i} v_{ix}v_{(i+1)y} - v_{(i+1)x}v_{iy}) @f]
/// where the vertices is covered in a circular fashion (@f$v_{N}=v_{0}@f$).
/// In addition, a target area, @f$A_{w}@f$ is given in a cell variable
/// and the update is multiplied by the factor
/// @f[ (1-\frac{A}{A_{w}}) @f]
/// and a flag for if area decrease is allowed is given as a second parameter,
/// i.e. if p_{1}=1, there will be no update leading to decreasing cell area.
/// The target area can represent an approximation of a water volume the cells
/// try to adapt to, which can be estimated from the relation between optimal
/// pressure and current estimate of the pressure as in the reaction
/// TargetAreaFromPressure. In a model file, the reaction is given by:
/// @verbatim
/// Pressure2D::AreaPotentialTargetArea 2 1 1
/// P flag_AreaDecrease(=0/1)]
/// AreaTarget_index
/// @endverbatim
/// where @f$p_{0}=P@f$ is the pressure magnitude, @f$p_{1}@f$ is a flag set to
/// 1 if a decrease in area is allowed.
///
/// @see Pressure2D::AreaPotential
/// @see TargetAreaFromPressure
/// @note Requires two dimensions.
///
class AreaPotentialTargetArea : public BaseReaction {
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
    AreaPotentialTargetArea(std::vector<double> &paraValue,
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
    /// @brief Calculates the area [should this be replaced with the
    /// Cell.calculateVolume(vertexData) ]?
    ///
    double polygonArea(std::vector<std::pair<double, double> > vertices);
};

}  // end namespace Pressure2D

#endif  // PRESSURE2D_H
