//
// Filename     : mechanical.h
// Description  : Classes describing mechanical updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef MECHANICAL_H
#define MECHANICAL_H

#include <cmath>

#include "baseReaction.h"
#include "tissue.h"

namespace CenterTriangulation {
///
/// @brief Updates vertices from a 2D cell pressure potential, i.e. perpendicular to edges
///
/// @details This function determines the direction of the pressure force term
/// from the position of the central mesh cell vertex to the center of the wall.
/// Applies a force proportional to a constant pressure (parameter(0)) and the
/// size of the wall. parameter(1) equal to 1 normalizes the force with cell
/// volume. (0 otherwise).
///
/// In a model file, the reaction is given by:
/// @verbatim
/// CenterTriangulation:VertexFromCellPressure 2 1 2
/// P V_normflag(=0/1)
/// startIndex
/// concentraion_index(if 0 pressure is constant)
/// @endverbatim
///
/// where the startindex is marking the start of internal edge varibales
/// (x,y,z,L_1,...).
///
/// @note Maybe it should rather be normal to the wall in the plane of the
/// triangle?
/// @note Assumes three dimensions as all CenterTriangulation functions.
/// @see VertexFromCellPressure
///
class VertexFromCellPressure : public BaseReaction {
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
    VertexFromCellPressure(std::vector<double> &paraValue,
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
/// @brief Updates vertices from a cell pressure potential linearly increasing
/// in a given time span
///
/// @details This function determines the direction of the pressure force term
/// from the position of the central mesh cell vertex to the center of the wall.
/// Applies a force proportional to the pressure (parameter(0)) and the
/// size of the wall.
///
/// In a model file the reaction is defined as:
/// @verbatim
/// CenterTriangulation:VertexFromCellPressureLinear 3 1 1
/// Pressure
/// V normalized_flag (0/1)
/// deltaT
///
/// COM index
/// @endverbatim
/// @note Maybe it should rather be normal to the wall in the plane of the
/// triangle?
///
class VertexFromCellPressureLinear : public BaseReaction {
   private:
    double timeFactor_;
    double totaltime;

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
    VertexFromCellPressureLinear(std::vector<double> &paraValue,
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
}  // end namespace CenterTriangulation

/// @brief This reaction estimates a 'water volume' change given a target
/// pressure from current pressure estimated from wall tension
///
/// @details This reaction updates a cell variable estimating a target area that
/// can be used in combination with updating vertex positions towards such a
/// target area, e.g. with reaction Pressure2D::AreaPotentialTargetArea. The
/// idea is to represent a constant target pressure given by
/// @f$P_{target}=p_{1}@f$, and estimate the current pressure, @f$P@f$ from wall
/// forces (as read from wall variables and calculated elsewhere):
/// @f[ P = k_{pp} \sum_w \frac{F_{w}}{L_{w}} @f]
/// where @f$p_{2}=k_{pp}@f$ is a scaling/normalisation factor for the forces,
/// @f$F_{w}@f$ the forces read from the wall variables and @f$L_{w}@f$ is the
/// wall length. The cell target area variable, @f$A_{c}@f$, is updated by
/// @f[ \frac{dA_{c}}{dt} = k_{p} (P_{target} - P) \sum_{w} L_{w} @f]
/// where @f$k_{p}=p_{0}@f$ is the rate of the update. In addition, a flag
/// (@f$p_{3}=1@f$) can be given to disallow the target area to shrink. In a
/// model file the reaction is defined as:
/// @verbatim
/// TargetAreaFromPressure 4 2 2 n
/// k_p  P_max  k_pp flag_allowShrink(=0/1)
/// WallLength_index  targetArea_index
/// Force indices
/// @endverbatim
/// or if The calculated pressure is saved
/// @verbatim
/// TargetAreaFromPressure 4 3 2 n 1
/// k_p  P_max  k_pp flag_allowShrink(=0/1)
/// WallLength_index  targetArea_index
/// Force indices
/// Cell_index_for_saving_the_pressure
/// @endverbatim
///
/// @see Pressure2D::AreaPotentialTargetArea
/// @note Find reference for the estimate of P and following taget area update.
///
class TargetAreaFromPressure : public BaseReaction {
   public:
    TargetAreaFromPressure(std::vector<double> &paraValue,
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
/// @brief Updates vertices from cells via a power diagram potential
///
/// @details Power diagrams is a version of voronoi-tesselation where different 'cells'
/// are allowed to have different sizes. See e.g.
/// F. AURENHAMMER (1987) POWER DIAGRAMS: PROPERTIES, ALGORITHMS AND APPLICATIONS
/// SIAM J Computing 16:78-96.
/// Given that cell sizes are stored, power diagrams can be used to calculate an 'optimal'
/// position of a vertex connected to threee cells. This function then updates the vertex
/// position towards this optimal position with a rate @f$p_0@f$. In the model file
/// the reaction is given by
/// @verbatim
/// VertexFromCellPowerDiagram 1 1 1
/// K_force
/// cellSizeIndex
/// @endverbatim
///
/// @note only works for 2D currently
///
class VertexFromCellPowerdiagram : public BaseReaction {
   public:
    VertexFromCellPowerdiagram(std::vector<double> &paraValue,
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
/// @brief A cell 'pressure' reaction providing forces perpendicular to walls
/// (with magnitude given by cell variable)
///
class PerpendicularWallPressure : public BaseReaction {
   public:
    PerpendicularWallPressure(std::vector<double> &paraValue,
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
/// @brief Collection of reactions or updating vertices via forces acting perpendicular
/// to faces.
///
/// @details Pressure forces acting perpendicular to faces (Cells in 2.5D). Forces
/// are proportional to the face area, and several versions of normalisation or
/// contribution from cell variables are provided. Importantly, except for the versions
/// using CenterTriangulation, or Triangular, the normal to the plane is approximated
/// via a PCA plane extraction and forces are perpendicular to this.
///
/// @see CalculatePCAPlane : required for those that applies forces on 3+vertex faces
/// @note The reactions in this namespace used to be called VertexFromCellPlane* (Pressure3D::*).
///
namespace Pressure3D {

///
/// @brief Updates vertices from an 'internal pressure' term defined to act in the cell
/// (face) normal direction.
///
/// @details This function calculates the area of a cell and then distribute a
/// force 'outwards' among the cell vertices. It relies on that the PCA cell
/// planes have been calculated. A cell contributes to a vertex update with
///
/// @f[ \frac{dx_{i}}{dt} = p_{0} A n_{i} / N_{vertex} @f]
///
/// where @f$p_{0}@f$ is a 'pressure' parameter, A is the cell area @f$n_{i}@f$
/// is the cell normal component and @f$N_{vertex}@f$ is the number of vertices
/// for the cell. An additional parameter @f$p_{2}@f$ can be used to not include
/// the area factor if set to zero (normally it is set to 1).
///
/// In a model file the reaction is defined as
/// @verbatim
/// Pressure3D::Constant 2 0
/// P A_flag
/// @endverbatim
///
/// @see Cell.getNormalToPCAPlane()
/// @note Used to be named VertexFromCellPlane
///
class Constant : public BaseReaction {
   public:
    Constant(std::vector<double> &paraValue,
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
/// @brief Updates vertices from a 'pressure' term defined to act in the cell
/// (face) normal direction. The pressure is applied increasingly (linear incease
/// in a given time span).
///
/// @details This function calculates the area of a cell and then distribute a
/// force 'outwards' among the cell vertices. It relies on that the PCA cell
/// planes have been calculated. A cell contributes to a vertex update with
/// @f[ \frac{dx_{i}}{dt} = p_{0} A n_{i} / N_{vertex} @f]
/// where @f$p_{0}@f$ is a 'pressure' parameter, A is the cell area @f$n_{i}@f$
/// is the cell normal component and @f$N_{vertex}@f$ is the number of vertices
/// for the cell. An additional parameter @f$p_{2}@f$ can be used to not include
/// the area factor if set to zero (normally it should be set to 1).
/// In a model file the reaction is defined as
/// @verbatim
/// Pressure3D::Linear 3 0
/// P A_flag deltaT
/// @endverbatim
///
/// @see Cell.getNormalToPCAPlane()
/// @see Pressure3D::Constant
/// @see Pressure3D
/// @note Used to be called VertexFromCellPlaneLinear
///
class Linear : public BaseReaction {
   private:
    double timeFactor_;

   public:
    Linear(std::vector<double> &paraValue,
           std::vector<std::vector<size_t> > &indValue);

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

///
/// @brief Updates vertices from a 'pressure' term defined to act in the cell
/// normal direction, similar to Pressure3D::Constant, but with the Force
/// depending on the distance from the maximal position in a given direction.
///
/// @details This function calculates the area of a cell and then distribute a
/// force 'outwards' among the cell vertices. It is similar to Pressure3D::Constant
/// but the magnitude of the force has a
/// spatial dependence, implemented as a (decreasing) Hill function as a function
/// of the distance to a vertex with a maximal position in a given direction (e.g.
/// the distance to the tip of a meristem, where the tip is defined by the vertex
/// located with maximal z-value). It relies on that the PCA cell planes have been
/// calculated. The force is described by:
/// \f[
/// F = \frac{A_{cell}}{N_{vertex}}(p_{0} +
/// p_{1} \frac{p_{2}^{p_{3}}}{p_{2}^{p_{3}}+D^{p_{3}}}
/// \f]
/// where the area factor \f$A_{cell}\f$ is present if \f$p_{4}=1\f$,
/// and D is the the Euclidean distance to the vertex maximal in a
/// specified dimension.
/// In a model file the reaction is defined as
/// @verbatim
/// Pressure3D::Spatial 5 1 1
/// P_const P_add K_H n_H A_flag
/// D_index
/// @endverbatim
///
/// @see Cell.getNormalToPCAPlane()
/// @see Pressure3D::Constant
/// @see Pressure3D
/// @note used to be called VertexFromCellPlaneSpatial
///
class Spatial : public BaseReaction {
   private:
    double Kpow_;

   public:
    Spatial(std::vector<double> &paraValue,
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
/// @brief Same as Pressure3D::Constant but with the Force magnitude
/// dependent on a molecular concentration
///
/// @details This class is the same 'pressure' from inside update as
/// Pressure3D::Constant with the difference that the strength of the
/// resulting force depends on the cellular concentration of a
/// molecule described with a Hill formalism (increasing with conc).
/// The force is described by:
/// \f[
/// F = \frac{A_{cell}}{N_{vertex}}(p_{0} +
/// p_{1} \frac{C^{p_{3}}}{p_{2}^{p_{3}}+C^{p_{3}}}
/// \f]
/// where the area factor \f$A_{cell}\f$ is present if \f$p_{4}=1\f$,
/// and C is the molecular concentration.
/// In a model file the reaction is defined as
/// @verbatim
/// Pressure3D::ConcentrationHill 5 1 1
/// P_const P_add K_H n_H A_flag
/// C_index
/// @endverbatim
///
/// @see Cell.getNormalToPCAPlane()
/// @see Pressure3D::Constant
/// @see Pressure3D
/// @note used to be called VertexFromCellPlaneConcentrationHill
///
class ConcentrationHill : public BaseReaction {
   private:
    double Kpow_;

   public:
    ConcentrationHill(std::vector<double> &paraValue,
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
/// @brief PCA-plane independent forces normal to cell face
///
/// @details In a model file the reaction is defined by:
/// @verbatim
/// Pressure3D::Normalized 2 0
/// P A_flag
/// @endverbatim
///
/// @see Pressure3D
/// @see Pressure3D::Constant
/// @note used to be called VertexFromCellPlaneNormalized
///
class Normalized : public BaseReaction {
   public:
    Normalized(std::vector<double> &paraValue,
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
/// @brief PCA-plane independent forces to cell face with spatial component
///
/// @details In a model file the reaction is defined by:
/// @verbatim
/// Pressure3D::NormalizedSpatial 2 0
/// P_const P_add K_Hill n_Hill A_flag
/// @endverbatim
///
/// @see Pressure3D
/// @see Pressure3D::Spatial
/// @note used to be called VertexFromCellPlaneNormalizedSpatial
///
class NormalizedSpatial : public BaseReaction {
   private:
    double Kpow_;

   public:
    NormalizedSpatial(std::vector<double> &paraValue,
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
/// @brief PCA-plane independent force in SphereCylinder 'direction'
///
/// @details In a model file the reaction is defined by:
/// @verbatim
/// Pressure3D::SphereCylinder 2 0
/// P A_flag
/// @endverbatim
///
/// @see Pressure3D
/// @note used to be called VertexFromCellPlaneSphereCylinder
///
class SphereCylinder : public BaseReaction {
   public:
    SphereCylinder(std::vector<double> &paraValue,
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
/// @brief PCA-indpendent Force in SphereCylinder direction with a
/// concentration input
///
/// @details In a model file the reaction is defined by:
/// @verbatim
/// Pressure3D::SphereCylinderConcentrationHill 5 1 1
/// P_const P_add K_Hill n_Hiil A_flag
/// Conc_index
/// @endverbatim
///
/// @see Pressure3D
/// @see Pressure3D::SphereCylinder
/// @see Pressure3D::ConcentrationHill
/// @note used to be called VertexFromCellPlaneSphereCylinderConcentrationHill
///
class SphereCylinderConcentrationHill : public BaseReaction {
   public:
    SphereCylinderConcentrationHill(std::vector<double> &paraValue,
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
/// @brief Updates vertices from a 'pressure' term defined to act in the cell
/// normal direction, for triangular cells only.
///
/// @details This function calculates the area of a cell and then distribute a
/// force 'outwards' among the cell vertices. It relies on that the cells are
/// triangular. A cell contributes to a vertex update with
/// @f[ \frac{dx_{i}}{dt} = p_{0} A n_{i} / N_{vertex} @f]
/// where @f$p_{0}@f$ is a 'pressure' parameter, A is the cell area @f$n_{i}@f$
/// is the cell normal component and @f$N_{vertex}@f$ is the number of vertices
/// for the cell. An additional parameter @f$p_{2}@f$ can be used to not include
/// the area factor if set to zero (normally it should be set to 1). As an
/// indication for equilibrium state in case of elastic deformations the
/// function calculates the volume between template and Z=Z0 plane.this volume
/// is not used in calculations so Z0 value can be choosen arbitrarily.
/// In a model file the reaction is defined as
/// @verbatim
/// Pressure3D::Triangular 3 0
/// P A_flag Z0
/// @endverbatim
/// or
/// @verbatim
/// Pressure3D::Triangular 6 0
/// P A_flag Z0 V_eq Vf_decrease Pf_increase
/// @endverbatim
///
/// @see Pressure3D
/// @see Pressure3D::CenterTriangulation::Linear
/// @note used to be called VertexFromCellPlaneTriangular
///
class Triangular : public BaseReaction {
   public:
    Triangular(std::vector<double> &paraValue,
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
namespace CenterTriangulation {
///
/// @brief Updates vertices from a 'pressure' term defined to act in the normal
/// direction  to  triangular elements of the cell.
///
/// @details The force for each triangular
/// element is calculated according to  the  element's  current area and in the
/// direction of normal  to each  triangular  element. The force is distributed
/// equally on  nodes (including the centeral node). The  pressure  is  applied
/// increasingly (... linear in a given time span).
/// It does  not rely  on  PCA  plane  in contrast with VertexFromCellPlane and
/// VertexFromCellPlaneLinear.
/// A cell contributes to a vertex update with
/// @f[ \frac{dx_{i}}{dt} = p_{0} A n_{i} / 3 @f]
/// where @f$p_{0}@f$ is a 'pressure' parameter, A is the triangular element
/// area @f$n_{i}@f$ is the triangular element normal vector . An additional
/// parameter @f$p_{2}@f$ can be used to not include the area factor if set to
/// zero (normally it should be set to 1).
/// In a model file the reaction is defined as
/// @verbatim
/// Pressure3D::CenterTriangulation::Linear 3 1 1
/// P (pressure)
/// Area_flag (0: no area, 1: area, 2: area and pressure only in z direction)
/// deltaT (time to increase force from 0 to P)
///
/// InternalVarStartIndex
/// @endverbatim
///
/// @see CalculatePCAPlane
/// @see Pressure3D
/// @see Pressure3D::Linear
/// @note used to be called VertexFromCellPlaneCenterTriangulation
///
class Linear : public BaseReaction {
   private:
    double timeFactor1, timeFactor2;

   public:
    Linear(std::vector<double> &paraValue,
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
}  // end namespace CenterTriangulation
}  // namespace Pressure3D

#endif
