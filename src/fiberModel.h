//
// Filename     : fiberModel.h
// Description  : Classes describing updates of mechanical properties, e.g. Young's modulii given input
//                from elsewhere
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//
#ifndef FIBERMODEL_H
#define FIBERMODEL_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

/// 
/// @brief Reactions updating the mechanical properties, e.g. Young's modulii, for anisotropic material using
/// input, e.g. stress anisotropy, from elsewhere.
///
/// @details These reactions uses input read from a cell variable to update the mechanical properties of a 2D
/// wall. It has been used e.g. when stress anisotropy feeds back to the material properties of the wall, and
/// updates the maximal (and minimal) values for the Young's modulii in different directions. The reactions only
/// reads scalar values and update scalars and do not deal with the directions themselves, which needs to be done
/// elsewhere. These reactions have for example been used to update the Young's modulii in Bozorg et al (2014).
///
/// @see Reactions in directionUpdate.h, e.g. UpdateMTDirectionEquilibrium.  
///
namespace FiberModel {

  ///
  /// @brief Updates Young modulus of cells given an anisotropy measure. Update is within update function
  /// (between derivs steps) and based on linear or nonlinear Fiber_model, i.e. feedback on mechanical
  /// anisotropy from input (e.g. stress) anisotropy. The update can be direct following the feedback
  /// equation or gradually responding to the input.
  ///
  /// @details This reaction updates the mechanical properties (Young's modulii) in response to an anisotropy
  /// (e.g. stress or strain) input. It uses (reads) an anisotropy value, @f$a@f$, as input (e.g. stress or
  /// strain anisotropy) that should be calculated elsewhere e.g. saved from VertexFromTRBS* functions.
  /// This is used to calculate an mechanical (Young's modulus) 'target' anisotropy value using a linear
  /// (p2=0):
  /// @f[ Y_L^{target} = Y_M + 0.5(1.0+a)Y_F@f]
  /// or Hill-form version (p2=1 or 2):
  /// @f[ Y_L^{target} = Y_M + 0.5(1.0 + \frac{a^n}{(1-a)^{n}K^{n} + a^{n}})Y_F}@f]
  /// @f$Y_M=p_5@f$ and @f$Y_F=p_6@f$ are given parameters for a matrix and fiber contribution to the Young's
  /// modulii in different directions and should be the same as for the mechanical model used. When the Hill
  /// version is used @f$K=p_3,n=p_4@f$ needs to be given (otherwise not used). When @f$p_2=0/1@f$, a smooth
  /// change of the Young' modulus is implemented via an Euler step
  /// @f[ Y_L^{new} = Y_L + k_{rate} (Y_L^{target}-Y_L)\Delta t@f]
  /// where @f$k_{rate}=p_0@f$ is the update rate, and @f$Y_L@f$ is read and stored in a given cell index.
  /// When @f$p_2=2@f$, the update is direct and @f$Y_L@f$ is read and the updated value is stored at another
  /// cell index (FiberLindex). In addition, the update will not happen if a cell variable (velocity_index)
  /// is larger than a threshold given as @f$p_1@f$ (in the case of @f$p_2=0,1@f$), or if @f$p_0=0.0@f$.
  /// This can e.g. be used to only update close to mechanical equilibrium by having a threshold on a vertex
  /// velocity calculated and stored by the Calculate::VertexVelocity function. Before the
  /// simulation starts, the anisotropy and @f$Y_L@f$ can be initiated: if @f$p_7=0@f$ no initiation; =1
  /// anisotropy is set to zero in all cells and @f$Y_L=Y_M+0.5Y_F@f$; =2 @f$Y_L@f$ is set to its target value
  /// given the anisotropy in the cell.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// FiberModel::General 8 3 1 1[2] 1
  ///  k_rate
  ///  velocity_threshold
  ///  linear-hill_flag (0=linear/1=Hill/2=Hill_direct)
  ///  k_hill
  ///  n_hill
  ///  Y_matrix
  ///  Y_fiber
  ///  initiate_flag (0=no initiation/1=initiate with isotropic/2=initiate with anisotropy from aniso_index)
  ///
  ///  anisotropy_index
  ///  Young_Longitudinal_index [FiberL_index (only given/needed if linear-hill_flag=2)]
  ///  velocity_index (e.g. from "Calculate::VertexVelocity"; only used if linear-hill_flag=0/1)
  /// @endverbatim
  ///
  /// @see Equations 8 and 9 and Fig. 6 in Bozorg et al (2014) PLoS Comp Biol for Hill version.
  /// @see VertexFromTRBScenterTriangulationMT and similar functions for calculating strain/stress anisotropy.
  /// @see Calculate::VertexVelocity (for example) to provide vertex velocity values that can block update.
  ///
  class General : public BaseReaction {
  
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
    General(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue );  
    ///
    /// @brief Reaction initiation applied before simulation starts
    ///
    /// @see BaseReaction::initiate(Tissue &T,...)
    ///
    void initiate(Tissue &T,
		  DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
		  DataMatrix &cellDerivs,
		  DataMatrix &wallDerivs,
		  DataMatrix &vertexDerivs );  
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///
    void derivs(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData, 
		double h); 
  };

  ///
  /// @brief Updates Young modulus of cells within "update" based on linear or nonlinear Fiber_model.
  ///
  /// @details
  /// ...
  /// 
  /// In a model file the reaction is defined as
  /// @verbatim
  /// 
  /// FiberModel::Deposition 6 1 6
  /// 
  ///  k_rate
  ///  initial uniform fiber 
  ///  velocity threshold
  ///  init flag
  ///  k_hill
  ///  n_hill
  ///
  ///  anisotropy index.
  ///  misses stress index
  ///  cell area index
  ///  Young Fiber index
  ///  Young Longitudinal Fiber index
  ///  velocity index
  ///
  ///
  /// @endverbatim
  ///
  class Deposition : public BaseReaction {
  
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
    Deposition(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue );  
    ///
    /// @brief Reaction initiation applied before simulation starts
    ///
    /// @see BaseReaction::initiate(Tissue &T,...)
    ///
    void initiate(Tissue &T,
		  DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
		  DataMatrix &cellDerivs,
		  DataMatrix &wallDerivs,
		  DataMatrix &vertexDerivs );  
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///
    void derivs(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData, 
		double h); 
  };

} // namespace FiberModel

#endif //FIBERMODEL_H
