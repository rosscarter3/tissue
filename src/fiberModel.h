//
// Filename     : fiberModel.h
// Description  : Classes describing updates of mechanical properties, e.g.
// Young's modulii given input
//                from elsewhere
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//
#ifndef FIBERMODEL_H
#define FIBERMODEL_H

#include <cmath>

#include "tissue.h"
#include "baseReaction.h"

///
/// @brief Reactions updating the mechanical properties, e.g. Young's modulii,
/// for anisotropic material using input, e.g. stress anisotropy, from
/// elsewhere.
///
/// @details These reactions uses input read from a cell variable to update the
/// mechanical properties of a 2D wall. It has been used e.g. when stress
/// anisotropy feeds back to the material properties of the wall, and updates
/// the maximal (and minimal) values for the Young's modulii in different
/// directions. The reactions only reads scalar values and update scalars and do
/// not deal with the directions themselves, which needs to be done elsewhere.
/// These reactions have for example been used to update the Young's modulii in
/// Bozorg et al (2014). The original version including all options is named
/// General, and simplified implementations for Linear and Hill functions for
/// the feedback is also provided. The Keywork Equilibrium in the reaction names
/// indicate a continuous update towards a target (equilibrium) value.
///
/// @see Reactions in directionUpdate.h, e.g. UpdateMTDirectionEquilibrium.
///
namespace FiberModel {

///
/// @brief Calculation of isotropic Young Modulii used for initiation of
/// material to isotropic.
/// @f$ Y_L = Y_T = Y_M + 0.5 Y_F @f$
///
double isotropicFunction(double YoungMatrix, double YoungFiber);

///
/// @brief Calculation of max Young's modulus used by FiberModel::Linear*
/// reactions.
/// @f$ Y_L = Y_M + 0.5 (1.0+a) Y_F @f$
///
double linearFunction(double input, double YoungMatrix, double YoungFiber);

///
/// @brief Calculation of max Young's modulus used by FiberModel::Hill*
/// reactions.
/// @f$ Y_L = Y_M + 0.5(1.0 + \frac{a^n}{(1-a)^{n}K^{n} + a^{n}})Y_F}@f$
///
double hillFunction(double input, double YoungMatrix, double YoungFiber);

///
/// @brief Calculation of Transverse (min) Young's modulus given Longitudinal,
/// Matrix and Fiber values
/// @f$ Y_T = 2 Y_M + Y_F - Y_L @f$
///
double transverseFromLongitudinal(double youngL, double youngMatrix,
                                  double youngFiber);

///
/// @brief Updates Young modulus of cells given an anisotropy measure. Update is
/// within update function (between derivs steps) and based on linear or
/// nonlinear Fiber_model, i.e. feedback on mechanical anisotropy from input
/// (e.g. stress) anisotropy. The update can be direct following the feedback
/// equation or gradually responding to the input.
///
/// @details This reaction updates the mechanical properties (Young's modulii)
/// in response to an anisotropy (e.g. stress or strain) input. It uses (reads)
/// an anisotropy value, @f$a@f$, as input (e.g. stress or strain anisotropy)
/// that should be calculated elsewhere e.g. saved from VertexFromTRBS*
/// functions. This is used to calculate an mechanical (Young's modulus)
/// 'target' anisotropy value using a linear (p2=0):
/// @f[ Y_L^{target} = Y_M + 0.5(1.0+a)Y_F@f]
/// or Hill-form version (p2=1 or 2):
/// @f[ Y_L^{target} = Y_M + 0.5(1.0 + \frac{a^n}{(1-a)^{n}K^{n} +
/// a^{n}})Y_F}@f]
/// @f$Y_M=p_5@f$ and @f$Y_F=p_6@f$ are given parameters for a matrix and fiber
/// contribution to the Young's modulii in different directions and should be
/// the same as for the mechanical model used. When the Hill version is used
/// @f$K=p_3,n=p_4@f$ needs to be given (otherwise not used). When
/// @f$p_2=0/1@f$, a smooth change of the Young' modulus is implemented via an
/// Euler step
/// @f[ Y_L^{new} = Y_L + k_{rate} (Y_L^{target}-Y_L)\Delta t@f]
/// where @f$k_{rate}=p_0@f$ is the update rate, and @f$Y_L@f$ is read and
/// stored in a given cell index. When @f$p_2=2@f$, the update is direct and
/// @f$Y_L@f$ is read and the updated value is stored at another cell index
/// (FiberLindex). In addition, the update will not happen if a cell variable
/// (velocity_index) is larger than a threshold given as @f$p_1@f$ (in the case
/// of @f$p_2=0,1@f$), or if @f$p_0=0.0@f$. This can e.g. be used to only update
/// close to mechanical equilibrium by having a threshold on a vertex velocity
/// calculated and stored by the Calculate::VertexVelocity function. Before the
/// simulation starts, the anisotropy and @f$Y_L@f$ can be initiated: if
/// @f$p_7=0@f$ no initiation; =1 anisotropy is set to zero in all cells and
/// @f$Y_L=Y_M+0.5Y_F@f$; =2 @f$Y_L@f$ is set to its target value given the
/// anisotropy in the cell. In a model file the reaction is defined as
/// @verbatim
/// FiberModel::General 8 3 1 1[2] 1
///  k_rate
///  velocity_threshold
///  linear-hill_flag (0=linear/1=Hill/2=Hill_direct)
///  k_hill
///  n_hill
///  Y_matrix
///  Y_fiber
///  initiate_flag (0=no initiation/1=initiate with isotropic/2=initiate with
///  anisotropy from aniso_index)
///
///  anisotropy_index
///  Young_Longitudinal_index [FiberL_index (only given/needed if
///  linear-hill_flag=2)] velocity_index (e.g. from "Calculate::VertexVelocity";
///  only used if linear-hill_flag=0/1)
/// @endverbatim
/// @note This version is the original FiberModel reaction.
/// @see Equations 8 and 9 and Fig. 6 in Bozorg et al (2014) PLoS Comp Biol for
/// Hill version. https://gitlab.com/slcu/teamhj/publications
/// @see VertexFromTRBScenterTriangulationMT and similar functions for
/// calculating strain/stress anisotropy.
/// @see Calculate::VertexVelocity (for example) to provide vertex velocity
/// values that can block update.
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
          std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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

///
/// @brief Updates Young's modulii of cells given an anisotropy input using a
/// linear function. Update is within update function and updates the mechanics
/// in a single step.
///
/// @details This reaction updates the mechanical properties (Young's modulii)
/// in response to an anisotropy (or other) input. It uses (reads) an anisotropy
/// value, @f$a@f$, as input (e.g. stress or strain anisotropy) that should be
/// calculated elsewhere e.g. saved from Mechanical::TRBS::* functions. This is
/// used to calculate an mechanical (Young's modulus) 'target' anisotropy value
/// using a linear function
/// @f[ Y_L = Y_M + 0.5 (1.0+a) Y_F @f]
/// where @f$Y_L@f$ is the Longitudinal (maximal) Young's modulus, which is
/// updated here, @f$Y_M=p_0,Y_F=p_1@f$ are the Matrix (isotropic) and Fiber
/// (anisotropic if @f$a@f$ not equal to 0) material parameters. These should be
/// the same as for the mechanical model used. In a model file, the most basic
/// form of the reaction is defined as
/// @verbatim
/// FiberModel::Linear 2 2 1 1
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
///
/// anisotropy_index # cell anisotropy variable input
/// Young_L_index    # cell variable for saving Y_L @endverbatim
/// or if initiation of and possible saving of the transverse (minimal) Young's
/// modulus is used
/// @verbatim
/// FiberModel::Linear 2/3 2 1 1[/2]
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index)
///
/// anisotropy_index              # cell anisotropy variable input
/// Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
/// @endverbatim where a second index can be used to store/update also the
/// transverse value of the material given by
/// @f[ Y_T = Y_M + 0.5 (1.0-a) Y_F @f]
/// Before the simulation starts, the anisotropy and @f$Y_L@f$ can be initiated:
/// if @f$p_2=0@f$ or not given - no initiation; @f$p_2=1@f$ - anisotropy is set
/// to zero in all cells and @f$Y_L=Y_M+0.5Y_F@f$;
/// @f$p_2=2@f$ -  @f$Y_L@f$ is set to its target value given the anisotropy
/// value in the cell.
///
/// In addition, the update will not happen if a cell variable (velocity_index)
/// is larger than a threshold given as an additional parameter @f$p_3@f$.
/// This can e.g. be used to only update close to mechanical equilibrium by
/// having a threshold on a vertex velocity calculated and stored by the
/// Calculate::VertexVelocity function. In a model file the reaction is defined
/// as
/// @verbatim
/// FiberModel::Linear 2/3/4 2/3 1 1[2] [1]
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index) [Threshold] # Threshold in read velocity for
/// when updates not done
///
///  anisotropy_index              # cell anisotropy variable input
///  Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
///  [velocity_index]              # e.g. from "Calculate::VertexVelocity"
///  @endverbatim
/// @see FiberModel::LinearEquilibrium for same function but with dynamic update
/// towards target value.
/// @see VertexFromTRBScenterTriangulationMT and similar functions for
/// calculating strain/stress anisotropy.
/// @see Calculate::VertexVelocity (for example) to provide vertex velocity
/// values that can block update.
///
class Linear : public BaseReaction {

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
  Linear(std::vector<double> &paraValue,
         std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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

///
/// @brief Updates Young's modulii of cells given an anisotropy input using a
/// linear function. Update is within update function and the material
/// parameters are updated continuosly towards target value.
///
/// @details This reaction updates the mechanical properties (Young's modulii)
/// in response to an anisotropy (or other) input. It uses (reads) an anisotropy
/// value, @f$a@f$, as input (e.g. stress or strain anisotropy) that should be
/// calculated elsewhere e.g. saved from Mechanical::TRBS::* functions. This is
/// used to calculate an mechanical (Young's modulus) 'target' anisotropy value
/// using a linear function
/// @f[ Y_L^{target} = Y_M + 0.5 (1.0+a) Y_F @f]
/// where @f$Y_L@f$ is the Longitudinal (maximal) Young's modulus, which is
/// updated here, @f$Y_M=p_1,Y_F=p_2@f$ are the Matrix (isotropic) and Fiber
/// (anisotropic if @f$a@f$ not equal to 0) material parameters. These should be
/// the same as for the mechanical model used. A smooth change of the Young'
/// modulus is implemented via an Euler step
/// @f[ Y_L^{new} = Y_L + k_{rate} (Y_L^{target}-Y_L)\Delta t@f]
/// where @f$k_{rate}=p_0@f$ is the update rate, and @f$Y_L@f$ is read and
/// stored in the given cell index. In a model file, the most basic form of the
/// reaction is defined as
/// @verbatim
/// FiberModel::LinearEquilibrium 3 2 1 1
/// k_rate      # Rate of update of value towards target
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
///
/// anisotropy_index # cell anisotropy variable input
/// Young_L_index    # cell variable for saving Y_L @endverbatim
/// or, if initiation and possibly transverse (minimal) Young's modulus is saved
/// @verbatim
/// FiberModel::LinearEquilibrium 3/4 2 1 1[/2]
/// k_rate      # Rate of update of value towards target
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index)
///
/// anisotropy_index              # cell anisotropy variable input
/// Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
/// @endverbatim where a second index can be used to store/update also the
/// transverse value of the material given by
/// @f[ Y_T^{target} = Y_M + 0.5 (1.0-a) Y_F @f]
/// Before the simulation starts, the anisotropy and @f$Y_L@f$ can be initiated:
/// if @f$p_3=0@f$ or not given - no initiation; @f$p_3=1@f$ - anisotropy is set
/// to zero in all cells and @f$Y_L=Y_M+0.5Y_F@f$;
/// @f$p_3=2@f$ -  @f$Y_L@f$ is set to its target value given the anisotropy
/// value in the cell.
///
/// In addition, the update will not happen if a cell variable (velocity_index)
/// is larger than a threshold given as an additional parameter @f$p_4@f$.
/// This can e.g. be used to only update close to mechanical equilibrium by
/// having a threshold on a vertex velocity calculated and stored by the
/// Calculate::VertexVelocity function. In a model file the reaction is defined
/// as
/// @verbatim
/// FiberModel::LinearEquilibrium 3/4/5 2/3 1 1[2] [1]
/// k_rate      # Rate of update of value towards target
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index) [Threshold] # Threshold in read velocity for
/// when updates not done
///
/// anisotropy_index              # cell anisotropy variable input
/// Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
/// [velocity_index]              # e.g. from "Calculate::VertexVelocity"
/// @endverbatim
/// @see FiberModel::LinearEquilibrium for same function but with dynamic update
/// towards target value.
/// @see VertexFromTRBScenterTriangulationMT and similar functions for
/// calculating strain/stress anisotropy.
/// @see Calculate::VertexVelocity (for example) to provide vertex velocity
/// values that can block update.
///
class LinearEquilibrium : public BaseReaction {

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
  LinearEquilibrium(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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

///
/// @brief Updates Young's modulii of cells given an anisotropy input using a
/// Hill function. Update is within update function and updates the mechanics in
/// a single step.
///
/// @details This reaction updates the mechanical properties (Young's modulii)
/// in response to an anisotropy (or other) input. It uses (reads) an anisotropy
/// value, @f$a@f$, as input (e.g. stress or strain anisotropy) that should be
/// calculated elsewhere e.g. saved from Mechanical::TRBS::* functions. This is
/// used to calculate an mechanical material (Young's modulus) 'target'
/// anisotropy value using a Hill function
/// @f[ Y_L = Y_M + 0.5(1.0 + \frac{a^n}{(1-a)^{n}K^{n} + a^{n}})Y_F}@f]
/// where @f$Y_L@f$ is the Longitudinal (maximal) Young's modulus, which is
/// updated here, @f$Y_M=p_0,Y_F=p_1@f$ are the Matrix (isotropic) and Fiber
/// (anisotropic if @f$a@f$ not equal to 0) material parameters.
/// @f$K=p_2,n=p_3@f$ are the Hill constant and coeffiecient, respectively.
/// These four parameters should be the same as for the mechanical model used.
/// In a model file, the most basic form of the reaction is defined as
/// @verbatim
/// FiberModel::Hill 4 2 1 1
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// K_hill      # Hill constant in feedback function
/// n_hill      # Hill coefficient
///
/// anisotropy_index # cell anisotropy variable input
/// Young_L_index    # cell variable for saving Y_L @endverbatim
/// or with added initiation and potential saving of transverse/minimal Young's
/// modulus
/// @verbatim
/// FiberModel::Hill 4/5 2 1 1[/2]
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// K_hill      # Hill constant in feedback function
/// n_hill      # Hill coefficient
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index)
///
/// anisotropy_index              # cell anisotropy variable input
/// Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
/// @endverbatim where a second index can be used to store/update also the
/// transverse value of the material given by
/// @f[ Y_T = Y_M + 0.5(1.0 - \frac{a^n}{(1-a)^{n}K^{n} + a^{n}})Y_F}@f]
/// Before the simulation starts, the anisotropy and @f$Y_L@f$ can be initiated:
/// if @f$p_4=0@f$ or not given - no initiation; @f$p_4=1@f$ - anisotropy is set
/// to zero in all cells and @f$Y_L=Y_M+0.5Y_F@f$;
/// @f$p_4=2@f$ - @f$Y_L@f$ is set to its target value given the anisotropy
/// value in the cell.
///
/// In addition, the update will not happen if a cell variable (velocity_index)
/// is larger than a threshold given as an additional parameter @f$p_5@f$.
/// This can e.g. be used to only update close to mechanical equilibrium by
/// having a threshold on a vertex velocity calculated and stored by the
/// Calculate::VertexVelocity function. In a model file the reaction is defined
/// as
/// @verbatim
/// FiberModel::Hill 2/3/4 2/3 1 1[2] [1]
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// K_hill      # Hill constant in feedback function
/// n_hill      # Hill coefficient
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index) [Threshold] # Threshold in read velocity for
/// when updates not done
///
///  anisotropy_index              # cell anisotropy variable input
///  Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
///  [velocity_index]              # e.g. from "Calculate::VertexVelocity"
///  @endverbatim
/// @see FiberModel::HillEquilibrium for same function but with dynamic update
/// towards target value.
/// @see VertexFromTRBScenterTriangulationMT and similar functions for
/// calculating strain/stress anisotropy.
/// @see Calculate::VertexVelocity (for example) to provide vertex velocity
/// values that can block update.
///
class Hill : public BaseReaction {

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
  Hill(std::vector<double> &paraValue,
       std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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

///
/// @brief Updates Young's modulii of cells given an anisotropy input using a
/// Hill feedback function. Update is within update function and the material
/// parameters are updated continuosly towards a target value.
///
/// @details This reaction updates the mechanical properties (Young's modulii)
/// in response to an anisotropy (or other) input. It uses (reads) an anisotropy
/// value, @f$a@f$, as input (e.g. stress or strain anisotropy) that should be
/// calculated elsewhere e.g. saved from Mechanical::TRBS::* functions. This is
/// used to calculate an mechanical (Young's modulus) 'target' anisotropy value
/// using a linear function
/// @f[ Y_L^{target} = Y_M + 0.5(1.0 + \frac{a^n}{(1-a)^{n}K^{n} +
/// a^{n}})Y_F}@f] where @f$Y_L@f$ is the Longitudinal (maximal) Young's
/// modulus, which is updated here, @f$Y_M=p_1,Y_F=p_2@f$ are the Matrix
/// (isotropic) and Fiber (anisotropic if @f$a@f$ not equal to 0) material
/// parameters. These should be the same as for the mechanical model used. A
/// smooth change of the Young' modulus is implemented via an Euler step
/// @f[ Y_L^{new} = Y_L + k_{rate} (Y_L^{target}-Y_L)\Delta t@f]
/// where @f$k_{rate}=p_0@f$ is the update rate, and @f$Y_L@f$ is read and
/// stored in the given cell index. In a model file, this basic form of the
/// reaction is defined as
/// @verbatim
/// FiberModel::HillEquilibrium 5/6 2 1 1[/2]
/// k_rate      # Rate of update of value towards target
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// K_hill      # Hill constant for feedback function
/// n_hill      # Hill coefficient
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index)
///
/// anisotropy_index              # cell anisotropy variable input
/// Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
/// @endverbatim where a second index can be used to store/update also the
/// transverse value of the material given by
/// @f[ Y_T^{target} = Y_M + 0.5(1.0 - \frac{a^n}{(1-a)^{n}K^{n} +
/// a^{n}})Y_F}@f] Before the simulation starts, the anisotropy and @f$Y_L@f$
/// can be initiated: if @f$p_5=0@f$ or not given - no initiation; @f$p_5=1@f$ -
/// anisotropy is set to zero in all cells and @f$Y_L=Y_M+0.5Y_F@f$;
/// @f$p_5=2@f$ -  @f$Y_L@f$ is set to its target value given the anisotropy
/// value in the cell.
///
/// In addition, the update will not happen if a cell variable (velocity_index)
/// is larger than a threshold given as an additional parameter @f$p_6@f$.
/// This can e.g. be used to only update close to mechanical equilibrium by
/// having a threshold on a vertex velocity calculated and stored by the
/// Calculate::VertexVelocity function. In a model file the reaction is defined
/// as
/// @verbatim
/// FiberModel::HillEquilibrium 5/6/7 2/3 1 1[2] [1]
/// k_rate      # Rate of update of value towards target
/// Y_matrix    # isotropic contribution to material
/// Y_fiber     # anisotropic contribution to materia
/// K_hill      # Hill constant for feedback function
/// n_hill      # Hill coefficient
/// [init_flag] # (0=no initiation, 1=initiate with isotropic, 2=initiate with
/// anisotropy from aniso_index) [Threshold] # Threshold in read velocity
/// variable for when updates not done
///
/// anisotropy_index              # cell anisotropy variable input
/// Young_L_index [Young_T_index] # cell variable for saving Y_L [Y_T]
/// [velocity_index]              # e.g. from "Calculate::VertexVelocity"
/// @endverbatim
/// @see FiberModel::Hill for same function but with direct update to target
/// value.
/// @see VertexFromTRBScenterTriangulationMT and similar functions for
/// calculating strain/stress anisotropy.
/// @see Calculate::VertexVelocity (for example) to provide vertex velocity
/// values that can block update.
///
class HillEquilibrium : public BaseReaction {

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
  HillEquilibrium(std::vector<double> &paraValue,
                  std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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

///
/// @brief Updates Young modulus of cells within "update" based on linear or
/// nonlinear Fiber_model.
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
             std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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

} // namespace FiberModel

#endif // FIBERMODEL_H
