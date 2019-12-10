//
// Filename     : mechanicalTRBS.h
// Description  : Classes describing mechanical updates for wall springs
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2007
// Revision     : $Id:$
//
#ifndef MECHANICALTRBS_H
#define MECHANICALTRBS_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangular faces (cells).
///
/// @details The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
/// @verbatim
/// VertexFromTRBS 2 1 1
/// Y_modulus P_coeff
/// L_ij-index
/// @endverbatim
///
class VertexFromTRBS : public BaseReaction {
  
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
  VertexFromTRBS(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue );  
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
};

///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
///
/// @details The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
/// @verbatim
/// VertexFromTRBScenterTriangulation 2 2 1 1
/// Y_modulus P_coeff
/// L_ij-index
/// InternalVarStartIndex
/// or
/// VertexFromTRBScenterTriangulation 2 4 1 1 1/0 1/0
/// Y_modulus P_coeff
/// L_ij-index
/// InternalVarStartIndex
/// Optional index for storing strain
/// Optional index for storing stress
/// @endverbatim
///
class VertexFromTRBScenterTriangulation : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulation(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );  
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
};

///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell
/// concentration(auxin)-dependent Young's modulus.
///
/// @details The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
/// @verbatim
/// VertexFromTRBScenterTriangulationConcentrationHill 5 2 2 1
/// Y_modulus_min Y_modulus_max P_coeff K_hill n_hill
/// L_ij-index  concentration-index
/// InternalVarStartIndex
/// @endverbatim
///
class VertexFromTRBScenterTriangulationConcentrationHill : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulationConcentrationHill(std::vector<double> &paraValue, 
						     std::vector< std::vector<size_t> > 
						     &indValue );  
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
};

///
/// @brief Triangular spring model with anisotropy for plates (2D walls) assuming
/// triangular walls/cells.
///
/// @details The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008), wich is developed into an
/// anisotropic model.
///
/// In a model file the reaction is defined as
/// @verbatim
/// VertexFromTRBSMT 10 1 10
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1) 
/// neighborWeight 
/// max stress for stress aniso 
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag 
/// 
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index
/// stressAnisotropy-Index 
/// areaRatioIndex 
/// isoEnergyIndex
/// anisoEnergyIndex 
/// YoungL-index
/// MTstress
/// start index for stress tensor
///
/// or
///
/// VertexFromTRBSMT 10 3 10 0/1/2/3 0/1/2
///
/// Y_matrix
/// Y_fiber
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1)
/// neighborWeight
/// unusedparameter
/// plane-strain/stress-flag
/// MT-angle
/// MT-feedback-flag
/// 
/// L_ij-index
/// MT_cellIndex
/// strainAnisotropy-Index
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex
/// anisoEnergyIndex
/// YoungL-index
/// MTstress
/// start index for stress tensor
///
/// optional index for storing strain(0: no strain,
///                                   1: strain, 
///                                   2: strain/perpendicular strain, 
///                                   3: strain/perpendicular strain/2nd strain)
///
/// optional index for storing stress(0: no stress, 
///                                   1: stress, 
///                                   2: stress/2nd stress)
///
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after  (3) components of vectors.  
/// The value for perpendicular strain is maximal strain value.
///
class VertexFromTRBSMT : public BaseReaction {
  
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
  VertexFromTRBSMT(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue );  
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
};



///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBScenterTriangulationMT 11 2 11 1 
/// 
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1/2 or <0: if heterogeneity is considered, the value for this 
///                     flag is the scale factor for heterogeneity) 
/// neighborWeight 
/// max_stress(if 0 absolute stress anisotropy is calculated)
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag 
/// unused parameter 
///
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index 
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex 
/// anisoEnergyIndex 
/// youngL-index/heterpogeneity_index 
/// MTstressIndex 
/// stressTensorIndex 
/// normalVectorIndex
///
/// InternalVarStartIndex
/// 
/// or
/// 
/// VertexFromTRBScenterTriangulationMT 11 4 11 1 0/1/2/3 0/1/2
///
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long  
/// Poisson_Trans 
/// MF_flag(0/1) 
/// neighborWeight 
/// unusedparameter 
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag
/// unused parameter
/// 
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index 
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex 
/// anisoEnergyIndex 
/// youngL-index 
/// MTstressIndex 
/// stressTensorIndex 
/// normalVectorIndex
///
/// InternalVarStartIndex
///
/// optional indices for storing strain(0: no strain, 
///                                     1: strain, 
///                                     2: strain/perpendicular strain, 
///                                     3: strain/perpendicular strain/2nd strain)
/// optional indices for storing stress(0: no stress, 
///                                     1: stress, 
///                                     2: stress/2nd stress)
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after (3) components of vectors.
/// The value for perpendicular strain is maximal strain value.  

class VertexFromTRBScenterTriangulationMT : public BaseReaction {
private:
  
  double timeC=0;
  bool lengthout=false;

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
  VertexFromTRBScenterTriangulationMT(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );  
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
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void update(Tissue &T,
              DataMatrix &cellData,
              DataMatrix &wallData,
              DataMatrix &vertexData, 
              double h); 
  
  void printState(Tissue *T,
                  DataMatrix &cellData,
                  DataMatrix &wallData,
                  DataMatrix &vertexData, 
                  std::ostream &os);
    

};


///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRLScenterTriangulationMT 11 2 11 1 
/// 
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1/2 or <0: if heterogeneity is considered, the value for this 
///                     flag is the scale factor for heterogeneity) 
/// neighborWeight 
/// max_stress(if 0 absolute stress anisotropy is calculated)
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag 
/// unused parameter 
///
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index 
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex 
/// anisoEnergyIndex 
/// youngL-index/heterpogeneity_index 
/// MTstressIndex 
/// stressTensorIndex 
/// normalVectorIndex
///
/// InternalVarStartIndex
/// 
/// or
/// 
/// VertexFromTRLScenterTriangulationMT 11 4 11 1 0/1/2/3 0/1/2
///
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long  
/// Poisson_Trans 
/// MF_flag(0/1) 
/// neighborWeight 
/// unusedparameter 
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag
/// unused parameter
/// 
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index 
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex 
/// anisoEnergyIndex 
/// youngL-index 
/// MTstressIndex 
/// stressTensorIndex 
/// normalVectorIndex
///
/// InternalVarStartIndex
///
/// optional indices for storing strain(0: no strain, 
///                                     1: strain, 
///                                     2: strain/perpendicular strain, 
///                                     3: strain/perpendicular strain/2nd strain)
/// optional indices for storing stress(0: no stress, 
///                                     1: stress, 
///                                     2: stress/2nd stress)
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after (3) components of vectors.
/// The value for perpendicular strain is maximal strain value.  

class VertexFromTRLScenterTriangulationMT : public BaseReaction {
private:
  
  
  
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
  VertexFromTRLScenterTriangulationMT(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );  
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
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  
// void initiate(Tissue &T,
// 		DataMatrix &cellData,
// 		DataMatrix &wallData,
// 		DataMatrix &vertexData,
// 		DataMatrix &cellDerivs,
// 		DataMatrix &wallDerivs,
// 		DataMatrix &vertexDerivs );  


};



///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
/// concentration(auxin)-dependent Young's modulus.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBScenterTriangulationConcentrationHillMT 8 2 3 1
/// Y_modulus_Longitudinal_min Y_modulus_Longitudinal_max P_coeff_Longitudinal 
/// Y_modulus_Transverse_min Y_modulus_Transverse_max P_coeff_Transverse  
/// K_hill n_hill
/// L_ij-index  concentration-index MT_cellIndex
/// InternalVarStartIndex
///or
/// VertexFromTRBScenterTriangulationConcentrationHillMT 8 6 3 1 1/0 1/0 1/0 1/0
/// Y_modulus_Longitudinal_min Y_modulus_Longitudinal_max P_coeff_Longitudinal 
/// Y_modulus_Transverse_min Y_modulus_Transverse_max P_coeff_Transverse  
/// K_hill n_hill
/// L_ij-index  concentration-index MT_cellIndex
/// InternalVarStartIndex
/// optional index for storing strain
/// optional index for storing 2nd strain
/// optional index for storing stress
/// optional index for storing 2nd stress
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after  2(3) components of vectors.  

class VertexFromTRBScenterTriangulationConcentrationHillMT : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulationConcentrationHillMT(std::vector<double> &paraValue, 
						     std::vector< std::vector<size_t> > 
						     &indValue );  
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
};

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
/// cell index (FiberLindex). In addition, the update will not happen if a cell variable (velocity_index,
/// e.g. a velocity calculated and stored by the UpdateMTdirectionEquilibrium function) is larger than a
/// threshold given as @f$p_1@f$ (in the case of @f$p_2=0,1@f$), or if @f$p_0=0.0@f$. Before the
/// simulation starts, the anisotropy and @f$Y_L@f$ can be initiated: if @f$p_7=0@f$ no initiation; =1
/// anisotropy is set to zero in all cells and @f$Y_L=Y_M+0.5Y_F@f$; =2 @f$Y_L@f$ is set to its target value
/// given the anisotropy in the cell.
/// In a model file the reaction is defined as
/// @verbatim
/// FiberModel 8 3 1 1[2] 1
/// 
///  k_rate
///  velocity_threshold
///  linear-hill_flag [0=linear/1=Hill/2=Hill_direct]
///  k_hill
///  n_hill
///  Y_matrix
///  Y_fiber
///  initiate_flag (0=no initiation/1=initiate with isotropic/2=initiate with anisotropy from aniso_index)
///
///  anisotropy_index
///  Young_Longitudinal_index [FiberL_index (only if linear-hill_flag=2)]
///  velocity_index (e.g. from "UpdateMTdirectionEquilibrium" only if linear-hill_flag=0/1)
/// @endverbatim
///
/// @see Equations 8 and 9 and Fig. 6 in Bozorg et al (2014) PLoS Comp Biol for Hill version.
/// @see VertexFromTRBScenterTriangulationMT and similar functions for calculating strain/stress anisotropy.
/// @see UpdateMTdirectionEquilibrium (for example) to provide velocity values that can block update.
class FiberModel : public BaseReaction {
  
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
  FiberModel(std::vector<double> &paraValue, 
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
/// FiberDeposition 6 1 6
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
class FiberDeposition : public BaseReaction {
  
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
  FiberDeposition(std::vector<double> &paraValue, 
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
/// @brief mechanical TRBS model with pressure via energy minimization
///
/// @details
/// ...
/// 
/// In a model file the reaction is defined as
/// @verbatim
///
/// VertexFromTRBScenterTriangulationMTOpt
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
class VertexFromTRBScenterTriangulationMTOpt : public BaseReaction {
private:
  std::vector<std::vector<std::vector<double> > > stateVector;
  double totalEnergy;
  double mechIsEn, mechAnEn, PEn;
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
  VertexFromTRBScenterTriangulationMTOpt(std::vector<double> &paraValue, 
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


#endif
