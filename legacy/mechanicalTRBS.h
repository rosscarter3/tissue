//
// Filename     : mechanicalTRBS.h
// Description  : Classes describing mechanical updates for wall springs
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2007
// Revision     : $Id:$
//
#ifndef MECHANICALTRBS_H
#define MECHANICALTRBS_H

#include "tissue.h"
#include "baseReaction.h"
#include <cmath>

///
/// @brief Mechanical models using the Triangular Biquadratic Spring (TRBS) for
/// triangular plates
///
/// @details These reactions update the vertex positions from the mechanical
/// feedback on 2D triangular elements. The theory of the mechanical model comes
/// from H. Delingette, Triangular springs for modelling non-linear membranes,
/// IEEE Trans Vis Comput Graph 14, 329-41 (2008)
///
/// @see https://gitlab.com/slcu/teamhj/publications/bozorg_etal_2014 (and
/// bozorg_etal_2016) for examples
///
namespace TRBS {

// Some functions for general triangle manipulations
///
/// @brief Calculates the triangular area from edges using Heron's formula
///
double AreaFromEdges(std::vector<double> &edgeLength);

///
/// @brief Calculates the cosine of the triangle angles from the edge lengths
///
void CosFromEdges(std::vector<double> &restingLength,
                  std::vector<double> &cosAngle);
///
/// @brief Calculates the Cotanges of the angles of the triangle from the
/// cosines
///
void CotanFromCos(std::vector<double> &cosAngle,
                  std::vector<double> &cotanAngle);
///
/// @brief Calculates the tensile and angular stiffness from cotan vector
///
void Stiffness(double lambda, double mio, double areaFactor,
               std::vector<double> &cotan, std::vector<double> &tensilStiffness,
               std::vector<double> &angularStiffness);
///
/// @brief Calculates the BiQuadratic strain
///
void BiQuadraticStrain(std::vector<double> &length,
                       std::vector<double> &restingLength,
                       std::vector<double> &Delta);
///
/// @brief Shape vector matrix in local coordinate system from shape factors (P
/// or Q = X)
///
/// @details This matrix is the inverse of coordinate matrix. Only first two
/// elements are used in calculations i.e. shapeVector[3][2] although
/// implemented as 3x3.
///
void ShapeVector(double Xa, double Xb, double Xc,
                 std::vector<std::vector<double> > &shapeVector);
///
/// @brief Calculates Forces from isotropic contribution of the material
/// (youngT,poissonT)
///
/// The calculation uses positions, tensile and angular stiffnesses and Delta
/// and adds to the Force matrix.
///
void AddForceIsotropic(std::vector<std::vector<double> > &position,
                       std::vector<double> &tensileStiffness,
                       std::vector<double> &angularStiffness,
                       std::vector<double> &Delta,
                       std::vector<std::vector<double> > &Force);
///
/// @brief Calculates the rotation matrix from positions
///
void RotationMatrix(std::vector<std::vector<double> > &position,
                    std::vector<std::vector<double> > &rotation);

///
/// @brief Rotates inVector into outVector using the rotation matrix
///
void Rotate(std::vector<double> &inVector,
            std::vector<std::vector<double> > &rotation,
            std::vector<double> &outVector);

///
/// @brief Rotates a tensor and saves the output in the input matrix
///
void RotTensorRot(std::vector<std::vector<double> > &rotation,
                  std::vector<std::vector<double> > &Tensor);

///
/// @brief Normalise a vector to length 1.0
///
void Normalise(std::vector<double> &inVector);

///
/// @brief Extract eigenvectors from a matrix using Jacobi method
///
/// @details Eigenvectors for the matrix are extracted using a Jacobi method.
/// The eigenvectors are normalised and stored in columns in the eigenVectors
/// output
///
void GetEigenVectors(std::vector<std::vector<double> > &eigenVectors,
                     std::vector<std::vector<double> > &inMatrix,
                     double epsilon = 1.0e-06);

///
/// @brief Sets the values of a matrix to the Identity matrix
///
void SetIdentity(std::vector<std::vector<double> > &matrix);

namespace CenterTriangulation {} // end namespace CenterTriangulation

} // end namespace TRBS

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
                 std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
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
  VertexFromTRBScenterTriangulation(
      std::vector<double> &paraValue,
      std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
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
  VertexFromTRBScenterTriangulationConcentrationHill(
      std::vector<double> &paraValue,
      std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
};

///
/// @brief Triangular spring model with anisotropy for plates (2D walls)
/// assuming triangular walls/cells.
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
                   std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
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
/// VertexFromTRBScenterTriangulationMT 11/13 2 11/12 1
///
/// Y_matrix
/// Y_fiber
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1/2/10) or < 0: if heterogeneity is considered, the value for this
///                           flag is the scale factor for heterogeneity)
/// neighborWeight
/// max_stress(if 0 absolute stress anisotropy is calculated)
/// plane-strain/stress-flag
/// MT-angle
/// MT-feedback-flag
/// unused parameter
/// (loosening_index)
///
/// L_ij-index
/// MT_cellIndex
/// strainAnisotropy-Index
/// stressAnisotropy-Index
/// areaRatioIndex
/// youngT-index
/// anisoEnergyIndex
/// youngL-index/heterpogeneity_index
/// MTstressIndex
/// stressTensorIndex
/// normalVectorIndex
/// (loosening_K) (loosening_n)
///
/// InternalVarStartIndex
///
/// or
///
/// VertexFromTRBScenterTriangulationMT 11/13 4 11/12 1 0/1/2/3 0/1/2
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
/// (loosening_index)
///
/// L_ij-index
/// MT_cellIndex
/// strainAnisotropy-Index
/// stressAnisotropy-Index
/// areaRatioIndex
/// youngT-index
/// anisoEnergyIndex
/// youngL-index
/// MTstressIndex
/// stressTensorIndex
/// normalVectorIndex
/// (loosening_K) (loosening_n)
///
/// InternalVarStartIndex
///
/// optional indices for storing strain(0: no strain,
///                                     1: strain,
///                                     2: strain/perpendicular strain,
///                                     3: strain/perpendicular strain/2nd
///                                     strain)
/// optional indices for storing stress(0: no stress,
///                                     1: stress,
///                                     2: stress/2nd stress)
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions,
/// strain/stress values will be stored after (3) components of vectors.
/// The value for perpendicular strain is maximal strain value.

class VertexFromTRBScenterTriangulationMT : public BaseReaction {
private:
  double timeC = 0;
  bool lengthout = false;

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
  VertexFromTRBScenterTriangulationMT(
      std::vector<double> &paraValue,
      std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, double h);

  void printState(Tissue *T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, std::ostream &os);
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
///                                     3: strain/perpendicular strain/2nd
///                                     strain)
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
  VertexFromTRLScenterTriangulationMT(
      std::vector<double> &paraValue,
      std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

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
/// or
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

class VertexFromTRBScenterTriangulationConcentrationHillMT
    : public BaseReaction {

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
  VertexFromTRBScenterTriangulationConcentrationHillMT(
      std::vector<double> &paraValue,
      std::vector<std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///
  void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
              DataMatrix &vertexData, DataMatrix &cellDerivs,
              DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
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
  VertexFromTRBScenterTriangulationMTOpt(
      std::vector<double> &paraValue,
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

#endif
