// Filename     : hypocotyl3D.h
// Description  : Classes describing special reactions for a 3D hypocotyl growth simulation
// Author(s)    : Behruz Bozorg
//              : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : August 2018
// Revision     : $Id:$
//
#ifndef HYPOCOTYL3D_H
#define HYPOCOTYL3D_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

///
/// @brief A collection of reaction used specifically for simulating growth of a 3D Hypocotyl
///
/// @details These reactions are specific cases of mechanical and growth reactions that where
/// altered to run specific simulations for a growing 3D hypocotyl used in:
/// Bou Daher F, Chen Y, Bozorg B, Heywood Clough J, Jönsson H, Braybrook SA (2018)
/// Anisotropic growth is achieved through the additive mechanical effect of material 
/// anisotropy and elastic asymmetry. eLife
///
namespace Hypocotyl3D {

  ///
  /// @brief scales the template by a factor via Initiate 
  /// copies vectors from one index to another in the cell vector
  /// ( 4 component after the indices will be copied)
  ///
  /// @details In the model file the reaction is defined as:
  /// @verbatim
  /// Hypocotyl3D::limitZdis 0 1 2
  /// copy_from_index
  /// copy_to_index
  /// @endverbatim 
  /// 
  class limitZdis : public BaseReaction
  {
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
    limitZdis(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
    
    ///
    /// @brief Initiation function for this reaction class
    ///
    /// @see BaseReaction::initiate(Compartment &compartment,size_t species,...)
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
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
    
  private:
    
    std::vector< size_t > topVertices; //holds the vertex indices at the top
    std::vector< size_t > bottomVertices; //holds the vertex indices at the bottom
  };

  ///
  /// @brief This reaction is currently using ad hoc additions/changes within the code and should
  /// only be used by an expert, i.e. Behruz? or when model file been provided by an expert.
  ///
  /// In a model file the reaction is defined as (parameters changed?) 
  ///
  /// @verbatim
  /// Hypocotyl3D::StrainTRBS 2 3 1 1 3
  /// k_growth s_threshold  
  ///    
  /// L_ij-index
  /// InternalVarStartIndex 
  ///
  /// strain1_index
  /// strain2_index
  /// strain_vector_index
  /// @endverbatim
  ///
  /// (strain value and direction calculated and updated from other (mechanical) reactions).
  ///
  class StrainTRBS : public BaseReaction {
    
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
    StrainTRBS(std::vector<double> &paraValue, 
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
    
    void update(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		double h );
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




} // end namespace Hypocotyl3D

#endif //HYPOCOTYL3D_H
