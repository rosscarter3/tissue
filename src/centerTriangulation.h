//
// Filename     : centerTriangulation.h
// Description  : Classes describing updates for tissues with cells storing a central point (and internal edges)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2012
// Revision     : $Id:$
//
#ifndef CENTERTRIANGULATION_H
#define CENTERTRIANGULATION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief center triangulation namespace
///


namespace CenterTriangulation
{

  ///
  /// @brief Initiates a center triangulation, either from information in the tissue or from scratch
  ///
  /// @details This reaction does not update the tissue. It only initiates a central triangulation of cells
  /// and add variables to cellData at initiation. This is mainly used for triangular biquadratic
  /// spring models. If no parameters are provided, the initiation will use information form the tissue
  /// if available (need to specify flag -centerTri_init at runtime to read CT information), or create
  /// from scratch if the tissue does not have the information. If one parameter (p0=overrideFlag) is
  /// provided, and its value is set to 1, the CT will be created from scratch (even if the tissue has one)
  /// A second parameter (p1=doubleEdgeFlag) can be provided, and if set to one a double edge CT will be
  /// created from information in the tissue if available (and p0=0) or from scratch. The double edge
  /// CT is used e.g. for growth rules connected to TRBS.
  /// One variable index is provided for compability with an old version (and useful to make sure all
  /// reactions using CT have the smae value for this variable), and
  /// should represents the end of the cellData vector (cell(i).numVariable()). 
  ///
  /// In a model file the reaction is given by
  /// @verbatim
  /// CeterTriangulation::Initiate 0/1/2 1 1
  /// [override_flag] #=0 follow init tissue if CT available, =1 always create from scratch
  /// [doubleEdge_flag] #=0 single edge,=1 double edge
  ///
  /// InternalVarStartIndex
  /// @endverbatim
  ///
  /// @note The order of reactions may be of importance when using this reaction, 
  /// e.g. if new vertices are added on the walls (then this should be done afterwards).
  /// @see Tissue::readInitCenterTri()
  /// @see Cell::centerPosition()
  ///
  class Initiate : public BaseReaction {
  
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
    Initiate(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > &indValue );
    
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
} // end namespace CenterTriangulation
  
#endif //CENTERTRIANGULATION_H
