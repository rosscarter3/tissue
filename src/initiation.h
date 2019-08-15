//
// Filename     : initiation.h
// Description  : Classes describing initiation rules for a tissue simulation
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : August 2019
// Revision     : $Id:$
//
#ifndef INITIATION_H
#define INITIATION_H

#include <cmath>

#include "myRandom.h"
#include "baseReaction.h"

///
/// @brief Collection of reactions describing tissue initiation rules (e.g. setting concentration variables)
///
/// @details Within this namespace rules for setting initial values of variables are included.
/// These reactions do not create the geometry, but rather only set variable values. For initiating geometries, see
/// tissue/tools/createInit/*, and more developed https://gitlab.com/slcu/teamhj/laura/tissue-init-generator.
///
/// Also, some initiations like CenterTriangulation, SisterVertex, edge resting lengths are not included in this namespace.
namespace Initiation {

  ///
  /// @brief This rule sets a cell variable to 1 with probability @f$p_{0}@f$, 0 otherwise
  ///
  /// @details For each cell the given cell index variable is set to
  /// restricting variable given by
  ///
  /// @f[ y_{ij}} = 1 if R<p_0 
  /// y_{ij} = 0 if R>=p_0 @f]
  ///
  /// where p_0 is the probability for each cell i to get value 1 for molecule j.
  /// In the model file the reaction is defined by
  /// @verbatim
  /// Initiation::RandomBoolean 1 1 1
  /// p
  /// cell_var_index
  /// @endverbatim
  ///
  class RandomBoolean : public BaseReaction {
  
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
    RandomBoolean(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > &indValue );
    
    ///
    /// @brief Derivative function for this reaction class (does nothing for this class).
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
    
    /// @brief Initiation made before simulation
    ///
    /// @see BaseReaction::initiate()
    ///
    void initiate(Tissue &T,
                  DataMatrix &cellData,
                  DataMatrix &walldata,
                  DataMatrix &vertexData,
                  DataMatrix &cellderivs,
                  DataMatrix &wallderivs,
                  DataMatrix &vertexDerivs );        
  };
} // end of namespace Initiation
#endif
