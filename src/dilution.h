//
// Filename     : dilution.h
// Description  : Classes describing dilution of moleculecular concentrations due to growth
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//

#ifndef DILUTION_H
#define DILUTION_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

///
/// @brief Reactions describing dilution of moleculecular concentrations due to growth
///
/// @details When cells are growing, the effect of dilution might be of importance. These reactions use the change
/// in 2D Cell Volumes to add a dilution contribution to cell variables.
///
/// @note Some of these reactions use the derivative contribution of growth and hence needs to be placed
/// after any reactions updating the vertex positions to get the correct values.
///
namespace Dilution {

    ///
  /// @brief Updates 'concentration' variables according to volume changes from
  /// derivatives of the vertex positions.
  ///
  /// @details The dilution of 'concentration' variables, C_i are calculated according to
  ///
  /// @f[ \frac{C_i}{dt} = - \frac{C_i}{V} \frac{dV}{dt} @f]
  ///
  /// where V is the volume, and dV/dt is calculated from vertex position
  /// derivatives. In a model file the reaction is defined as
  /// @verbatim
  /// DilutionFromVertexDerivs 0 1 k
  /// c1_index [c2_index...ck_index] @endverbatim
  /// @note Since this function uses the derivatives of the vertex positions it
  /// needs to be applied after all other derivatives applied to the vertices.
  ///
  class FromVertexDerivs : public BaseReaction
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
    FromVertexDerivs(std::vector<double> &paraValue,
                             std::vector< std::vector<size_t> > &indValue);
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
                DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T,
		       DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       DataMatrix &cellDerivs,
		       DataMatrix &wallDerivs,
		       DataMatrix &vertexDerivs,
		       DataMatrix &sdydtCell,
		       DataMatrix &sdydtWall,
		       DataMatrix &sdydtVertex );
  };

}

#endif // DILUTION_H
