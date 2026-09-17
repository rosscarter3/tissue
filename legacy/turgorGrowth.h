//
// Filename     : turgorGrowth.h
// Description  : Classes describing growth updates based on active updates of turgor pressure
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//
#ifndef TURGORGROWTH_H
#define TURGORGROWTH_H

#include <cmath>

#include "baseReaction.h"
#include "tissue.h"

///
/// @brief Reactions updates of turgor pressure and growth resulting from dynamic turgor pressure
/// and water movement (using water potantial etc)
///
/// @details Reactions updating turgor pressure e.g. via water potential and water volume growth reactions
/// following this.
/// Vertex positions are updated to generate increase (or possibly decrease) of 2D cell volumes, and turgor is
/// calculated by water potential and water movement can be calculated. Additional cell
/// variables are updated and used for dynamic variables.
///
/// @note These reactions have never been used in publications and should be checked before used in a
/// real application.
///
namespace TurgorGrowth {

///
/// @brief Updates the water volume variable given osmotic and turgor potentials
///
/// @details This function uses a constant osmotic potential and calculates the turgor
/// potential to calculate water intake into the cell according to
/// @f[ \frac{V_w}{dt} = p_0 A (p_1-p_2 T) @f]
/// where @f$V_w@f$ is the water volume, T is the turgor, @f$p_0@f$ is the update rate, @f$p_1@f$ is the
/// osmotic potential and @f$p_2@f$ is an scaling factor. Also @f$p_3@f$ = denyShrinkFlag
/// and @f$p_4@f$ = allowNegTurgorFlag can be set to restrict the behavior.
/// The turgor, T, is calculated as @f$T=(V_w-V)/V@f$, the difference between the current water volume
/// and cell (face) volume. In a model file the reaction is defined as
/// @verbatim
/// TurgorGrowth::WaterVolume 5 1/2 1 [1]
/// k_p
/// P_max
/// k_pp
/// denyShrink_flag
/// allowNegativeTurgor_flag
///
/// WaterVolume
/// [Turgor_save] @endverbatim
/// @note Used to be called WaterVolumeFrom Turgor
class WaterVolume : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    /// @param indValue vector of vectors with variable indices
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    WaterVolume(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> > &indValue);
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
};

}  // namespace TurgorGrowth

#endif  // TURGORGROWTH_H
