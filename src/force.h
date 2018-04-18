//
// Filename     : force.h
// Description  : Classes describing mechanical updates from external forces
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : April 2018
// Revision     : $Id:$
//
#ifndef FORCE_H
#define FORCE_H

#include <cmath>
#include "baseReaction.h"
#include "tissue.h"

///
/// @brief Reactions describing updates coming from forces applied in different scenarios.
///
/// @details External forces can be applied to the tissue to represent walls (not allowed to pass)
/// or attraction to specific geometries to help simulations generate specific shapes etc.
/// These reaction are collected within this namespace.
///
namespace Force {

  ///
  /// @brief Applies a force perpendicular to a defined wall of infinite size, defined for a specific coordinate
  ///
  /// @details A spring force in a perpendicular direction to a specific coordinate is applied.
  /// The update is given by
  /// @f[\frac{dx[/y/z]_i}{dt} -= p_0 (x_i-p_1) @f]
  /// if @f$x_i>p_1@f$, i.e. the vertex has crossed the 'wall', and for the coordinate specified.
  /// This is when @f$p_2=1@f$. When @f$p_2=-1@f$, the update is
  /// @f[\frac{dx[/y/z]_i}{dt} -= p_0 (x_i-p_1) @f]
  /// if @f$x_i<p_1@f$, i.e. crossing the 'wall' from the right.
  /// In a model file, the reaction is given by
  /// @verbatim
  /// Force::InfiniteWall 3 1 1
  /// K_force threshold(wall position) direction_flag
  /// position(coordinate)
  /// @endverbatim
  ///
  class InfiniteWall : public BaseReaction {
  public:
    InfiniteWall(std::vector<double> &paraValue,
		      std::vector<std::vector<size_t>> &indValue);
    
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };
  
  ///
  /// @brief Applies a force on epidermal vertices in direction (coordinate) given by index
  ///
  /// @details A force in a specific coordinate direction is applied to epidermal vertices. 
  /// It will update epidermal indices according to
  /// @f[\frac{dx[y,z]_i}{dt} += p_0 p_1 @f]
  /// where the coordinate is given as first index and @f$p_0@f$ is the force and @f$p_1@f$
  /// is a flag +/-1 to set the direction (1->outwards, -1->inwards). 
  /// In a model file this is given by
  /// @verbatim
  /// Force::EpidermalCoordinate 2 1 1
  /// K_force direction(=1/-1)
  /// coordinate_index
  /// @endverbatim
  ///
  class EpidermalCoordinate : public BaseReaction {
  public:
    EpidermalCoordinate(std::vector<double> &paraValue,
			std::vector<std::vector<size_t>> &indValue);
    
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };
  
  ///
  /// @brief Applies a force on epidermal vertices in a radial direction
  ///
  /// @details A force in a radial direction is applied to epidermal vertices. 
  /// It will update epidermal indices according to
  /// @f[\frac{dx[y,z]_i}{dt} -= p_0 \frac{x_i}{R_i} @f]
  /// where each coordinate is updated @f$p_0@f$ is the force and @f$R_i@f$
  /// is the radius (distance to origo) of the vertex. Note the negative sign. i.e
  /// it moves vertices inwards if the parameter is positive (and outwards if negative).
  /// In a model file this is given by
  /// @verbatim
  /// Force::EpidermalRadial 1 0
  /// K_force
  /// @endverbatim
  ///
  class EpidermalRadial : public BaseReaction {
  public:
    EpidermalRadial(std::vector<double> &paraValue,
		    std::vector<std::vector<size_t>> &indValue);
    
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };
} // end namespace Force

#endif
