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
    InfiniteWall(std::vector<double> &paraValue,
		 std::vector<std::vector<size_t>> &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///    
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
    EpidermalCoordinate(std::vector<double> &paraValue,
			std::vector<std::vector<size_t>> &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };

  ///
  /// @brief Growth via vertex forces radially outwards
  ///
  /// The tissue grows from vertex movement radially outwards. The update is given by
  /// @f[ \frac{dr}{dt} = p_{0} @f] (if @f$ p_1=0 @f$) or
  /// @f[ \frac{dr}{dt} = p_{0} r @f] (if @f$ p_{1}=1 @f$)
  /// where @f$ p_{0} @f$ is the force/rate (@f$ k_{growth} @f$),
  /// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used (@f$ r_{pow} @f$).
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Force::Radial 2 0
  /// p_0 p_1
  /// @endverbatim
  ///
  /// @note Used to be named MoveVertexRadially (still allowed)
  ///
  class Radial : public BaseReaction {
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
    Radial(std::vector<double> &paraValue,
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
    /// @brief For this function it just returns the derivs (i.e. no noise added)
    ///
    /// @see BaseReaction::derivsWithAbs(Tissue &T,...)
    ///
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
  
  ///
  /// @brief Applies a force on epidermal vertices in a radial direction
  ///
  /// @details A force in a radial direction is applied to epidermal vertices. 
  /// It will update epidermal indices according to
  /// @f[\frac{dx[y]_i}{dt} += p_0 \frac{x_i}{R_i} @f]
  /// where x and y coordinates for vertex i are updated. @f$p_0@f$ is the force and @f$R_i@f$
  /// is the radius (distance to origo (x,y)) of the vertex.
  /// It moves vertices outwards if the parameter is positive (and inwards if negative).
  /// In a model file this is given by
  /// @verbatim
  /// Force::EpidermalRadial 1 0
  /// K_force
  /// @endverbatim
  ///
  /// @see Force::Radial
  ///
  class EpidermalRadial : public BaseReaction {
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
    EpidermalRadial(std::vector<double> &paraValue,
		    std::vector<std::vector<size_t>> &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };

  ///
  /// @brief Applies a force in a radial direction on vertices listed
  ///
  /// @details A force in a radial direction is applied to vertices listed. 
  /// It will update epidermal indices according to
  /// @f[\frac{dx[y,z]_i}{dt} += p_0 p_1 x_i @f]
  /// where each coordinate is updated @f$p_0@f$ is the force.
  /// Note, @f$p_1@f$ is a flag setting the direction i.e
  /// if 1, the reaction moves vertices outwards (and inwards if -1).
  /// In a model file this is given by
  /// @verbatim
  /// Force::IndexRadial 1 1 N
  /// K_force
  /// v_1 ... v_N
  /// @endverbatim
  ///
  /// @note Used to be named VertexForceOrigoFromIndex.
  /// @see Force::Radial
  /// @see Force::EpidermalRadial (similar update but vertices are selected to be updated if they are at boundary)
  /// @see Force::CellIndexRadial (similar but cell indices are listed to set vertices to be updated)
  ///
  class IndexRadial : public BaseReaction {
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
    IndexRadial(std::vector<double> &paraValue,
		std::vector<std::vector<size_t>> &indValue);  
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };

  ///
  /// @brief Applies a force in a radial direction on vertices of cells listed
  ///
  /// @details A force in a radial direction is applied to vertices of cells listed. 
  /// It will update epidermal indices according to
  /// @f[\frac{dx[y,z]_i}{dt} += p_0 p_1 x_i @f]
  /// where each coordinate is updated @f$p_0@f$ is the force.
  /// Note, @f$p_1@f$ is a flag setting the direction i.e
  /// if 1, the reaction moves vertices outwards (and inwards if -1).
  /// In a model file this is given by
  /// @verbatim
  /// Force::CellIndexRadial 1 1 N
  /// K_force direction_flag(=+/-1)
  /// c_1 ... c_N
  /// @endverbatim
  /// where the list of cells indices whose vertices should be updated.
  ///
  /// @note Used to be named CellForceOrigoFromIndex.
  /// @see Force::Radial
  /// @see Force::EpidermalRadial (similar update but vertices are selected to be updated if they are at boundary)
  /// @see Force::CellIndexRadial (similar but cell indices are listed to set vertices to be updated)
  ///
  class CellIndexRadial : public BaseReaction {
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
    CellIndexRadial(std::vector<double> &paraValue,
		    std::vector<std::vector<size_t>> &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };
  
  ///
  /// @brief Applies a constant axial force on vertices in two intervals in z (opposing force direction),
  /// used e.g on a cylinderical template to approximate internal tissue forces generating axial forces/growth.
  ///
  /// @details A constant force, F, is applied axially to regions in z in [z_0+a, z_0+a+d] (upward force) and
  /// [z_0-a, z_0-a-d] (downward) resembling forces generated by internal tissues in e.g.
  /// Hypocotyls that can provide forces for axial growth.
  /// In a model file the reaction is defined as:
  /// @verbatim
  /// Force::Axial 4 0
  /// z_0 a d F
  /// @endverbatim
  /// where z_0 is the 'center' of the tissue, a and d defines the two regions, and F is the magnitude of the F. 
  ///
  /// @note Used to be called VertexFromHypocotylGrowth
  /// @note For developers: A more advanced version (8 parmeters) generating axial and radial forces is commented out in force.h(.cc) 
  ///
  class Axial : public BaseReaction {
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
    Axial(std::vector<double> &paraValue,
	  std::vector<std::vector<size_t>> &indValue);    
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
		DataMatrix &vertexData, DataMatrix &cellDerivs,
		DataMatrix &wallDerivs, DataMatrix &vertexDerivs);
  };
  
} // end namespace Force

#endif
