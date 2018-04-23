//
// Filename     : calculate.h
// Description  : Classes for calculating and updating variables
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : April 2018
// Revision     : $Id:$
//
#ifndef CALCULATE_H
#define CALCULATE_H

#include <cmath>
#include "baseReaction.h"
#include "tissue.h"

namespace Calculate {
  ///
  /// @brief Calculates abs(cos(...)) of angle between two 3d vectors
  /// (starting from given indices) in cellData vector and stores it in the given
  /// index in cellData vector.
  ///
  /// @details This reaction uses no parameters. In a model file the
  /// reaction is defined as
  /// @verbatim
  /// Calculate::AngleVectors 0 2 2 1
  /// start-index(1st vector)   start-index(2nd vector)
  /// store-index(angle-deg)
  /// @endverbatim
  ///
  class AngleVectors : public BaseReaction {
  private:
    DataMatrix vertexDataRest;
  public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AngleVectors(std::vector<double> &paraValue,
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
  /// @brief Calculates abs(cos(...)) of angle between a 3d vector
  /// (starting from given indices) in cellData vector and XY plane and stores it
  /// in the given index in cellData vector.
  ///
  /// @details This reaction uses no parameters. In a model file the
  /// reaction is defined as:
  /// @verbatim
  /// Calculate::AngleVectorXYplane 0 2 1 1
  /// start-index(vector)
  /// store-index(angle-deg)
  /// @endverbatim
  ///
  class AngleVectorXYplane : public BaseReaction {
  private:
    DataMatrix vertexDataRest;
  public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AngleVectorXYplane(std::vector<double> &paraValue,
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
  /// @brief Calculates the angle between a 3d vector (starting from given
  /// indices) in cellData vector and a given axes (x,y,z).
  ///
  /// @details Uses one parameter for specifying the axis, and two variable
  /// indices. The first index specifies the start of the vector and the second
  /// where the angle is stored.
  /// In a model file the reaction is defined as:
  /// @verbatim
  /// Calculate::AngleVector 1 2 1 1
  /// axis_flag (0:X, 1:Y, 2:Z)
  /// start-index(the vector)
  /// store-index(angle-deg)
  /// @endverbatim
  ///
  /// @note Needs to be updated to be used for 3D and y or z axes.
  ///
  class AngleVector : public BaseReaction {
  private:
    DataMatrix vertexDataRest;
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
    AngleVector(std::vector<double> &paraValue,
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
  /// @brief Extracts the maximum velocity of vertices (from vertexDerivs)
  /// and stores it in a given index in cellData vector
  /// to check the closeness to mechanical equilibrium
  ///
  /// @details Uses no parameter
  /// The first index specifies the index for storage
  /// In a model file the reaction is defined as:
  /// @verbatim
  /// Calculate::MaxVelocity 0 1 1
  /// velocity-store-index(angle-deg)
  /// @endverbatim
  ///
  class MaxVelocity : public BaseReaction {
  private:
  DataMatrix vertexDataRest;
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
  MaxVelocity(std::vector<double> &paraValue,
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
  /// @brief Calculates change in the tissue volume and its time derivative
  /// and total Derivative and stores them in the given indices in cellData
  /// vector.
  ///
  /// @details This reaction calculates the total change in volume as the sum
  /// of movement of the vertices from time 0. It also calculates the change
  /// between time points, and also the sum of the positional derivatives. 
  /// In a model file the reaction is defined as
  /// @verbatim
  /// TemplateVolumeChange 0 1 6
  /// cell-index-VolumeChange       component-index-VolumeChange
  /// cell-index-deltaVolumeChange  component-index-deltaVolumeChange
  /// cell-index-totalDerivative    component-index-totalDerivative
  /// @endverbatim
  ///
  /// @note Used to be called TemplateVolumeChange (still allowed).
  /// @note Since these are 'global' variables the storage in specific
  /// components of the cellData matrix is somewhat ad hoc.
  ///
  class TissueVolumeChange : public BaseReaction {
  private:
    DataMatrix vertexDataRest;
    double VolumeChange;
    double deltaVolumeChange;
    double totalDerivative;    
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
    TissueVolumeChange(std::vector<double> &paraValue,
		       std::vector<std::vector<size_t>> &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
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
  
} // end namespace Calculate

#endif
