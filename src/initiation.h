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
  /// @f[ y_{ij} = 1 if R<p_0 y_{ij} = 0 if R>p_0 @f]
  ///
  /// where p_0 is the probability for each cell i to get value 1 for molecule j.
  /// seed is an integer to seed the random number. A @f$\Delta T@f$ parameter can be given
  /// and then the random initiation is redone with a period (but where only ones can be given to
  /// cells currently off/zero). If set to 0, the rule will not repeat.
  /// In the model file the reaction is defined by
  /// @verbatim
  /// Initiation::RandomBoolean 2[3] 1 1
  /// p seed [Delta_T]
  /// cell_var_index
  /// @endverbatim
  ///
  class RandomBoolean : public BaseReaction {

    public:

      double localTime_=0.0; //variable to update time using h in update function
      double nextTime_=0.0; //variable updated to hold the next time it will launch the procedure
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

      ///
      /// @see BaseReaction::update()
      ///
      void update(Tissue &T,
          DataMatrix &cellData,
          DataMatrix &walldata,
          DataMatrix &vertexData,
          double h);
  };

    ///
  /// @brief This rule sets a cell variable to a random value at initiation
  ///
  /// @details For each cell the given cell index variable is set to
  /// restricting variable given by
  ///
  /// @f[ y_{ij} = p_[0] * R @f]
  ///
  /// where p_0 is multiplied with a random number 2f$ R \in [0:1] @f$.
  /// seed is an integer to seed the random number.
  /// In the model file the reaction is defined by
  /// @verbatim
  /// Initiation::Random 2 1 1
  /// p seed
  /// cell_var_index
  /// @endverbatim
  ///
  class Random : public BaseReaction {

    public:

      double localTime_=0.0; //variable to update time using h in update function
      double nextTime_=0.0; //variable updated to hold the next time it will launch the procedure
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
      Random(std::vector<double> &paraValue, 
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

      ///
      /// @see BaseReaction::update()
      ///
      void update(Tissue &T,
          DataMatrix &cellData,
          DataMatrix &walldata,
          DataMatrix &vertexData,
          double h);
  };

  ///
  /// @brief This rule sets a cell variable to 1 with probability @f$p_{0}@f$, 0 otherwise with a bias disallowing setting cells where
  /// more than @$N_t@$ neighbour cells are on
  ///
  /// @details For each cell the given cell index variable is set to
  /// restricting variable given by
  ///
  /// @f[ y_{ij} = 1 if R<p_0 and N_{on}<p_1, else y_{ij} = 0 @f]
  ///
  /// where p_0 is the probability for each cell i to get value 1 for molecule j.
  /// probability will be 0 if the number of neighbours that are 'on' ar more than @$N_{t}@$ or more.
  /// seed is an integer to seed the random number. A @f$\Delta T@f$ parameter can be given
  /// and then the random initiation is redone with a period (but where only ones can be given to
  /// off cells). If set to 0, it will not repeat.
  /// In the model file the reaction is defined by
  /// @verbatim
  /// Initiation::RandomBoolean 2[3] 1 1
  /// p N_t seed [\DeltaT]
  /// cell_var_index
  /// @endverbatim
  ///
  class RandomBooleanBiased : public BaseReaction {

    public:

      double localTime_=0.0; //variable to update time using h in update function
      double nextTime_=0.0; //variable updated to hold the next time it will launch the procedure
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
      RandomBooleanBiased(std::vector<double> &paraValue, 
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

      ///
      /// @see BaseReaction::update()
      ///
      void update(Tissue &T,
          DataMatrix &cellData,
          DataMatrix &walldata,
          DataMatrix &vertexData,
          double h);
  };

  ///
  /// @brief Reaction that initiates a cell variable to the cell size
  ///
  /// @details Sets a cell variable to the size of the cell (area in 2D is assumed). If a flag parameter
  /// is given the variable is updated using the update function. The cell size is calculated using the function
  /// cell.calculateVolume().
  /// In a model file the reaction is defined as:
  /// @verbatim
  /// Initiate::FaceArea2D 0/1 1 1
  /// [flag] # 0 for no update and 1 for updating
  /// varIndex
  /// @endverbatim
  ///
  class FaceArea2D : public BaseReaction {

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
      FaceArea2D(std::vector<double> &paraValue, 
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

      ///
      /// @see BaseReaction::update()
      ///
      void update(Tissue &T,
          DataMatrix &cellData,
          DataMatrix &walldata,
          DataMatrix &vertexData,
          double h);
  };

  ///
  /// @brief Set individual cell labels for cells by adding integer values to a cell variable
  ///
  /// @details Defines a cell label by initiating each cell with an integer label. Since it is not
  /// updated, the label can be used to follow lineages.
  /// In a model file the reaction is defined as:
  /// @verbatim
  /// Initiation::CellLabel 0 1 1
  /// varIndex
  /// @endverbatim
  ///
  class CellLabel : public BaseReaction {

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
      CellLabel(std::vector<double> &paraValue, 
          std::vector< std::vector<size_t> > &indValue );

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

  };

} // end of namespace Initiation
#endif
