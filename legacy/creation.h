//
// Filename     : creation.h
// Description  : Classes describing molecular production/creation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#ifndef CREATION_H
#define CREATION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

/// @brief Namespace collecting reactions for molecular production
/// @details The reactions describe production (creation) of molecules in
/// terms of concentrations. They may include advanced rules, while transcription
/// is taken care of by Grn.
/// @see Grn
namespace Creation {
  /// @brief In each cell a molecule is produced/created with a constant rate
  /// @details The variable update is for each cell given by 
  /// @f[ \frac{dc}{dt} = k_c @f]
  /// where @f$ k_c @f$ is a constant parameter, and @f$ c @f$ is the variable to be updated.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::Zero 1 1 1
  /// k_c
  /// c_index
  /// @endverbatim
  class Zero : public BaseReaction {    
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
    Zero(std::vector<double> &paraValue, 
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
    /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
    ///
    /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
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

  /// @brief In each cell a molecule is produced/created with a rate dependent on another molecule.
  /// @details The variable update is for each cell given by 
  /// @f[ \frac{dc}{dt} = - k_c X @f]
  /// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the variable to be updated,
  /// and @f$ X @f$ is the concentration of the production-dependent molecule.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::One 1 2 1 1
  /// k_c
  /// c_index
  /// X_index
  /// @endverbatim
  class One : public BaseReaction {    
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
    One(std::vector<double> &paraValue, 
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
    /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
    ///
    /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
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
  
  /// @brief In each cell a molecule is produced/created with a rate dependent on two other molecules.
  /// @details The variable update is for each cell given by 
  /// @f[ \frac{dc}{dt} = - k_c X Y @f]
  /// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the variable to be updated,
  /// and @f$ X Y @f$ are the concentrations of the production-dependent molecules.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::Two 1 2 1 2
  /// k_c
  /// c_index
  /// X_index
  /// Y_index
  /// @endverbatim
  class Two : public BaseReaction {
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
    Two(std::vector<double> &paraValue, 
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
  };

  /// @brief In each cell a molecule is produced/created with a rate dependent on three other molecules.
  /// @details The variable update is for each cell given by 
  /// \f[ \frac{dc}{dt} = - k_c X Y \f]
  /// where \f$ k_c \f$ is a constant parameter, \f$ c \f$ is the variable to be updated,
  /// and \f$ X Y Z \f$ are the concentrations of the production-dependent molecules.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::Three 1 2 1 3
  /// k_c
  /// c_index
  /// X_index
  /// Y_index
  /// Z_index
  /// @endverbatim
  class Three : public BaseReaction {
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
    Three(std::vector<double> &paraValue, 
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
  };

  /// @brief In each cell a molecule is produced/created with rate dependent on the distance of the cell from the center
  /// @details The variable update is for each cell given by ( SIGN= -1, production inside the sphere)
  /// @f[ \frac{dc}{dt} = V \frac{R^n}{r^n + R^n} @f]
  /// or (SIGN = +1, production outside the sphere),
  /// @f[ \frac{dc}{dt} = V \frac{r^n}{r^n + R^n} @f]
  /// where @f$ V, R, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to be updated
  /// and @f$ r @f$ the distance of the cell to the center of the template.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::SpatialSphere 4 1 1
  /// V R n SIGN
  /// c_index
  /// @endverbatim
  class SpatialSphere : public BaseReaction {
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
    SpatialSphere(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief In each cell a molecule is produced/created with rate dependent on the distance of the cell
  /// from the center of a cylinder along the z axis
  ///
  /// @details The variable update is for each cell given by ( SIGN= -1, production inside the sphere)
  /// @f[ \frac{dc}{dt} = V \frac{r^n + R^n}{R^n} @f]
  /// or (SIGN = +1, production outside the sphere),
  /// @f[ \frac{dc}{dt} = V \frac{r^n + R^n}{R^n} @f]
  /// where @f$ V, R, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to be updated
  /// and @f$ r @f$ the distance of the cell to the center of a cylinder along the z axis with radius @f$ R @f$.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::SpatialCylinder 4 1 1
  /// V R n SIGN
  /// c_index
  /// @endverbatim
  ///
  class SpatialCylinder : public BaseReaction {
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
    SpatialCylinder(std::vector<double> &paraValue, 
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
  };
  
  ///
  /// @brief In each cell a molecule is produced/created with rate dependent on the distance of the cell from a ring
  ///
  /// @details The variable update is for each cell given by ( SIGN= -1, production inside the ring)
  /// @f[ \frac{dc}{dt} = V \frac{K^n}{K^n + d^n} @f]
  /// or (SIGN = +1, production away from the ring),
  /// @f[ \frac{dc}{dt} = V \frac{d^n}{K^n + d^n} @f]
  /// where @f$ V, K, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to be updated,
  /// and @f$ d @f$ the distance of the cell to ring:
  /// @f[ d = |r_{cell} - r_{ring}| @f]
  /// In the above definition, @f$ r_{cell} @f$ is the distance of the cell to the center of the template,
  /// and @f$ r_{ring} @f$ is a constant parameter
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::SpatialRing 5 1 1
  /// V R r_ring n SIGN
  /// c_index
  /// @endverbatim
  ///
  class SpatialRing : public BaseReaction {
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
    SpatialRing(std::vector<double> &paraValue, 
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
  };
  
  ///
  /// @brief In each cell a molecule is produced/created with rate dependent on one of
  /// the coordinates of each cell
  ///
  /// @details The variable update is for each cell given by ( SIGN= -1, higher production for
  /// lower values of the coordinate)
  /// @f[ \frac{dc}{dt} = V \frac{X^n}{X^n + x^n} @f]
  /// or (SIGN = +1, higher production at larger coordinate values),
  /// @f[ \frac{dc}{dt} = V \frac{x^n}{x^n + X^n} @f]
  /// where @f$ V, X, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to
  /// be updated and @f$ x @f$ is the cell coordinate used to set the rate-dependent
  /// production.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::SpatialCoordinate 4 2 1 1
  /// V X n SIGN
  /// c_index
  /// x_index
  /// @endverbatim
  ///
  class SpatialCoordinate: public BaseReaction {
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
    SpatialCoordinate(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief In each cell a molecule is produced/created if it is on one side of
  /// a plane
  ///
  /// @details The variable update is for each cell given by ( SIGN= -1, production for
  /// lower values of the coordinate)
  /// @f[ \frac{dc}{dt} = V if x <= X, otherwise 0 @f]
  /// or (SIGN = +1, production at larger coordinate values),
  /// @f[ \frac{dc}{dt} = V if x >= X, otherwise 0 @f]
  /// where @f$ V, X,SIGN@f$ are constant parameters, @f$ c @f$ is the variable to
  /// be updated and @f$ x @f$ is the cell coordinate used to set the rate-dependent
  /// production.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::SpatialPlane 3 2 1 1
  /// V X 
  /// c_index
  /// x_index
  /// @endverbatim
  ///
  class SpatialPlane: public BaseReaction {
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
    SpatialPlane(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief In the cells with indices given in a list a molecule is created with constant rate.
  ///
  /// @details The variable update is for each cell given by
  /// @f[ \frac{dc}{dt} = k_c @f]
  /// if cell index is in the given list
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::FromList 1/2 2 1 n
  /// k_c
  /// [numberFlag] # if 1 const num instead of conc added 
  /// c_index
  /// a list of n cell indices 
  /// @endverbatim
  /// If the second parameter (numberFlag) is provided and the given value is one, the
  /// constant production is in number of molecules, i.e. production constant is divided by
  /// cell size.
  ///
  class FromList: public BaseReaction {
  private: 
    size_t proCells;
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
    FromList(std::vector<double> &paraValue, 
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
  };
  
  ///
  /// @brief In each cell a molecule is produced/created with a rate dependent on another molecule and cell volume.
  ///
  /// @details The variable update is for each cell given by 
  /// @f[ \frac{dc}{dt} = volume* k_c X @f]
  /// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the variable to be updated,
  /// and @f$ X @f$ is the concentration of the production-dependent molecule.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::OneGeometric 1 2 1 1
  /// k_c
  /// c_index
  /// X_index
  /// @endverbatim
  ///
  class OneGeometric : public BaseReaction {
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
    OneGeometric(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief Molecular creation oscillating with a sinusodial shape
  ///
  /// @details Creation reaction describing an oscillatory behavior described as 
  /// @f[\frac{dc}{dt} = p_{0}/2 (1+sin( \frac{(p_{1}*time+p_{2})}{2\pi}) @f]
  /// where p0 is the amplitude, p1 the period, p3 the phase and time (updated
  /// locally during simulation.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// Creation::Sinus 3 1
  /// amp period phase
  /// index
  /// @endverbatim
  ///
  class Sinus : public BaseReaction {
  private:
    double time_;
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
    /// @see BaseReaction::createReaction(std::vector<double>&,...)
    ///
    Sinus(std::vector<double> &paraValue, 
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
    ///
    /// @see BaseReaction::update()
    ///
    void update(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &walldata,
		DataMatrix &vertexData,
		double h) ;
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
} // end namespace Creation
#endif


