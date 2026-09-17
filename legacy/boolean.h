//
// Filename     : boolean.h
// Description  : Classes describing some boolean reactions applied as updates
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//
#ifndef BOOLEAN_H
#define BOOLEAN_H

#include <cmath>

#include "baseReaction.h"
#include "myRandom.h"
#include "tissue.h"

/// @brief Reactions implementing boolean rules via the update function
///
/// @details This collection of reactions updates cell variables according to boolean
/// rules. It uses the update function and directly sets the values in the output variable
/// instead of implementing an update of the derivative of the variables. A set of the
/// reactions have 'Count' in the name, which indicates that instead of updating the output
/// variable (to 0/1), the reactions add 1 everytime the boolean expression is fulfilled.
///
namespace Boolean {

/// @brief This logical gate function makes a downstream species reversibly or
/// irreversibly switch from 0 to 1 if the two input variables are 1.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndGate 1 2 2 1  # number of parameters is one
/// gate_type	     # the gate_type parameter takes the values 0 and 1 for
///                  # defining the reversible and irreversible gate, respectively.
/// index_var1       # index of the fist variable upstream the gate.
/// index_var2       # index of the second variable upstream the gate.
/// index_var_out    # updated index where the output of the gate is written. @endverbatim
/// 
class AndGate : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    /// @param indValue vector of vectors with variable indices
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    AndGate(std::vector<double> &paraValue,
            std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species reversibly or
/// irreversibly switch from 0 to 1 if a first input variable is 1 and a
/// second input variable is 0.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndNotGate 1 2 2 1 # number of parameters is set to one
/// gate_type          # the gate_type parameter takes the values 0 and 1 for
///                    # defining the reversible and irreversible gate, respectively.
/// index_var1         # index of the fist variable upstream the gate.
/// index_var2         # index of the second variable upstream the gate.
/// index_var_out      # updated index where the output of the gate is written. @endverbatim
/// 
class AndNotGate : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AndNotGate(std::vector<double> &paraValue,
               std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species switch from 0 to 1 if the
/// following conditions are met: the first input variable is 1 the second input variable is 0
/// the third input variable is 0.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndSpecialGate 0 2 3 1 # number of parameters is set to zero
/// index_var1   	   # index of the fist variable upstream the gate.
/// index_var2   	   # index of the second variable upstream the gate.
/// index_var3   	   # index of the third variable upstream the gate.
/// index_var_out  	   # updated index where the output of the gate is written. @endverbatim
/// 
class AndSpecialGate : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AndSpecialGate(std::vector<double> &paraValue,
                   std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species switch from 0 to 1 if the
/// following conditions are met: the first input variable is 1, the second input variable is 1,
/// the third input variable is 0.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndSpecialGate2 0 2 3 1 # number of parameters is set to zero
/// index_var1              # index of the fist variable upstream the gate.
/// index_var2   	    # index of the second variable upstream the gate.
/// index_var3   	    # index of the third variable upstream the gate.
/// index_var_out  	    # updated index where the output of the gate is written. @endverbatim
/// 
class AndSpecialGate2 : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AndSpecialGate2(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @note This logical gate function makes a downstream species
/// switch from 0 to 1 if the folowing conditions are met:
/// the first input variable is higher than a threshold
/// the second input variable is 1
/// the third input variable is 0.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndSpecialGate3 1 2 3 1 # number of parameters is set to zero
/// thresh      	    # threshold variable
/// index_var1   	    # index of the fist variable upstream the gate.
/// index_var2   	    # index of the second variable upstream the gate.
/// index_var3   	    # index of the third variable upstream the gate.
/// index_var_out  	    # updated index where the output of the gate is written. @endverbatim
/// 
class AndSpecialGate3 : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AndSpecialGate3(std::vector<double> &paraValue,
                    std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///

    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species add +1
/// if the two input variables are 1.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndGateCount 0 2 2 1 # number of parameters is set to zero
/// index_var1   	 # index of the first variable upstream the gate.
/// index_var2   	 # index of the second variable upstream the gate.
/// index_var_out  	 # updated index where the output of the gate is written. @endverbatim
/// 
class AndGateCount : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AndGateCount(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species add +1
/// if one of the two input variables is 1.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// OrGateCount 0 2 2 1 # number of parameters is set to zero
/// index_var1   	 # index of the fist variable upstream the gate.
/// index_var2   	 # index of the second variable upstream the gate.
/// index_var_out  	 # updated index where the output of the gate is written. @endverbatim
/// 
class OrGateCount : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    OrGateCount(std::vector<double> &paraValue,
                std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species add +1
/// if the first input variables is 1 or if the second input variable is
/// larger than 0.
///
/// @brief In the model file, the reaction is specified as:
/// @verbatim
/// OrSpecialGateCount 0 2 2 1 # number of parameters is set to zero
/// index_var1   	       # index of the fist variable upstream the gate.
/// index_var2   	       # index of the second variable upstream the gate.
/// index_var_out  	       # index where the output of the gate is written. @endverbatim
/// 
class OrSpecialGateCount : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    OrSpecialGateCount(std::vector<double> &paraValue,
                       std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species irreversibly
/// switch from 0 to 1 if the two input variables are larger than their respective thresholds.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// AndThresholdsGate 2 2 2 1 # 2 parameters, 2 index types, 2 inputs 1 output
/// thresh_var1   	      # threshold of the fist variable upstream the gate.
/// thresh_var2   	      # threshold of the second variable upstream the gate.
/// index_var1   	      # index of the first variable upstream the gate.
/// index_var2   	      # index of the second variable upstream the gate.
/// index_var_out  	      # index where the output of the gate is written. @endverbatim
/// 
class AndThresholdsGate : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    AndThresholdsGate(std::vector<double> &paraValue,
                      std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This function makes a downstream species add +1.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// Count 0 1 1   # number of parameters is set to zero
/// index_var_out # index where the output of the gate is written. @endverbatim
///
class Count : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    Count(std::vector<double> &paraValue,
          std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

/// @brief This logical gate function makes a downstream species add +1
/// if the input variable is 1.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// FlagCount 0 2 1 1 # number of parameters is set to zero
/// index_var         # index of the variable upstream the gate.
/// index_var_out     # index where the output of the gate is written. @endverbatim
///
class FlagCount : public BaseReaction {
   public:
    ///
    /// @brief Main constructor
    ///
    /// This is the main constructor which checks and sets the parameters and
    /// variable indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    ///
    /// @param indValue vector of vectors with variable indices
    ///
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    FlagCount(std::vector<double> &paraValue,
              std::vector<std::vector<size_t> > &indValue);

    ///
    /// @brief This class does not use derivatives for updates.
    ///
    /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
    ///
    void derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, DataMatrix &cellDerivs,
                DataMatrix &wallDerivs, DataMatrix &vertexDerivs);

    void derivsWithAbs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs,
                       DataMatrix &sdydtCell, DataMatrix &sdydtWall,
                       DataMatrix &sdydtVertex);
    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(double h, double t, ...)
    ///
    void update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                DataMatrix &vertexData, double h);
};

}  // end namespace Boolean

#endif  // BOOLEAN_H
