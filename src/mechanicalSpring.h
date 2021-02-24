//
// Filename     : mechanicalSpring.h
// Description  : Classes describing mechanical updates for wall springs
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2007
// Revision     : $Id:$
//
#ifndef MECHANICALSPRING_H
#define MECHANICALSPRING_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief Namespace containing mechanical 'reactions' connected to walls, i.e. 1D edges between 2D cells
///
/// Classes describing updates related to 'spring'-like updates, i.e. the vertex positions connected by an edge are
/// updated depending on their position, the resting length of the edge and possible other inputs such as molecular
/// concentration in a cell.
///
namespace WallMechanics {
  ///
  /// @brief Updates vertices from an asymmetric wall spring potential
  ///
  /// The update (in all dimensions) are given by
  ///
  /// @f[ \frac{dx_i}{dt} = (x_{i}-x_{j}) \frac{K_{force}}{L_{ij}}(1-\frac{L_{ij}}{d}) @f]
  ///
  /// where @f$ d @f$ is the distance between vertices,
  /// @f$ x_i,x_j @f$ is the vertex position,
  /// @f$ L_{ij} @f$ is the variable for storing the resting length of the wall.
  ///
  /// The parameters are @f$ K_{force} @f$ (parameter(0)), which sets the strength
  /// of the spring (spring constant), and @f$ K_{adh} @f$ (parameter(1)), which
  /// sets the relative strength of attractive forces compared to repressive
  /// forces (when attractive forces, the two parameters are multiplied 
  /// (@f$ K=K_{force}K_{adhFrac} @f$), and hence @f$K_{adhFrac}@f$ would usually
  /// be set to one or a lower value to make wall 'shrinkage' less probable. 
  /// The update needs the index of the wall length variable at the 
  /// first level (variableIndex(0,0)), and 
  /// optionally a wall variable index for storing the total wall Force (variableIndex(1,0)). 
  ///
  /// In a model file the reaction is defined as:
  ///
  /// @verbatim
  /// WallMechanics::Spring 2 1 1
  /// K_force K_adh
  /// L_ij-index
  /// @endverbatim
  /// or, when the force is saved in wall variable:
  /// @verbatim
  /// WallMechanics::Spring 2 2 1 1
  /// K_force K_adh
  /// L_ij-index
  /// Forcesave-index
  /// @endverbatim
  /// A third alternative is available for setting a different spring constant ( @f$ K_{force2} @f$
  /// , parameter(2)) for walls where a wall variable is set exactly to 1, and then the third index-layer
  /// holds the index of the 'flag' variable.
  /// @verbatim
  /// WallMechanics::Spring 3 3 1 1/0 1
  /// K_force K_adh K_force2
  /// L_ij-index
  /// [Forcesave-index]
  /// wall_type_index
  /// @endverbatim
  /// @note This reaction used to be called VertexFromWallSpring
  ///
  class Spring : public BaseReaction {

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
      Spring(std::vector<double> &paraValue, 
          std::vector< std::vector<size_t> > 
          &indValue);

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
  /// @brief Updates vertices of epidermal walls from an asymmetric wall spring potential
  ///
  /// Mechanical spring with possibility to have stronger repressive spring constants compared to
  /// attractive.
  /// Same as WallMechanics::Spring (basic version), except that only walls connected to the tissue
  /// boundary/background are updated.
  ///
  /// @see WallMechanics::Spring
  ///
  class SpringEpidermal : public BaseReaction {

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
      SpringEpidermal(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief Updates vertices of epidermal cells from an asymmetric wall spring potential
  ///
  /// Mechanical spring with possibility to have stronger repressive spring constants compared to
  /// attractive.
  /// Same as WallMechanics::Spring (basic version), except that only walls connected to epidermal
  /// cells are updated. Epidermal cells are defined to be neighbor with the tissue boundary/background.
  ///
  /// @see WallMechanics::Spring
  ///
  class SpringEpidermalCell : public BaseReaction {

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
      SpringEpidermalCell(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief Updates vertices from an asymmetric wall spring potential depending
  /// on cell concentrations via a (decreasing) Hill function
  ///
  /// Same as basic WallMechanics::Spring, but where the spring constant depends on a cell
  /// variable in the form (contributions from cells on either side)
  ///
  /// @f[ K_{spring} = p_{0} + p_{1} (\frac{p_{2}^{p_{3}}}{p_{2}^{p_{3}} + c_{1}^{p_{3}}} +
  /// \frac{p_{2}^{p_{3}}}{p_{2}^{p_{3}}+c_{2}^{p_{3}}} ) @f]
  ///
  /// where @f$ p_{0},p_{1} @f$ sets the range of spring constant values, and @f$p_{2},p_{3}@f$
  /// are the Hill constant and coefficient, respectively.
  ///
  /// This function is implemented as stiffness decreasing with the concentration as used for
  /// e.g. auxin contribution. To be used for increasing stiffness ith concentration, replace
  /// @f$p_{2} -> -p_{2}, p_{1} -> p_{1}+2p_{2}@f$ (verify before you use).
  ///
  /// As for WallMechanics::Spring, an additional parameter (@f$p_{4}@f$) can set a ratio of
  /// attractive vs repressive force.
  ///
  /// @see WallMechanics::Spring for spring force calculation.
  ///
  class SpringConcentrationHill : public BaseReaction {

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
      SpringConcentrationHill(std::vector<double> &paraValue, 
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
  };

  ///
  /// @brief Updates vertices from a wall spring potential depending on 1)
  /// the cell concentration, and 2) whether the cell is epidermal or not.
  ///
  /// Similar to the basic reactions such as WallMechanics::Spring, but takes
  /// into account the cell concentration (additive contributions from both 
  /// sides) and whether the cell is epidermal or not.
  /// 
  /// @f[
  /// K_{spring}=\begin{cases}
  /// p_{0}&c<threshold~\text{and~cell~epidermal} \\
  /// p_{1}&c>threshold~\text{and~cell~epidermal} \\
  /// p_{2}&c<threshold~\text{and~cell~not~epidermal} \\
  /// p_{3}&c>threshold~\text{and~cell~not~epidermal}
  /// \end{cases}
  /// @f]
  ///
  /// Here @f$p_i@f$ denotes the spring constant added for every cell
  /// concentration.
  /// As for WallMechanics::Spring, a sixth parameter (fraction/adhesion) can set a ratio of
  /// attractive vs repressive force.
  /// @verbatim
  ///   SpringInternalExternalThreshold 6 1 2  
  ///   threshold p0 p1 p2 p3
  ///   frac_adh
  ///   wall_length_index
  ///   species_index
  /// @endverbatim
  /// @note This process is done additively for both cells connected to a cell
  /// wall.
  ///
  /// @see WallMechanics::Spring for spring force calculation.
  ///
  class SpringInternalExternalThreshold : public BaseReaction {

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
      SpringInternalExternalThreshold(std::vector<double> &paraValue, 
          std::vector< std::vector<size_t> > 
          &indValue);
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

  ///
  /// @brief A viscoelastic update of edges
  ///
  /// The update (in all dimensions) are given by
  ///
  /// @f[ \frac{dx_i}{dt} = @f] 
  ///
  /// In a model file the reaction is defined as:
  ///
  /// @verbatim
  /// WallMechanics::ViscoElastic 2 1 1
  /// eta        # viscosity
  /// k          # spring constant
  /// L_ij-index # index of the resting lenth
  /// @endverbatim
  /// or
  /// @verbatim
  /// WallMechanics::ViscoElastic 2 2 1 1
  /// eta          # viscosity
  /// k            # spring constant
  /// L_ij-index   # index of the resting lenth
  /// strain_index # optional wall index to store strain
  /// @endverbatim
  
  ///
  class ViscoElastic : public BaseReaction {

  private:
    DataMatrix historyData_;
    std::vector<double> historyTime_;
      
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
    ViscoElastic(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > 
		 &indValue);
    
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
    
    void initiate(Tissue &T,
		  DataMatrix &cellData,
		  DataMatrix &walldata,
		  DataMatrix &vertexData,
		  DataMatrix &cellderivs,
		  DataMatrix &wallderivs,
		  DataMatrix &vertexDerivs );
    void update(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData, 
		double h);       
  };
  
} // end namespace WallMechanics

///
/// @brief Updates vertices from an asymmetric wall spring potential
///
/// 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromWallSpringMTnew 3/4 2 1 1 
/// K_matrix K_adh K_fiber (double_length_flag)
/// L_ij-index
/// MT_index
///
/// or
///
/// VertexFromWallSpringMTnew 3/4 3 1 1 1
/// K_matrix K_adh K_fiber (double_length_flag)
/// L_ij-index
/// MT_index
/// forceSave_index
///
/// @endverbatim
///
/// 
class VertexFromWallSpringMTnew : public BaseReaction {

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
    VertexFromWallSpringMTnew(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > 
        &indValue );

    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );

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
};



///
/// @brief Updates vertices from an asymmetric wall spring potential
///
/// Similar to VertexFromWallSpring but only acts on the boundary walls,
/// where boundary is defined as the wall being connected to the
/// outside/background of the tissue.
///
///
class VertexFromWallBoundarySpring : public BaseReaction {

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
    VertexFromWallBoundarySpring(std::vector<double> &paraValue, 
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
};

namespace CenterTriangulation {
  ///
  /// @brief Updates vertices from an asymmetric spring potential on internal edges
  ///
  /// The update (in all dimensions) are given by
  ///
  /// @f[ \frac{dx_i}{dt} = (x_{i}-x_{j}) \frac{K_{force}}{L_{ij}}(1-\frac{L_{ij}}{d}) @f]
  ///
  /// where @f$ d @f$ = distance between vertices (vertex and center point),
  /// where @f$ x_i,x_j @f$ = vertex position in specific dimension,
  /// @f$ L_{ij} @f$ = resting length of internal edge
  ///
  /// The parameters are @f$ K_{force} @f$ (parameter(0)), which sets the strength
  /// of the spring (spring constant), and @f$ K_{adh} @f$ (parameter(1)), which
  /// sets the relative strength of adhesive forces compared to repressive
  /// forces (when adhesive forces, the two parameters are multiplied 
  /// (@f$ K=K_{force}K_{adhFrac} @f$). 
  /// The column index for the cell additional variables of the central mesh 
  /// (x,y,z,L_1,...,L_n) should be given in the first level of indices.
  ///
  /// In a model file the reaction is defined as
  ///
  /// @verbatim
  /// CenterTriangulation::EdgeSpring 2 1 1
  /// K_force K_adh
  /// index
  /// @endverbatim
  ///
  /// @see VertexFromWallSpring
  ///
  class EdgeSpring : public BaseReaction {

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
      EdgeSpring(std::vector<double> &paraValue, 
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
  };
}

///
/// @brief Updates vertices from an asymmetric wall spring potential
/// where the spring constants are taken from two wall variables
///
/// This class updates according to a spring potential where two
/// spring constants are set from elsewhere, i.e. values are taken
/// from two variables representing the values for the TWO parts of
/// the wall connected to the two neighboring cells. 
///
class VertexFromDoubleWallSpring : public BaseReaction {

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
    VertexFromDoubleWallSpring(std::vector<double> &paraValue, 
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
};

///
/// @brief Updates vertices from an asymmetric wall spring potential with a spatial factor
///
class VertexFromWallSpringSpatial : public BaseReaction {

  private:

    double Kpow_;

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
    VertexFromWallSpringSpatial(std::vector<double> &paraValue, 
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
};

///
/// @brief Updates vertices from an asymmetric wall spring potential with a spatial MT factor
///
/// This function has a base level for the spring constant and then adds a
/// factor dependent on the MT direction in the cells as well as a spatial
/// factor. The spring constant of the wall is determined by:
///
/// @f[ k = p_0 + (p_1 + p_2*spatialFactor)*mtFactor @f]
///
/// where the spatial factor is determined by a Hill function and the MT
/// factor comes from the absolute value of the scalar products of the wall
/// vector and cell direction vectors.
///
class VertexFromWallSpringMTSpatial : public BaseReaction {

  private:

    double Kpow_;

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
    VertexFromWallSpringMTSpatial(std::vector<double> &paraValue, 
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
};

///
/// @brief Updates vertices from an spatially asymmetric wall spring potential
/// given by microtubule directions
///
/// This function has a base level for the spring constant and then adds a
/// factor dependent on the MT direction in the cells. The spring constant 
/// of the wall is determined by:
///
/// @f[ k = p_0 + p_1*|n_{wall} n_{MT}| @f]
///
/// where the MT factor comes from the absolute value of the scalar products 
/// of the wall vector and cell direction (microtubular) vector. An additional
/// parameter @f$p_{2}@f$ sets an additional multiplied contribution if the spring
/// is contracted.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromWallSpringMT 3 2 2 0/1
/// K_homogeneous K_MT f_adh
/// L_ij-index MT-index
/// [Forcesave-index]
/// @endverbatim
///
/// Alternatively if no force save index is supplied the first line
/// can be replaced by 'VertexFromWallSpringMT 3 1 2'.
///
class VertexFromWallSpringMT : public BaseReaction {

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
    VertexFromWallSpringMT(std::vector<double> &paraValue, 
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
};

///
/// @brief Updates vertices from an spatially asymmetric wall spring potential
/// given by microtubule directions and also updates the spring constants
/// slowly.
///
class VertexFromWallSpringMTHistory : public BaseReaction {

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
    VertexFromWallSpringMTHistory(std::vector<double> &paraValue, 
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
    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );
};

class VertexFromWallSpringMTConcentrationHill : public BaseReaction {

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
    VertexFromWallSpringMTConcentrationHill(std::vector<double> &paraValue, 
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
};

///
/// @brief Sets mechanical properties for wall segments from MT
/// directions and (auxin) concentrations
///
/// Uses a square scalar product between wall and MT directions and a
/// Hill formalism for auxin dependent wall mechanics. The walls are
/// divided into two segments 'belonging to' respective neighboring
/// cell.
///
class VertexFromDoubleWallSpringMTConcentrationHill : public BaseReaction {

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
    VertexFromDoubleWallSpringMTConcentrationHill(std::vector<double> &paraValue, 
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
};




///
/// @brief Sets external  breakable springs between vertices at given indices
///
/// This function has a spring constant and adhision factor 
///
/// @f[ F=-K*f_{adh}*(L_0-X) @f]
///
/// when the distance between vertices is greater than Lmax Force becomes zero.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromExternalSpring 4 3 1 1/.../n 1/.../n
///
/// K f_adh Lmaxfactor growth_rate 
/// growth_flag
/// vertex-index11 ... vertex-index1n 
/// vertex-index21 ... vertex-index2n  
/// 
/// @endverbatim
///
///
class VertexFromExternalSpring : public BaseReaction {

  private: 
    std:: vector<double> restinglength;
    std:: vector<double> Kspring;
    size_t Npairs;

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
    VertexFromExternalSpring(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs);

    void derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        double h);


};



///
/// @brief Sets external  breakable springs between vertices which see eachother 
/// in a given angle intervall respect to the normal direction to their membrane.
///
/// This function has a spring constant and adhision factor 
///
/// @f[ F=-K*f_{adh}*(L_0-X) @f]
///
/// when the distance between vertices is greater than Lmax Force becomes zero.
///
/// In a model file the reaction is defined as
///
/// @verbatim
///
/// VertexFromExternalSpringFromPerpVertex 8 1 4 
/// K   
/// f_adh 
/// Lmaxfactor 
/// growth_rate  
/// intraction_angle 
/// corner_angle
/// growth_rate_decay_rate
/// growthStress
///
/// growth_flag (0:non ,1: ,2: ,3: ,4: ,5: ,6: 7: 8: )
/// connection_flag (1: constraint on the first node only, 2: for both)
/// exclude_corner_flag (0:include corners, 1: exclude corners)
/// initiate_flag       (0: for all vertices, 
///                      1: only the vertices in the sisterVertex list)
/// @endverbatim
///
///
class VertexFromExternalSpringFromPerpVertex : public BaseReaction {

  private: 

    //size_t Npairs;
    std:: vector<std::vector<std::vector<double> > >  connections;
    std:: vector<std::vector<double> >  vertexVec;

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
    VertexFromExternalSpringFromPerpVertex(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs);

    void derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        double h);
    ///
    /// @brief Prints internal variables for plotting
    ///
    /// Plots the internal edges for plotting using gnuplot. The format is 
    /// t_i i x1 y1 [z1] edgeLength
    /// t_i i x2 y2 [z2] edgeLength
    /// for each internal edge.
    ///
    /// @note Requires the BaseSolver::print() function to call printState()
    /// @see BaseReaction::printState()
    ///
    void printState(Tissue *T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        std::ostream &os=std::cout);
};


////////////////////////////////////////////////////////////

///
/// @brief Sets external  breakable springs between vertices which see eachother 
/// in a given angle intervall respect to the normal direction to their membrane.
///
/// This function has a spring constant and adhision factor 
///
/// @f[ F=-K*f_{adh}*(L_0-X) @f]
///
/// when the distance between vertices is greater than Lmax Force becomes zero.
///
/// In a model file the reaction is defined as
///
/// @verbatim
///
/// VertexFromExternalSpringFromPerpVertexDynamic 8 1 4 
/// K   
/// f_adh 
/// Lmaxfactor 
/// growth_rate  
/// intraction_angle 
/// corner_angle
/// growth_rate_decay_rate
/// growthStress
///
/// growth_flag (0:non ,1: ,2: ,3: ,4: ,5: ,6: 7: 8: )
/// connection_flag (1: constraint on the first node only, 2: for both)
/// exclude_corner_flag (0:include corners, 1: exclude corners)
/// initiate_flag       (0: for all vertices, 
///                      1: only the vertices in the sisterVertex list)
/// @endverbatim
///
///
class VertexFromExternalSpringFromPerpVertexDynamic : public BaseReaction {

  private: 

    //size_t Npairs;
    std:: vector<std::vector<std::vector<double> > >  connections;
    std:: vector<std::vector<double> >  vertexVec;

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
    VertexFromExternalSpringFromPerpVertexDynamic(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs);

    void derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        double h);
    ///
    /// @brief Prints internal variables for plotting
    ///
    /// Plots the internal edges for plotting using gnuplot. The format is 
    /// t_i i x1 y1 [z1] edgeLength
    /// t_i i x2 y2 [z2] edgeLength
    /// for each internal edge.
    ///
    /// @note Requires the BaseSolver::print() function to call printState()
    /// @see BaseReaction::printState()
    ///
    void printState(Tissue *T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        std::ostream &os=std::cout);
};



///
/// @brief A repulstion spring force with a given constant
/// between nodes on the vertices of two cells
/// in the places that membranes intersect.
///
/// @details In a model file the reaction is defined as
/// @verbatim
/// cellcellRepulsion 1 0 
/// K   
/// checking_radius
/// @endverbatim
///

class cellcellRepulsion : public BaseReaction {

  private: 

    std:: vector<std::vector<double> >  vertexVec;
    std:: vector<std::vector<int> > grid;
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
    cellcellRepulsion(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs);

    void derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        double h);
};

///
/// @brief Assumes the template on a stretching substrate 
/// connected to the vertices with springs with a given constant
/// @details In a model file the reaction is defined as
/// @verbatim
/// vertexFromSubstrate 7 0 
/// K_sp
/// stretch_rate
/// x_direction
/// y_direction
/// x_center
/// y_center
/// vertex_percentage
/// @endverbatim
///

class vertexFromSubstrate : public BaseReaction {

  private: 
    std::vector<size_t> list;
    std::vector<std::vector<std::vector<double> > >  vertexVec;
    //std:: vector<std::vector<double>>  vertexVec;
    size_t numAttachedCells;
    //size_t numAttachedVertices;
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
    vertexFromSubstrate(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
    void initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs);

    void derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs );

    ///
    /// @brief Update function for this reaction class
    ///
    /// @see BaseReaction::update(Tissue &T,...)
    ///
    void update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData, 
        double h);
};




#endif
