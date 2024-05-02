//
// Filename     : growth.h
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef WALLGROWTH_H
#define WALLGROWTH_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

/// 
/// @brief Reactions describing wall growth
///
/// @details Different reactions for wall growth implemented as increase of the resting length of the walls.
/// Either a constant or strain/stress based approach is used, and in some cases the possibility to
/// truncate the growth at specific lengths can be applied. Also versions reading cell variables
/// (typically a concentration) that affects the growth rate are implemented and named ConcentrationHill.
/// In all cases the index of the resting length is given and this is typically index 0.
/// In some cases named stress, strain or stress (read from wall variables) can be used and selected by a flag.
/// In some cases one can choose between const (no) or linear dependence on the current length of the wall,
/// again given by a flag.
///
namespace WallGrowth {

  ///
  /// @brief Constant/exponential wall growth which can be truncated at threshold length
  ///
  /// @details Constant or exponential (proportional to length) growth that can be truncated at 
  /// a threshold length. 
  /// The wall lengths, L, are updated according to
  /// @f[ \frac{dL}{dt} = p_0 L (1-\frac{L}{p_2}) @f]
  /// if linearFlag (@f$p_1@f$) is set, or 
  /// @f[ \frac{dL}{dt} = p_0 (1-\frac{L}{p_2}) @f]
  /// otherwise. @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$), @f$ p_1 @f$ is the flag setting
  /// the update version (0 for constant and 1 for linear dependence on the current length), and 
  /// the optional parameter @f$ p_2 @f$ is the threshold for growth truncation (@f$ L_{trunc} @f$).
  /// If this parameter is not given, no truncation is applied.
  /// The column index for the wall length should be given at the first level.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// WallGrowth::Constant 2/3 1 1
  /// k_growth linear_flag [L_trunc]
  /// L @endverbatim
  ///
  class Constant : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// @details This is the main constructor which sets the parameters and variable
  /// indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  /// @param indValue vector of vectors with variable indices
  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  ///  
  Constant(std::vector<double> &paraValue, 
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
  /// @brief Stress/strain-driven wall growth dependent on a threshold.
  ///
  /// @details Constant (constant mode) or exponential (linear mode) growth driven by a 
  /// streched/stressed wall. The wall 
  /// lengths are updated only if the length plus a threshold value is shorter 
  /// than the distance between the vertices of the wall (in strain mode), and if 
  /// the total stress is above the threshold (in stress mode). Hence the update follows:
  ///  @f[ \frac{dL}{dt} = p_0 L (S-p_1) @f]
  /// if @f$ S > p_1 @f$ and 0 otherwise.
  /// where S is the stress/strain and L is the wall length.
  /// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$). 
  /// @f$ p_1 @f$ is a threshold (@f$ s_{threshold} @f$) in stress or strain depending 
  /// on @f$ p_2 @f$. If set to zero, shinkage is allowed. 
  /// @f$ p_2 @f$ is a flag (@f$ strain_{flag} @f$) for using stretch/strain (@f$p_2=1 @f$) 
  /// or stress (@f$p_2=0 @f$). Strain is calculated by @f$ S=(d-L) / L @f$, where d is the distance
  /// between the vertices. Stress is read from the second layer of variable indices.  
  /// @f$ p_3 @f$ is a flag (@f$ linear_{flag} @f$) for using growth proportional to 
  /// wall length (@f$ p_3=1 @f$, as in equation above. If @f$ p_3=0 @f$ the above equation
  /// will not include the length L on the right-hand side.
  /// If a 5th parameter is given, the growth is truncated at a maximal length by multiplying the 
  /// rate with a factor @f$ (1 - L/p_4) @f$.
  /// The column index for the wall length should be given in the first layer, and if stress
  /// is to be used, a second layer with calculated wall stress variables has to be given.
  /// In a model file the reaction is defined as
  /// @verbatim
  /// WallGrowthStress 4/5 1/2 1 [N]
  /// k_growth s_threshold stretch_flag linear_flag [L_trunc]  
  /// L
  /// [stress1 ... stressN] @endverbatim
  /// If stress is used (stretch_flag=0) a second level of wall stresses has to be read
  /// (calculated and updated from other (mechanical) reactions).
  ///
  /// @note If s_threshold is set to zero, also shrinkage is allowed. To avoid shrinkage set small value.
  ///
  class Stress : public BaseReaction {
    
  public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    /// @param indValue vector of vectors with variable indices
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///  
    Stress(std::vector<double> &paraValue, 
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
  /// @brief almansi-strain-driven wall growth dependent on a threshold via "update".
  ///
  /// @details In a model file the reaction is defined as
  /// @verbatim
  /// WallGrowthStrain 4 1 1 
  /// k_growth s_threshold strain_flag linear_flag   
  /// L0-index @endverbatim
  ///
  /// If Almansi strain is used (strain_flag=0) .
  ///
  /// @note If s_threshold is set to zero, also shrinkage is allowed. To avoid shrinkage set small value.
  ///
  class Strain : public BaseReaction {
    
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
    Strain(std::vector<double> &paraValue, 
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

    void update(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                double h );
  };
      
  ///
  /// @brief Constant stress/strain-driven wall growth dependent on a
  /// distance to an maximal vertex coordinate (e.g. tip) in specified dimension.
  ///
  /// @details Constant growth driven by a streched wall. The wall lengths, L, are
  /// updated only if the length is shorter than the distance between the
  /// vertices of the wall and then according to
  /// @f[\frac{dL}{dt} = p_{0}(d_{v}-L-p_{1}) \frac{p_{2}^{p_{3}}}{(p_{2}^{p_{3}}+d^{p_{3}})}@f]
  /// iff @f$(d_{v}-L) > p_{1}@f$
  /// @f$p_0@f$ is the growth rate.
  /// @f$p_1@f$ is a stress/strain threshold
  /// @f$p_2@f$ is the K_Hill of the spatial factor
  /// @f$p_3@f$ is the n_Hill of the spatial factor
  /// @f$p_4@f$ is a flag for using stretch/strain instead of stress
  /// @f$p_5@f$ is a flag for using growth proportional to wall length (not constant)
  /// @f$d_v@f$ is the distance between the two wall vertices.
  /// @f$d@f$ is the total (all dimensions) distance between the max value and wall.
  /// In addition, the column index for the wall length, the distance
  /// coordinate should be given at first level and stress index in second.
  /// @note The difference to StressSpatialSingle is that the distance is the total (all dimensions) distance.
  class StressSpatial : public BaseReaction {
    
  private:
    
    double Kpow_;
    
  public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    /// @param indValue vector of vectors with variable indices
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    StressSpatial(std::vector<double> &paraValue,
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
  /// @brief Constant stress/strain-driven wall growth dependent on a
  /// distance to an maximal vertex coordinate (e.g. tip) in specified dimension
  ///
  /// @details Constant growth driven by a streched wall. The wall lengths, L, are
  /// updated only if the length is shorter than the distance between the
  /// vertices of the wall and then according to
  /// @f[ \frac{dL}{dt} = p_{0} (d_{v}-L-p_{1}) \frac{p_{2}^{p_3}}{(p_{2}^{p_{3}}+d^{p_{3}})} @f]
  /// iff @f$(d_{v}-L) > p_{1}@f$.
  /// @f$p_0@f$ is the growth rate.
  /// @f$p_1@f$ is a stress/strain threshold
  /// @f$p_2@f$ is the K_Hill of the spatial factor
  /// @f$p_3@f$ is the n_Hill of the spatial factor
  /// @f$p_4@f$ is a flag for using stretch/strain instead of stress
  /// @f$p_5@f$ is a flag for using growth proportional to wall length (not constant)
  /// @f$d_v@f$ is the distance between the two wall vertices in the specified dimension.
  /// @f$d@f$ is the distance between the max value and wall.
  /// In addition, the column index for the wall length, the distance
  /// coordinate should be given at first level and stress index in second.
  /// @note The difference to Stress Spatial is that the distance is measured in the specified dimension only.
  class StressSpatialSingle : public BaseReaction {
    
  private:
    
    double Kpow_;
    
  public:
    ///
    /// @brief Main constructor
    ///
    /// @details This is the main constructor which sets the parameters and variable
    /// indices that defines the reaction.
    ///
    /// @param paraValue vector with parameters
    /// @param indValue vector of vectors with variable indices
    /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
    ///
    StressSpatialSingle(std::vector<double> &paraValue,
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
  /// @brief Constant stress/strain-driven wall growth dependent on a
  /// concentration level in the cell
  ///
  /// @details Constant growth driven by a streched wall. The wall lengths, L, are
  /// updated only if the length is shorter than the distance between the
  ///  vertices of the wall and then according to
  /// @f[\frac{dL}{dt} = (p_{0} + p_{1} f(c,p_2,p_3)) (d_{v}-L)@f]
  /// iff @f$(d_{v}-L) > p_{4}@f$
  /// @f$p_{0}@f$ is the constant growth rate,
  /// @f$p_{1}@f$ is the maximal added growth rate (V_max in the Hill function),
  /// @f$p_{2}@f$ is the K_Hill,
  /// @f$p_{3}@f$ is the n_Hill
  /// @f$p_{4}@f$ is a threshold
  /// @f$p_{5}@f$ is a flag for using stretch/strain instead of stress
  /// @f$p_{6}@f$ is a flag for using growth proportional to wall length (not constant)
  /// @f$ f = \frac{c^{n}}{c^{n}+K^{n}}@f$ is an increasing hill function
  /// @f$c@f$ is the concentration
  /// @f$d_{v}@f$ is the distance between the two wall vertices.
  /// In addition, the column index for the wall length and the cell
  /// concentration should be given and wall stress measures if used.
  /// In a model fle the reaction should be defined as
  /// @verbatim
  /// WallGrowth::StressConcentrationHill 7 2 2 n
  /// k_growth k_Growth_hill K_hill n_hill threshold strain_flag lin_flag
  /// L(0) c
  /// [S1 ...] @endverbatim
  ///
  /// where L is the wall index for length (should usually be 0), c is the cell
  /// variable index for the concentration, and S are the wall stresses calculated
  /// elsewhere (only required if stress is used as growth signal (p_5=0).
  ///
  class StressConcentrationHill : public BaseReaction {
    
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
    StressConcentrationHill(std::vector<double> &paraValue,
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
  /// @brief Constant strech-driven wall growth with epidermal walls treated specially
  ///
  /// @details Constant growth driven by a streched wall. The wall lengths, L, are
  /// updated only if the length is shorter than the distance between the
  /// vertices of the wall and then according to
  ///
  /// @f[ \frac{dL}{dt} = p_0 (d_v-L) f_e @f]
  ///
  /// @f$p_0@f$ is the growth rate.
  /// @f$f_e@f$ fraction for epidermal walls
  /// @f$d_v@f$ is the distance between the two wall vertices.
  ///
  /// In addition, the column index for the wall length should be given.
  ///
  class ConstantStressEpidermalAsymmetric : public BaseReaction {
    
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
    ConstantStressEpidermalAsymmetric(std::vector<double> &paraValue,
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
  /// @brief Reads Force(s) from wall variables and increase resting length if larger than
  /// threshold
  ///
  /// @details Implements growth as
  /// @f[ \frac{dL}{dt} = p_0 (F-p_1) @f]
  /// if @f$ F > p_1 @f$ and 0 otherwise.
  /// where F is added up from wall variables given as indices. In a model file the reaction
  /// is given as:
  /// @verbatim
  /// WallGrowth::Force 2 2 1 N
  /// k_growth F_threshold
  /// L
  /// F_1 ... F_N @endverbatim
  ///
  /// This is an old reaction that is replaced by the more general WallGrowth::Stress(), but
  /// hangs around since it was used in publications.
  ///
  /// @note used to be called WallLengthGrowExperimental.
  /// @see WallGrowth::Stress()
  ///
  class Force : public BaseReaction {
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
    Force(std::vector<double> &paraValue,
	  std::vector< std::vector<size_t> >
	  &indValue );
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
  
  namespace CenterTriangulation {
    
    ///
    /// @brief Constant internal edge growth which can be truncated at threshold length
    ///
    /// @details Constant (constant mode) or exponential (proportional to length, linear mode) 
    /// growth of the internal edges 
    /// in a central cell vertex meshed description. The internal edges 
    /// lengths are updated only if the length plus a threshold value is shorter 
    /// than the distance between the vertices (vertex and central cell vertex). 
    ///
    /// The internal edge lengths, L, are updated according to
    ///
    /// @f[ \frac{dL}{dt} = p_0 L (1-\frac{L}{p_2}) @f]
    ///
    /// if linearFlag (@f$p_1@f$) is set, or 
    ///
    /// @f[ \frac{dL}{dt} = p_0 (1-\frac{L}{p_2}) @f]
    ///
    /// otherwise. @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$), @f$ p_1 @f$ is the flag setting
    /// the update version (0 for constant and 1 for linear dependence on the current length), and 
    /// the optional parameter @f$ p_2 @f$ is the threshold for growth truncation (@f$ L_{trunc} @f$).
    /// If this parameter is not given, no truncation is applied.
    ///
    /// The column index for the cell additional variables of the central mesh 
    /// (x,y,z,L_1,...,L_n) should be given in the first level of indices.
    ///
    /// In a model file the reaction is defined as
    ///
    /// @verbatim
    /// CenterTriangulation::WallGrowth::Constant 2/3 1 1
    /// k_growth linear_flag [L_trunc]
    /// index @endverbatim
    /// @see WallGrowth::Constant (for same update of edges/2D walls)
    ///
    class Constant : public BaseReaction {
      
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
      Constant(std::vector<double> &paraValue, 
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
    /// @brief Constant stress/strain-driven internal edge growth dependent on a threshold
    ///
    /// @details Constant (constant mode) or exponential (linear mode) growth of the internal edges 
    /// in a central cell vertex meshed description and driven by a 
    /// streched/stressed wall. The internal edges 
    /// lengths are updated only if the length plus a threshold value is shorter 
    /// than the distance between the vertices (vertex and central cell vertex) 
    /// (in strain mode), and if 
    /// the total stress is above the threshold (in stress mode). Hence the update follows:
    ///
    /// @f[ \frac{dL}{dt} = p_0 L (S-p_1) @f]
    /// if @f$ S > p_1 @f$ and 0 otherwise.
    ///
    /// where S is the stress/strain and L is the wall length.
    /// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$). 
    /// @f$ p_1 @f$ is a threshold (@f$ s_{threshold} @f$) in stress or strain depending 
    /// on @f$ p_2 @f$. 
    /// @f$ p_2 @f$ is a flag (@f$ strain_{flag} @f$) for using stretch/strain (@f$p_2=1 @f$) 
    /// or stress (@f$p_2=0 @f$). Strain is calculated by @f$ S=(d-L)/L @f$, where d is the distance
    /// between the vertices. Stress is read from the second layer of variable indices.  
    /// @f$ p_3 @f$ is a flag (@f$ linear_{flag} @f$) for using growth proportional to 
    /// wall length (@f$ p_3=1 @f$, as in equation above. If @f$ p_3=0 @f$ the above equation
    /// will not include the length L on the right-hand side.
    /// If a 5th parameter is given, the growth is truncated at a maximal length by multiplying the 
    /// rate with a factor @f$ (1 - L/p_4) @f$.
    ///
    /// The column index for the cell additional variables of the central mesh (x,y,z,L_1,...,L_n) 
    /// should be given in the first layer, and if stress
    /// is to be used, a second layer with calculated wall stress variables has to be given.
    ///
    /// In a model file the reaction is defined as
    ///
    /// @verbatim
    /// CenterTriangulation::WallGrowth::Stress 4/5 1/2 1 [N]
    /// k_growth s_threshold stretch_flag linear_flag [L_trunc] 
    /// L
    /// [stress1 ... stressN] @endverbatim
    /// If stress is used (stretch_flag=0) a second level of wall stresses has to be read
    /// (calculated and updated from other (mechanical) reactions).
    ///
    /// @see WallGrowth::Stress (for same updatre of 1D walls)
    ///
    class Stress : public BaseReaction {
      
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
      Stress(std::vector<double> &paraValue, 
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
    /// @brief ADD DOCUMENTAION
    ///
    class StressConcentrationHill : public BaseReaction {
      
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
      StressConcentrationHill(std::vector<double> &paraValue, 
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
    /// @brief This reaction is currently using ad hoc additions/changes within the code and should
    /// only be used by an expert, i.e. Behruz?
    ///
    /// @details In a model file the reaction is defined as
    ///
    /// @verbatim
    /// CenterTriangulation::WallGrowth::StrainTRBS 2 3 1 1 3
    /// k_growth s_threshold  
    ///    
    /// L_ij-index
    /// InternalVarStartIndex 
    ///
    /// strain1_index
    /// strain2_index
    /// strain_vector_index @endverbatim
    ///
    /// (strain value and direction calculated and updated from other (mechanical) reactions).
    ///
    class StrainTRBS : public BaseReaction {
      
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
      StrainTRBS(std::vector<double> &paraValue, 
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
      
      void update(Tissue &T,
		  DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
                  double h );
    };

    class StrainTRBSConcentrationHill : public BaseReaction {
      
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
      StrainTRBSConcentrationHill(std::vector<double> &paraValue, 
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
      
      void update(Tissue &T,
		  DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
                  double h );
    };
    
    
    ///
    /// @brief This reaction reads a vector representing e.g. strain or stress calculated elsewhere
    /// as input for growth 
    /// 
    /// @details A vector t for example can be strain or stress rections and magnitudes are given as input to
    /// this reaction, where the information is read from cell variables hence assuming that another
    /// reaction is updating the values unless they are supposed to be constant.
    ///
    /// The reaction checks if the values are above a threshold value and then update edge elements
    /// resting lengths accordingly. The reactions use the global input and for each triangle in the
    /// cell 'project' the main directions down to the edge directions of the triangles to calculate
    /// the contribution per edge. It follows the description in
    /// @verbatim Bozorg, Krupinski and Jonsson (2016) A continuous growth model for plant tissue.
    /// Phys Biol 13:065002 @endverbatim
    /// and is an implementation of the update in Eq. 32:
    /// @f[ \frac{dL_i}{dt} = k_g R(g_i-g_t) L_i @f]
    /// where @f$ L_i @f$ is the resting length of edge i, @f$k_g@f$ growth rate, @f$ R @f$ is the ramp
    /// function (linearly increasing if the argument is above zero (zero otherwise). @f$g_t@f$ is the
    /// given threshold value and @f$g_i@f$ is the (e.g. strain) value in the direction of the edge.
    /// The update is done in the update function (not derivs), such that it can be controlled
    /// to only be done the elastic mechancs is in equilibrium. This will be done if the velocityThreshold
    /// is parameter is larger than any vertex movement from mechanical updates (as stored in a cellData
    /// variable provided). An option is to set the velocityThreshold>100, and then the relaxation requirement
    /// is not applied.
    /// In a model file the reaction is defined as
    /// @verbatim
    /// CenterTriangulation::WallGrowth::VectorTRBS 4 4 1 1 1 3
    /// k_growth              # growth rate
    /// g_threshold           # signal threshold above which there will be growth
    /// velocity_threshold    # vertex movement threshold for when growth will be applied (>100 = each time)
    /// double_edge_flag      # 0 single edge elements and 1 if double, i.e. indep edges for connected tri's
    /// L_ij-index            # edge resting length wor wall-edge element (index in wallData, usually=0)
    /// InternalVarStartIndex # center triangulation data (index in cellData)
    /// VelocityStoreIndex    # cell index for velocity data (to check for equilibrium)
    /// strain1_index         # 'signal' magnitue in first principal direction (index in cellData)
    /// strain2_index         # 'signal' magnitue in second principal direction
    /// strain_vector_index   # start index where principal direction stored (in cellData) @endverbatim
    /// @note Strain value and direction calculated and updated from other (mechanical) reactions.
    /// @see namespace (to come) TRBS, triangular spring plate elements.
    /// @see VertexFromTRBScenterTriangulation
    ///
    class VectorTRBS : public BaseReaction {
      
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
      VectorTRBS(std::vector<double> &paraValue, 
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
                  DataMatrix &vertexDerivs);

      void update(Tissue &T,
		              DataMatrix &cellData,
		              DataMatrix &wallData,
		              DataMatrix &vertexData,
                  double h );
    };
  } // namespace CenterTriangulation
} // namespace WallGrowth

#endif // WALLGROWTH_H
