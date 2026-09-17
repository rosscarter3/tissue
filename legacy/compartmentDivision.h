//
// Filename     : compartmentDivision.h
// Description  : Classes describing compartmentDivision updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : July 2006
// Revision     : $Id:$
//
#ifndef COMPARTMENTDIVISION_H
#define COMPARTMENTDIVISION_H

#include <cmath>
#include <vector>

#include "tissue.h"
#include "baseCompartmentChange.h"

/// @brief Namespace for classes describing cell division rules.
///
/// All divisions are described by a flag and an update function where the first one flags when it is time to
/// divide, e.g. when a cell size passes a threshold value, and the second describes the update, e.g. puts the
/// new wall along the shortest path dividing the cell in two equally sized daughters.
///
/// @ see BaseCompartmentChange
namespace Division {

  /// @brief Divides a cell when volume above a threshold,
  /// with new wall perpendicular to the longest wall segment.
  ///
  /// Divides a cell when volume above a threshold. New wall is created
  /// prependicular to the longest cell wall. In a model file it is defined as
  /// 
  /// @verbatim
  /// Division::VolumeViaLongestWall 3 1 [1]
  /// V_th //L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class VolumeViaLongestWall : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWall(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > 
				 &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };

  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular
  /// to the longest wall segment.
  ///
  /// Divides a cell when volume above a threshold. New wall is created
  /// prependicular to the longest cell wall. In a model file it is defined as
  /// @verbatim
  /// Division::VolumeViaLongestWallCenterTriangulation 3 2 [1] 2
  /// V_th L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// com index, restinglengthIndex
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class VolumeViaLongestWallCenterTriangulation : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWallCenterTriangulation(std::vector<double> &paraValue, 
						    std::vector< std::vector<size_t> > 
						    &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular
  /// to the longest wall segment.
  ///
  /// Divides a cell when volume above a threshold in 3D with centerTriangulation. New wall is created
  /// prependicular to the longest cell wall. In a model file it is defined as
  /// @verbatim
  /// Division::VolumeViaLongestWall3DCenterTriangulation 3 2 [1] 2
  /// V_th L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// com index, restinglengthIndex
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class VolumeViaLongestWall3DCenterTriangulation : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWall3DCenterTriangulation(std::vector<double> &paraValue, 
						      std::vector< std::vector<size_t> > 
						      &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular
  /// to the longest wall segment.
  ///
  /// Creates a branch as
  /// @verbatim
  /// Branching 3 1 [1]
  /// V_th L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class Branching : public BaseCompartmentChange {
    
  public:
    
    Branching(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  
  /// @brief Divides a cell when volume above a threshold and cell close enough to 'apex'
  ///
  /// Divides a cell when volume above a threshold and cell within a distance
  /// from max. New wall is created prependicular to the longest cell wall.
  ///
  class VolumeViaLongestWallSpatial : public BaseCompartmentChange {
    
  private:
    
    double sMax_;
    
  public:
    
    VolumeViaLongestWallSpatial(std::vector<double> &paraValue, 
					std::vector< std::vector<size_t> > 
					&indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold in 3D
  /// Divides a cell when volume above a threshold. Same as
  /// Division::VolumeViaLongestWall but used for surfaces in 3D.  
  class VolumeViaLongestWall3D : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWall3D(std::vector<double> &paraValue, 
				   std::vector< std::vector<size_t> > 
				   &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold in 3D
  /// Divides a cell when volume above a threshold. Same as
  /// Division::VolumeViaLongestWall but used for surfaces in 3D.
  class VolumeViaShortestWall3D : public BaseCompartmentChange {
    
  public:
    
    VolumeViaShortestWall3D(std::vector<double> &paraValue, 
           std::vector< std::vector<size_t> > 
           &indValue );
    
    int flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
    DataMatrix &cellData,
    DataMatrix &wallData,
    DataMatrix &vertexData,
    DataMatrix &cellDerivs,
    DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs );  
  };


  /// @brief Divides a cell when volume above a threshold in 3D and distance to apex smaller than th
  ///
  /// Divides a cell when volume above a threshold. Same as
  /// Division::VolumeViaLongestWall but used for surfaces in 3D.
  ///
  class VolumeViaLongestWall3DSpatial : public BaseCompartmentChange {
    
  private:
    
    double sMax_;
    
  public:
    
    VolumeViaLongestWall3DSpatial(std::vector<double> &paraValue, 
					  std::vector< std::vector<size_t> > 
					  &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold
  /// Divides a cell when volume above a threshold. New wall is created
  ///  prependicular to maximal strain rate.  
  class VolumeViaStrain : public BaseCompartmentChange {
    
  public:
    
    VolumeViaStrain(std::vector<double> &paraValue, 
			    std::vector< std::vector<size_t> > 
			    &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold
  /// Divides a cell when volume above a threshold. New wall is created
  /// prependicular to direction given as cell variable.
  class VolumeViaDirection : public BaseCompartmentChange {
    
  public:
    
    VolumeViaDirection(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold where the threshold
  /// is determined as an increasing Hill function of a concentration inside the
  /// cell.
  ///
  /// Divides a cell when volume above a concentration-determined threshold.
  /// New wall is created
  /// in a random direction through center of mass.
  class VolumeRandomDirectionConcentration : public BaseCompartmentChange {
    
  public:
    
    VolumeRandomDirectionConcentration(std::vector<double> &paraValue, 
				       std::vector< std::vector<size_t> > 
				       &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };

  /// @brief Divides a cell when volume above a threshold
  ///
  /// Divides a cell when volume above a threshold. New wall is created
  ///  in a random direction through center of mass.
  class VolumeRandomDirection : public BaseCompartmentChange {
    
  public:
    
    VolumeRandomDirection(std::vector<double> &paraValue, 
				  std::vector< std::vector<size_t> > 
				  &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief  Divides a cell when volume above a threshold in a random direction
  ///
  /// Divides a cell when volume above a threshold. New wall is created
  /// from a random (wall central) vertex and the second vertex is
  /// chosen such that the two (daughter) areas are as equal as
  /// possible.
  ///
  /// @note This class relies on that a central vertex is defined.
  ///
  class VolumeRandomDirectionCenterTriangulation : public BaseCompartmentChange {
    
  public:
    
    VolumeRandomDirectionCenterTriangulation(std::vector<double> &paraValue, 
						     std::vector< std::vector<size_t> > 
						     &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  class ForceDirection : public BaseCompartmentChange
  {
  public:
    ForceDirection(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
  };
  
  /// @brief Divides a cell when volume above a threshold
  ///
  /// Divides a cell when volume above a threshold. New wall is created at shortest
  /// path that divides the volume in equal parts. 
  class VolumeViaShortestPath : public BaseCompartmentChange {
    
  public:
    
    VolumeViaShortestPath(std::vector<double> &paraValue, 
				  std::vector< std::vector<size_t> > 
				  &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };

  /// @brief Divides a cell (in 2D) along the shortest path through center of mass (or random point).
  ///
  /// @details Divides a cell in 2D when volume above a threshold, @f$V_{threshold}@f$
  /// with new wall created at shortest path that divides the cell through COM
  /// (almost Volume in equal parts), or through a random internal point where p3 can set a weight
  /// for combining both (@f$ p_{3} x_{c} + (1-p_{3}) x_{r} @f$.
  /// Parameters and variable indices are standard and in addition a time
  /// variable can be set and 'read' at division to measure time since previous division.
  /// It is also optional to provide two indices to store the size of the two sisters after division.
  /// In a model file, the reaction is given by 
  /// @verbatim
  /// Division::ShortestPath2D 4 1[2] K [1/3] 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM, [0:1] weighted in between)
  ///
  /// I_k (optional volume (and other variables that should be divided with size) related index to be updated)
  ///
  /// cell_time_index (optional)
  /// sizeSister1Index sizeSister2Index (optional)
  /// @endverbatim
  /// @see Division::ShortestPath for 3D version also applicable for CenterTriangulation
  ///
  class ShortestPath2D : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    ShortestPath2D(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<ShortestPath2D::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
  };

  /// @brief Divides a cell (in 2D) along the shortest path through center of mass (or random point).
  ///
  /// @details Divides a cell in 2D when volume above a threshold, @f$V_{threshold}@f$
  /// with new wall created at shortest path that divides the cell through COM
  /// (almost Volume in equal parts), or through a random internal point.
  /// Parameters and variable indices are standard and in addition a time
  /// variable can be set and 'read' at division to measure time since previous division.
  /// In a model file, the reaction is given by 
  /// @verbatim
  /// Division::ShortestPath2D 7 2 K 1 
  /// V_{threshold}
  /// V_{threshold-max}
  /// K_{hill}
  /// n_{hill}
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  ///
  /// I_k (optional volume (and other variables that should be divided with size) related index to be updated)
  ///
  /// cell concentration index (optional)
  /// @endverbatim
  /// @see Division::ShortestPath for 3D version also applicable for CenterTriangulation
  ///
  class ShortestPath2DConcentration : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    ShortestPath2DConcentration(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<ShortestPath2DConcentration::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
  };

  /// @brief Divides a cell along the shortest path through center of mass (or random point).
  ///
  /// @details Divides a cell when volume above a threshold, with new wall created at shortest
  ///  path that divides the cell through COM (almost Volume in equal parts), or through a
  /// random internal point where p3 can set the weight in between (@f$ p_{3} x_{c} + (1-p_{3}) x_{r} @f$).
  /// Using centerTriangulation and doubleLength 
  /// formats are optional and can be done by setting the corresponding flags. 
  /// @verbatim
  /// Division::ShortestPath 4 1[2] 0/1 1/3 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM, [0:1] weighted in between)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  /// @endverbatim
  /// or
  /// @verbatim
  /// Division::ShortestPath 6 3 0/1 1/3 2 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM, [0:1] weighted inbetween)
  /// centerTriangulation flag (0/1)
  /// double length flag (0/1)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell_time_index(optional)
  /// sister1SizeIndex sister2SizeIndex (optional)
  /// com index 
  /// restinglengthIndex
  ///
  /// @endverbatim
  ///
  /// @see Division::ShortestPath2D for 2D version
  /// @note Should also work for 2D, but needs to be checked
  ///
  class ShortestPath : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    ShortestPath(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<ShortestPath::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
  };

  /// @brief Divides a cell given a Sizer/Timer/Adder rule to decide when to divide and a shortest path for
  /// deciding which plane of division.
  ///  
  /// Divides according to the STA rule and at division selects the shortest path through the cell for 
  /// the new wall.
  ///
  /// @verbatim
  ///
  /// Division::STAViaShortestPath 4 2 0/1 0/1 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index (optional)
  ///
  /// @endverbatim
  ///
  /// or
  ///
  /// @verbatim
  ///
  /// Division::STAViaShortestPath 6 3 0/1 0/1 2 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  /// centerTriangulation flag (0/1)
  /// double length flag (0/1)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index (optional)
  ///
  /// com index 
  /// restinglengthIndex
  ///
  /// @endverbatim  
  class STAViaShortestPath : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    STAViaShortestPath(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<STAViaShortestPath::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
  };

 class ShortestPathGiantCells : public BaseCompartmentChange
 {
 public:
   struct Candidate {
     double distance;
      size_t wall1;
     size_t wall2;
     double px, py;
     double qx, qy;
   };
   
   ShortestPathGiantCells(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > &indValue);
   
   int flag(Tissue *T, size_t i,
	    DataMatrix &cellData,
	    DataMatrix &wallData,
	    DataMatrix &vertexData,
	    DataMatrix &cellDerivs,
	    DataMatrix &wallDerivs,
	    DataMatrix &vertexDerivs);
   
   void update(Tissue* T, size_t i,
	       DataMatrix &cellData,
	       DataMatrix &wallData,
	       DataMatrix &vertexData,
	       DataMatrix &cellDerivs,
	       DataMatrix &wallDerivs,
	       DataMatrix &vertexDerivs);  
   
   std::vector<ShortestPathGiantCells::Candidate> 
     getCandidates(Tissue* T, size_t i,
		   DataMatrix &cellData,
		   DataMatrix &wallData,
		   DataMatrix &vertexData,
		   DataMatrix &cellDerivs,
		   DataMatrix &wallDerivs,
		   DataMatrix &vertexDerivs);
   
   double astar(double sigma, double A, double B);
   double f(double a, double sigma, double A, double B);
 };

 class Random : public BaseCompartmentChange
  {
  public:
    Random(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    // Returns an integer between 0 and n - 1. 
    int random(int n);
  };
  
  class VolumeRandomDirectionGiantCells : public BaseCompartmentChange
  {
  public:
    
    VolumeRandomDirectionGiantCells(std::vector<double> &paraValue, 
					    std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
  };
  
  class MainAxis : public BaseCompartmentChange
  {
  public:
    
    MainAxis(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<double> getMainAxis(Cell &cell, DataMatrix &vertexData);
    
  private:
    struct Candidate
    {
      double s;
      size_t index;
      std::vector<double> p;
    };
    
    class CompareCandidate
    {
    public:
      bool operator()(const Candidate &a, const Candidate &b)
      {
	return std::abs(a.s) > std::abs(b.s);
      }
    };
  };


  /// @brief Division rule that is being used for the ATML1 model, so it is not meant to be used for generic models.
  /// Adapted from Division::ShortestPath
  ///  Divides a cell when a certain flag has value 1, with New wall created at shortest
  ///  path that divides the volume (not!) in equal parts. Using centerTriangulation and doubleLength 
  ///  formats are optional and can be done by setting the corresponding flags. 
  ///  
  ///  Other variables are reset, according to what needs to be done for the ATML1 model.
  ///
  /// @verbatim
  /// Division::FlagResetShortestPath 6 4 0 0 2 1 # 6 pars, 4 types of indices
  /// (ShortestPath was 6 3 0/1 1 2, pending to check if hence FlagResetShortestPath should be 6 4 0/1 1 2) 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  /// centerTriangulation flag (0/1)
  /// double length flag (0/1)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// com index 
  /// restinglengthIndex
  /// flagindex
  /// @endverbatim
  ///
  /// or it might be possible to use (although it has not been tested)

  /// @verbatim
  ///
  /// FlagResetShortestPath 4 3 0/1 1 1 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  ///
  /// @endverbatim
  class FlagResetShortestPath : public BaseCompartmentChange {
    
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
        
    FlagResetShortestPath(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
    std::vector<FlagResetShortestPath::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
  };

  /// @brief UNDER CONSTRUCTION, DO NOT USE YET!!!
  ///
  class FlagResetViaLongestWall : public BaseCompartmentChange {
    
  public:
    FlagResetViaLongestWall(std::vector<double> &paraValue, 
			    std::vector< std::vector<size_t> > 
			    &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
} //end namespace Division
#endif
  
  
