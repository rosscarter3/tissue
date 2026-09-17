//
// Filename     : directionReaction.h
// Description  : Classes describing some reaction updates related to directions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : May 2008
// Revision     : $Id:$
//
#ifndef DIRECTIONREACTION_H
#define DIRECTIONREACTION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the derivative and moves the
/// direction as a function of the differnce in angle.
///
/// In a model file the reaction is defined as
/// @verbatim
/// ContinousMTDirection 1 2 1 1
/// k_rate
/// target index
/// MT index
/// @endverbatim
///
/// @note It is only implemented for two dimensions.
///
class ContinousMTDirection : public BaseReaction
{
public:
  ContinousMTDirection(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};


///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the derivative and moves the
/// direction based on velocity vector proportional to the subtraction of 
/// target vector and the direction which is going to be updated.
/// Note: it is only implemented for three dimensions.
///
/// In a model file the reaction is defined as
/// @verbatim
/// ContinousMTDirection3d 1 2 1 1
///
/// k_rate
///
/// target index
///
/// MT index
///
/// @endverbatim
///
class ContinousMTDirection3d : public BaseReaction
{
public:
  ContinousMTDirection3d(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the update and moves the
/// direction toward a given direction using an Euler step directly on the direction vector.
/// It normalizes the vector after the update.
///
/// In a model file the reaction is defined as
/// @verbatim
/// UpdateMTDirection 1 2 1 1
/// k_rate
/// target index
/// MT index
/// @endverbatim
///
class UpdateMTDirection : public BaseReaction
{
 public:
  UpdateMTDirection(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
	void initiate(Tissue &T,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
		      
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
		    double h);
};

///
/// @brief Updates a direction (MT) continuosly towards a target direction with a given rate, where the update
/// can be blocked for cells using a vertex velocity threshold or a stress difference threshold.
///
/// @details This reaction updates an direction (MT) towards a given target direction using an Euler step directly
/// on the direction vector. In its simples form it only requires a rate parameter for the update and
/// indices for where the target and updated directions are stored in the cellData, e.g.
/// @verbatim
/// UpdateMTDirectionEquilibrium 1 2 1 1 
/// k_rate
///
/// target_index # (read)
/// MT_index     # (updated)
/// @endverbatim
/// where it assumed the indices marks the start index of vectors (same size as dimension).
/// It normalizes the MT vector after the update.
/// Before the simulation starts, the MT direction is initiated to the target direction.
/// In addition, a threshold can be set such that the update is done only in the cells that are close "enough" to
/// equilibrium. Equilibrium is considered as a state in which a cell's vertex velocities (e.g. from reaction
/// Calculate::VertexVelocity) are less than a defined threshold. 
/// Then the reaction is defined as
/// @verbatim
/// UpdateMTDirectionEquilibrium 2 3 1 1 1 
/// k_rate
/// velocity_threshold
///
/// target_index
/// MT_index
/// velocity_index
/// @endverbatim
/// where the update only happens in cells where the value read at the velocityIndex is less then the velocityThreshold.
/// Similarly, a threshold for only updating if the (absolute) relative difference between two stresses
/// (e.g. Principal and vonMises) in cells is larger than a threshold. Then the reaction is defined as
/// @verbatim
/// UpdateMTDirectionEquilibrium 3 4 1 1 1 2 
/// k_rate
/// velocity_threshold
/// stressDiff_threshold
///
/// target_index
/// MT_index
/// velocity_index
/// Stress1_index Stress2_index
/// @endverbatim
/// where the alternative stresses are read from cell variables specified.
///
/// @note This function applies the update within the update function (i.e. between derivative steps).
/// @see Calculate::VertexVelocity calculates average vertex velocities for cells and stores them.
/// @see VertexFromTRBScenterTriangulationMT and similar for calculating and storing different stresses.
///
class UpdateMTDirectionEquilibrium : public BaseReaction
{
 public:
  UpdateMTDirectionEquilibrium(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
	void initiate(Tissue &T,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
	
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
                    double h);
};


///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the update and moves the
/// direction toward a given direction using an Euler step directly on the direction vector.
/// the update rate depends on both a user defined rate and a Hill function of a concentration 
/// that can be stress or strain anisotropy. It normalizes the vector after the update.
///
/// In a model file the reaction is defined as
/// @verbatim
/// UpdateMTDirectionConcenHill 3 3 1 1 1
/// k_rate
/// k_Hill
/// n_hill
/// target index
/// MT index
/// concentration(anisotropy) index
/// @endverbatim
///
class UpdateMTDirectionConcenHill : public BaseReaction
{
 public:
  UpdateMTDirectionConcenHill(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
	void initiate(Tissue &T,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
		      
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
		    double h);
};




class RotatingDirection : public BaseReaction
{
 public:
  RotatingDirection(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

#endif //DIRECTIONREACTION_H
