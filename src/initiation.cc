//
// Filename     : initiation.cc
// Description  : Classes describing initiation rules for a tissue simulation
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : August 2019
// Revision     : $Id:$
//
#include"tissue.h"
#include "baseReaction.h"
#include "initiation.h"
#include<cstdlib>

namespace Initiation {

  RandomBoolean::
    RandomBoolean(std::vector<double> &paraValue,
        std::vector< std::vector<size_t> >
        &indValue )
    {
      // Do some checks on the parameters and variable indeces
      //
      if( paraValue.size() != 2 && paraValue.size() != 3) {
        std::cerr << "Initiation::RandomBoolean::"
          << "RandomBoolean() "
          << "Uses first parameter p (probability to set cell variable to one)." << std::endl
          << "The second parameter sets the seed of the random generator." << std::endl
          << "A third parameter can be given to set a deltaT for repeat of the procedure, " << std::endl
          << "but where only off cells can possibly be set to one." << std::endl;
        exit(EXIT_FAILURE);
      }
      if( paraValue[0]<0.0 || paraValue[0]>1.0) {
        std::cerr << "Initiation::RandomBoolean::"
          << "RandomBoolean() "
          << "First parameter a probability and have to be in [0:1]." << std::endl;
        exit(EXIT_FAILURE);
      }

      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
        std::cerr << "Initiation::RandomBoolean::"
          << "RandomBoolean() "
          << "Index for cell variable to be initiated given." << std::endl;
        exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("Initiation::RandomBoolean");
      setParameter(paraValue);
      setVariableIndex(indValue);

      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp[0] = "p";
      tmp[1] = "seed";
      if (numParameter()==3)
        tmp[2] = "DeltaT";

      setParameterId( tmp );

      // Set the next time if deltaT given (assumes starting at t=0)
      localTime_ = 0.0;
      if (numParameter()==3 && parameter(2)>0.0)
        nextTime_ = nextTime_+parameter(2);

    }

  void RandomBoolean::
    derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs )
    {
      // nothing
    }

  void RandomBoolean::
    initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs)
    {
      // Initiate with given seed
      long int idum = long(parameter(1));
      myRandom::sran3(idum);

      //Do the initiation for each cell
      size_t numCells = T.numCell();

      size_t cIndex = variableIndex(0,0);
      double prob = parameter(0);
      //For each cell
      for (size_t cellI = 0; cellI < numCells; ++cellI) {
        cellData[cellI][cIndex] = 0.0;
        if (myRandom::ran3()<prob)
          cellData[cellI][cIndex] = 1.0;
      }
    }

  void RandomBoolean::
    update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &walldata,
        DataMatrix &vertexData,
        double h) 
    {
      localTime_ += h;
      if (numParameter()==3 && parameter(2) > 0.0 && localTime_>nextTime_) {

        //Do the check for each cell
        size_t numCells = T.numCell();

        size_t cIndex = variableIndex(0,0);
        double prob = parameter(0);
        //For each cell
        for (size_t cellI = 0; cellI < numCells; ++cellI) { 
          if (cellData[cellI][cIndex]<0.5 && myRandom::ran3()<prob) { //only off cells checked
            cellData[cellI][cIndex] = 1.0;
          }
        }
        nextTime_ += parameter(2);
      }
    }

    Random::
    Random(std::vector<double> &paraValue,
        std::vector< std::vector<size_t> >
        &indValue )
    {
      // Do some checks on the parameters and variable indeces
      //
      if( paraValue.size() != 2 ) {
        std::cerr << "Initiation::Random::"
          << "Random() "
          << "Uses first parameter maxVal to be multiplied with Random number in [0:1] as initiation." << std::endl
		  << "The second parameter sets the seed of the random generator." << std::endl;
	  exit(EXIT_FAILURE);
      }
      if( paraValue[0]<0.0 ){
        std::cerr << "Initiation::Random::"
          << "Random() "
          << "First parameter needs to be a positive number." << std::endl;
        exit(EXIT_FAILURE);
      }

      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
        std::cerr << "Initiation::Random::"
          << "Random() "
          << "Index for cell variable to be initiated given." << std::endl;
        exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("Initiation::Random");
      setParameter(paraValue);
      setVariableIndex(indValue);

      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp[0] = "p";
      tmp[1] = "seed";

      setParameterId( tmp );
    }

  void Random::
    derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs )
    {
      // nothing
    }

  void Random::
    initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs)
    {
      // Initiate with given seed
      long int idum = long(parameter(1));
      myRandom::sran3(idum);

      //Do the initiation for each cell
      size_t numCells = T.numCell();

      size_t cIndex = variableIndex(0,0);
      double maxRate = parameter(0);
      //For each cell
      for (size_t cellI = 0; cellI < numCells; ++cellI) {
        cellData[cellI][cIndex] = maxRate*myRandom::ran3();
      }
    }

  void Random::
    update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &walldata,
        DataMatrix &vertexData,
        double h) 
    {
      // Nothing
    }

  RandomBooleanBiased::
    RandomBooleanBiased(std::vector<double> &paraValue,
        std::vector< std::vector<size_t> >
        &indValue )
    {
      // Do some checks on the parameters and variable indeces
      //
      if( paraValue.size() != 3 && paraValue.size() != 4) {
        std::cerr << "Initiation::RandomBooleanBiased::"
          << "RandomBooleanBiased() "
          << "Uses first parameter p (probability to set cell variable to one)." << std::endl
          << "The second parameter set the bias ." << std::endl
          << "The third parameter sets the seed of the random generator." << std::endl
          << "A fourth parameter can be given to set a deltaT for repeat of the procedure, " << std::endl
          << "but where only off cells can possibly be set to one." << std::endl;
        exit(EXIT_FAILURE);
      }
      if( paraValue[0]<0.0 || paraValue[0]>1.0) {
        std::cerr << "Initiation::RandomBooleanBiased::"
          << "RandomBooleanBiased() "
          << "First parameter a probability and have to be in [0:1]." << std::endl;
        exit(EXIT_FAILURE);
      }

      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
        std::cerr << "Initiation::RandomBooleanBiased::"
          << "RandomBooleanBiased() "
          << "Index for cell variable to be initiated given." << std::endl;
        exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("Initiation::RandomBooleanBiased");
      setParameter(paraValue);
      setVariableIndex(indValue);

      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp[0] = "p";
      tmp[1] = "N_t";
      tmp[2] = "seed";
      if (numParameter()==4)
        tmp[3] = "DeltaT";

      setParameterId( tmp );

      // Set the next time if deltaT given (assumes starting at t=0)
      localTime_ = 0.0;
      if (numParameter()==4 && parameter(3)>0.0)
        nextTime_ = nextTime_+parameter(3);

    }

  void RandomBooleanBiased::
    derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs )
    {
      // nothing
    }

  void RandomBooleanBiased::
    initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs)
    {
      // Initiate with given seed
      long int idum = long(parameter(2));
      myRandom::sran3(idum);

      //Do the initiation for each cell
      size_t numCells = T.numCell();

      size_t cIndex = variableIndex(0,0);
      double prob = parameter(0);
      //For each cell
      for (size_t cellI = 0; cellI < numCells; ++cellI) {
        // Count number of neighs on
        size_t numOn=0;
        for (size_t k=0; k<T.cell(cellI).numWall(); k++) {
          size_t neighI = T.cell(cellI).cellNeighbor(k)->index();
          if (T.cell(cellI).cellNeighbor(k) != T.background() && cellData[neighI][cIndex]>0.5)
            numOn++;
        }
        cellData[cellI][cIndex] = 0.0;
        if (numOn<=parameter(1) && myRandom::ran3()<prob)
          cellData[cellI][cIndex] = 1.0;
      }
    }

  void RandomBooleanBiased::
    update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &walldata,
        DataMatrix &vertexData,
        double h) 
    {
      localTime_ += h;
      if (numParameter()==4 && parameter(3) > 0.0 && localTime_>nextTime_) {

        //Do the check for each cell
        size_t numCells = T.numCell();

        size_t cIndex = variableIndex(0,0);
        double prob = parameter(0);
        //For each cell
        for (size_t cellI = 0; cellI < numCells; ++cellI) { 
          if (cellData[cellI][cIndex]<0.5) {
            // Count number of neighs on
            size_t numOn=0;
            for (size_t k=0; k<T.cell(cellI).numWall(); k++) {
              size_t neighI = T.cell(cellI).cellNeighbor(k)->index();
              if (T.cell(cellI).cellNeighbor(k) != T.background() && cellData[neighI][cIndex]>0.5)
                numOn++;
            }
            if (numOn<=parameter(1) && myRandom::ran3()<prob) { //only off cells checked
              cellData[cellI][cIndex] = 1.0;
            }
          }
        }
        nextTime_ += parameter(3);
      }
    }

  FaceArea2D::
    FaceArea2D(std::vector<double> &paraValue,
        std::vector< std::vector<size_t> >
        &indValue )
    {
      // Do some checks on the parameters and variable indeces
      if (paraValue.size() != 0 && paraValue.size() != 1) {
        std::cerr << "Initiation::FaceArea2D::"
          << "FaceArea2D() "
          << "takes a single optional parameter." << std::endl;
        exit(EXIT_FAILURE);
      }

      if (indValue.size() != 1) {
        std::cerr << "Initiation::FaceArea2D::"
          << "FaceArea2D() "
          << "Index for cell variable representing cell face area "
          << "to be initiated given." << std::endl;
        exit(EXIT_FAILURE);
      }

      // Set the variable values
      setId("Initiation::FaceArea2D");
      setParameter(paraValue);
      setVariableIndex(indValue);

      // Set the parameter identities

      std::vector<std::string> tmp(numParameter());
      if (numParameter() == 1) {
        tmp[0] = "updateFlag";
      }
      setParameterId(tmp);
    }

  void FaceArea2D::
    derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs )
    {
      // Do nothing
    }

  void FaceArea2D::
    initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs)
    {
      size_t numCells = T.numCell();
      size_t cellVolumeIdx = variableIndex(0, 0);

      for (size_t ii = 0; ii < numCells; ii++) {
        cellData[ii][cellVolumeIdx] = T.cell(ii).calculateVolume();
      }
    }

  void FaceArea2D::
    update(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &walldata,
        DataMatrix &vertexData,
        double h) 
    {
      if (numParameter() == 1 && parameter(0) > 0.0) {
        size_t numCells = T.numCell();
        size_t cellVolumeIdx = variableIndex(0, 0);

        for (size_t ii = 0; ii < numCells; ii++) { 
          cellData[ii][cellVolumeIdx] = T.cell(ii).calculateVolume();
        }
      }
    }

  CellLabel::
    CellLabel(std::vector<double> &paraValue,
        std::vector< std::vector<size_t> >
        &indValue )
    {
      // Do some checks on the parameters and variable indeces
      if (paraValue.size() != 0) {
        std::cerr << "Initiation::CellLabel::"
          << "CellLabel() "
          << "uses no parameters." << std::endl;
        exit(EXIT_FAILURE);
      }

      if (indValue.size() != 1) {
        std::cerr << "Initiation::CellLabel::"
          << "CellLabel() "
          << "Index for cell variable set to integer index "
          << "to be initiated given." << std::endl;
        exit(EXIT_FAILURE);
      }

      // Set the variable values
      setId("Initiation::CellIndex");
      setParameter(paraValue);
      setVariableIndex(indValue);
    }

  void CellLabel::
    initiate(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs)
    {
      size_t numCells = T.numCell();
      size_t cellLabelIndex = variableIndex(0, 0);

      for (size_t ii = 0; ii < numCells; ii++) {
        cellData[ii][cellLabelIndex] = ii;
      }
    }

    void CellLabel::
    derivs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs )
    {
      // Do nothing
    }

} // end namespace Initiation

