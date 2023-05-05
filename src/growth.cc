//
// Filename     : growth.cc
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include <vector>
#include <cmath>
#include"growth.h"
#include"baseReaction.h"

namespace WallGrowth {
  
  Constant::
  Constant(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 && paraValue.size() !=3) {
      std::cerr << "WallGrowth::Constant::"
		<< "Constant() "
		<< "Two or three parameters used  k_growth, linearFlag, [L_trunc]"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "WallGrowth::Const::"
		<< "Const() "
		<< "One variable index is used for specifying resting length."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Constant");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "linearFlag";
    if (numParameter()>2) {
      tmp[2] = "L_trunc";
    }
    setParameterId( tmp );
  }
  
  void Constant::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      double arg = parameter(0);
      if (parameter(1)==1) {//linearFlag
	arg *= wallData[i][lengthIndex];
      }
      if (numParameter()>2) {//truncated at maximal length
	arg *= (1 - wallData[i][lengthIndex]/parameter(2));
      }
      wallDerivs[i][lengthIndex] += arg;
    }
  }

  Stress::
  Stress(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=4 && paraValue.size()!=5 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "Uses four or five parameters k_growth, stress_threshold "
		<< "stretch_flag (0 for stress (read from wall variable),1 for strain),"
		<< " linear_flag (0 const, 1 prop to wall length), and [L_threshold]." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[2] != 0.0 && paraValue[2] != 1.0 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 1 
	|| (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "One variable index is used (wall length index) at first "
		<< "level, and stress variable index/indices at second (if strain_flag not set)."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Stress");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "s_threshold";
    tmp[2] = "strain_flag";
    tmp[3] = "linear_flag";
    if (numParameter()>4) {
      tmp[4] = "L_threshold";
    }
    setParameterId( tmp );
  }
  
  void Stress::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(2)) {//Stress used, read from saved data in the wall
	      for (size_t k=0; k<numVariableIndex(1); ++k)
	        stress += wallData[i][variableIndex(1,k)];
      }
      else { //Strain/stretch used
      	double distance=0.0;
      	for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
      	  distance += (vertexData[v1][d]-vertexData[v2][d])*
      	    (vertexData[v1][d]-vertexData[v2][d]);
      	  distance = std::sqrt(distance);
      	  stress = (distance-wallData[i][lengthIndex]) /
      	    wallData[i][lengthIndex];
      }
      if (parameter(1)==0.0 || stress > parameter(1)) {
      	double growthRate = parameter(0)*(stress - parameter(1));
      	if (parameter(3))
      	  growthRate *= wallData[i][lengthIndex];
      	if (numParameter()>4) {
      	  growthRate *= (1.0 - wallData[i][lengthIndex]/parameter(4));
	      }
	      wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }


  void Stress::
  derivsWithAbs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs,
        DataMatrix &sdydtCell,
        DataMatrix &sdydtWall,
        DataMatrix &sdydtVertex){
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(2)) {//Stress used, read from saved data in the wall
        for (size_t k=0; k<numVariableIndex(1); ++k)
          stress += wallData[i][variableIndex(1,k)];
      }
      else { //Strain/stretch used
        double distance=0.0;
        for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
          distance += (vertexData[v1][d]-vertexData[v2][d])*
            (vertexData[v1][d]-vertexData[v2][d]);
        distance = std::sqrt(distance);
        stress = (distance-wallData[i][lengthIndex]) /
          wallData[i][lengthIndex];
      }
      if (parameter(1)==0.0 || stress > parameter(1)) {
        double growthRate = parameter(0)*(stress - parameter(1));
        if (parameter(3))
          growthRate *= wallData[i][lengthIndex];
        if (numParameter()>4) {
          growthRate *= (1.0 - wallData[i][lengthIndex]/parameter(4));
        }
      wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }

  Strain::
  Strain(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=4 && paraValue.size()!=5) {
      std::cerr << "WallGrowth::Strain::"
		<< "Strain() "
		<< "Uses four parameters k_growth, strain_threshold "
		<< "strain_flag (0 for almansi strain),"
		<< "linear_flag (0 const, 1 prop to wall length) "
                << "fifth parameter if set is the velocity threshold for update "
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[2] != 0.0 ) {
      std::cerr << "WallGrowth::Strain::"
		<< "Strain() "
		<< "strain_flag parameter must be 0 (almansi strain) or " 
		<< "1 (other strain measures are not implemented)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
      std::cerr << "WallGrowth::Strain::"
		<< "Strain() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( indValue.size()!=1 || (indValue[0].size() != 1 && indValue[0].size() != 2) ) {
      std::cerr << "WallGrowth::Strain::"
		<< "Strain() "
		<< "One variable index is used (wall length index) at first level"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Strain");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "s_threshold";
    tmp[2] = "strain_flag";
    tmp[3] = "linear_flag";
    if (numParameter()==5)
      tmp[3] = "velocity_threshold";
    
    setParameterId( tmp );
  }
  
  void Strain::
  derivs(Tissue &T,
         DataMatrix &cellData,
         DataMatrix &wallData,
         DataMatrix &vertexData,
         DataMatrix &cellDerivs,
         DataMatrix &wallDerivs,
         DataMatrix &vertexDerivs ) {
  }
  
  void Strain::
  update(Tissue &T,
         DataMatrix &cellData,
         DataMatrix &wallData,
         DataMatrix &vertexData,
         double h ) {

    size_t numWalls = T.numWall();
    size_t numCells = T.numCell();
    double velocityThreshold=parameter(4); 
    size_t velocityStoreIndex =variableIndex(0,1);
    size_t lengthIndex = variableIndex(0,0);

    static size_t numEquil=0;
    static size_t numSteps=0;
    numSteps++;
    
    static double totaltime=0;
    static double deltat=0;
    deltat +=h;
    totaltime+=h;  
    bool equil=true;
    for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
      if(cellData[cellIndex][velocityStoreIndex]>velocityThreshold)
        equil=false;    
    
    if(equil) {// if close to mechanical equilibrum
      double totalStrain=0;
      
      if (parameter(2)==0) //Almansi Strain used
        for( size_t i=0 ; i<numWalls ; ++i ) {
	  
          size_t v1 = T.wall(i).vertex1()->index();
          size_t v2 = T.wall(i).vertex2()->index();
          double restingL=wallData[i][variableIndex(0,0)];
          double strain=0.0;
          double distance2=0.0;
          for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
            distance2 += (vertexData[v1][d]-vertexData[v2][d])*
              (vertexData[v1][d]-vertexData[v2][d]);
          
          distance2=std::sqrt(distance2);
          strain = 0.5*(1-(restingL/distance2)*(restingL/distance2));
          totalStrain+=strain;
	  
          if (numParameter()!=5)              
            if (strain > parameter(1)) {
	      double growthRate = h*restingL*parameter(0)*(strain - parameter(1));
              wallData[i][lengthIndex] += growthRate;
            }
          if (numParameter()==5){ 
            size_t c1,c2;          
            if(T.wall(i).cell1()!=T.background())
              c1= T.wall(i).cell1()->index();
            else
              c1= T.wall(i).cell2()->index();
            if(T.wall(i).cell2()!=T.background())   
              c2= T.wall(i).cell2()->index();
            else
              c2= T.wall(i).cell1()->index();
            
            //double velocity=0.5*(cellData[c1][variableIndex(0,1)]+cellData[c2][variableIndex(0,1)]);
            if (strain > parameter(1) ) {	      
              double factor1=parameter(0)*h*(strain-parameter(1))/(1-2*strain);
              factor1=factor1+factor1*factor1;// first and second terms in taylor expansion              
              //double factor1=parameter(0)*h*(strain-parameter(1));
              wallData[i][lengthIndex] += factor1*restingL;              
            }
          }          
        }
      numEquil++;
      //std::cout<<"walls  "<<totalStrain/numWalls<<std::endl;          
      deltat =0;       
    }
  }
  
  namespace CenterTriangulation {

    Constant::
    Constant(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=2 && paraValue.size() !=3) {
	      std::cerr << "WallGrowth::CenterTriangulation::Constant::"
		    << "Constant() "
		    << "Two or three parameters used  k_growth, linearFlag, [L_trunc]"
		    << std::endl;
	      exit(EXIT_FAILURE);
      }
      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
	      std::cerr << "WallGrowth::CenterTriangulation::Const::"
		    << "Const() "
		    << "Start of additional Cell variable indices (center(x,y,z) "
		    << "L_1,...,L_n, n=num vertex) is given in first level." 
		    << std::endl;
	      exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowth::CenterTriangulation::Constant");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "linearFlag";
      if (numParameter()>2) {
	      tmp[2] = "L_trunc";
      }
      setParameterId( tmp );
    }
    
    void Constant::
    derivs(Tissue &T, 
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t lengthIndex = variableIndex(0,0);
      size_t lengthStartIndex = lengthIndex+3;//assuming 3D
      
      //std::cerr << "====================" << std::endl;

      //std::cerr << "wallData[0]: " << wallData[0][0] << std::endl;
      //for (size_t j=0; j < T.cell(0).numVertex(); ++j) {
      //  std::cerr << j+lengthStartIndex << ": " << cellData[0][j+lengthStartIndex] << std::endl;
      //}


      for (size_t i=0; i<numCells; ++i) {
	           for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
	        double arg = parameter(0);
	        if (parameter(1)==1) {//linearFlag (prop to length)
	          arg *= cellData[i][k+lengthStartIndex];
	        }
	        if (numParameter()>2) {//truncated at maximal length
	          arg *= (1 - cellData[i][k+lengthStartIndex]/parameter(2));
	        }
	        cellDerivs[i][k+lengthStartIndex] += arg;
          // std::cerr << "arg: " << arg << std::endl;
	      }
      }
    }
    
    Stress::
    Stress(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=4 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "Uses four parameters k_growth, stress_threshold "
		  << "stretch_flag and linear_flag (0 const, 1 prop to wall length)" 
		  << std::endl;
	exit(0);
      }
      if( paraValue[2] != 0.0 && paraValue[2] != 1.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "stretch_flag parameter must be 0 (stress used) or " 
		  << "1 (stretch used)." << std::endl;
	exit(0);
      }
      if( paraValue[2] == 0.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "stretch_flag parameter must be 1 (stretch used) (not implemented for" 
		  << " stress yet..." 
		  << std::endl;
	exit(0);
      }
      if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "linear_flag parameter must be 0 (constant growth) or " 
		  << "1 (length dependent growth)." << std::endl;
	exit(0);
      }
      
      if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 1 
	  || (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in first level, " 
		  << "and stress variable indices at second (if strain_flag not set)."
		  << std::endl;
	exit(0);
      }
      //Set the variable values
      //
      setId("WallGrowthStresscenterTriangulation");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "s_threshold";
      tmp[2] = "strain_flag";
      tmp[3] = "linear_flag";
      setParameterId( tmp );
    }
    
    void Stress::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t posStartIndex = variableIndex(0,0);
      size_t lengthStartIndex = posStartIndex+3;
      
      for (size_t i=0; i<numCells; ++i) {
	for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
	  size_t v = T.cell(i).vertex(k)->index();
	  double stress=0.0;
	  if (!parameter(2)) {//Stress used, read from saved data in the wall
	    std::cerr << "WallGrowthStresscenterTriangulation::derivs() " << std::endl
		      << "Strain (and not stress) is the only implemented version sofar."
		      << std::endl;
	    //for (size_t k=0; k<numVariableIndex(1); ++k)
	    //stress += wallData[i][variableIndex(1,k)];
	  }
	  else { //Strain/stretch used
	    double distance=0.0;
	    for( size_t d=0 ; d<vertexData[v].size() ; d++ )
	      distance += (vertexData[v][d]-cellData[i][d+posStartIndex])*
		(vertexData[v][d]-cellData[i][d+posStartIndex]);
	    distance = std::sqrt(distance);
	    stress = (distance-cellData[i][k+lengthStartIndex]) /
	      cellData[i][k+lengthStartIndex];
	  }
	  if (stress > parameter(1)) {
	    double growthRate = parameter(0)*(stress - parameter(1));
	    if (parameter(3))
	      growthRate *= cellData[i][k+lengthStartIndex];
	    cellDerivs[i][k+lengthStartIndex] += growthRate;
	  }
	}
      }
    }

    StressConcentrationHill::
    StressConcentrationHill(std::vector<double> &paraValue, 
			    std::vector< std::vector<size_t> > 
			    &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=7 ) {
	std::cerr << "WallGrowth::CenterTriangulation::StressConcentrationHill"
		  << "StressConcentrationHill() "
		  << "Uses seven parameters k_growthConst, k_growthHill, K_Hill, n_Hill and "
		  << "stress_threshold, stretch_flag (1 for now), linear_flag (0 const, 1 prop to wall length)" 
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
	std::cerr << "WallGrowth::CenterTriangulation::StressConcentrationHill"
      << "StressConcentrationHill() "
		  << "stretch_flag parameter must be 0 (stress used) or " 
		  << "1 (stretch used)." << std::endl;
	exit(EXIT_FAILURE);
      }
      if( paraValue[5] == 0.0 ) {
	std::cerr << "WallGrowth::CenterTriangulation::StressConcentrationHill"
      << "StressConcentrationHill() "
		  << "stretch_flag parameter must be 1 (stretch used) (not implemented for" 
		  << " stress yet..." 
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      if( paraValue[6] != 0.0 && paraValue[6] != 1.0 ) {
	std::cerr << "WallGrowth::CenterTriangulation::StressConcentrationHill"
      << "StressConcentrationHill() "
		  << "linear_flag parameter must be 0 (constant growth) or " 
		  << "1 (length dependent growth)." << std::endl;
	exit(EXIT_FAILURE);
      }
      
      if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 2 
	  || (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
	std::cerr << "WallGrowth::CenterTriangulation::StressConcentrationHill"
      << "StressConcentrationHill() "
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) and cellIndex for concentration are given in first level, " 
		  << "and stress variable indices at second (if strain_flag not set)."
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowth::CenterTrangulation::StressConcentrationHill");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growthConst";
      tmp[1] = "k_growthHill";
      tmp[2] = "K_Hill";
      tmp[3] = "n_Hill";
      tmp[4] = "s_threshold";
      tmp[5] = "strain_flag";
      tmp[6] = "linear_flag";
      
      setParameterId( tmp );
    }
    
    void StressConcentrationHill::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t posStartIndex = variableIndex(0,0);
      size_t lengthStartIndex = posStartIndex+3;
      size_t concIndex = variableIndex(0,1); 
      double Kpow = std::pow(parameter(2),parameter(3));

      
      for (size_t i=0; i<numCells; ++i) {
      	// Get Hill factor for concentration
	      double concpow = std::pow(cellData[i][concIndex],parameter(3));
	      double hillFactor = concpow/(Kpow+concpow); 

      	for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
      	  size_t v = T.cell(i).vertex(k)->index();
      	  double stress=0.0;
      	  if (!parameter(5)) {//Stress used, read from saved data in the wall
      	    std::cerr << "WallGrowth::CenterTriangulation::StressConcentrationHill::derivs() " << std::endl
      		      << "Strain (and not stress) is the only implemented version sofar."
      		      << std::endl;
      	    //for (size_t k=0; k<numVariableIndex(1); ++k)
      	    //stress += wallData[i][variableIndex(1,k)];
      	  }
      	  else { //Strain/stretch used
      	    double distance=0.0;
      	    for( size_t d=0 ; d<vertexData[v].size() ; d++ )
      	      distance += (vertexData[v][d]-cellData[i][d+posStartIndex])*
      		      (vertexData[v][d]-cellData[i][d+posStartIndex]);
      	      distance = std::sqrt(distance);
      	      stress = (distance-cellData[i][k+lengthStartIndex]) /
      	        cellData[i][k+lengthStartIndex];
      	  }
      	  if (stress > parameter(4)) {
      	    double growthRate = (parameter(0)+parameter(1)*hillFactor)*(stress - parameter(4));
      	    if (parameter(6))
      	      growthRate *= cellData[i][k+lengthStartIndex];
      	    cellDerivs[i][k+lengthStartIndex] += growthRate;
      	  }
      	}
      }
    }

    StrainTRBS::
    StrainTRBS(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {

      //Do some checks on the parameters and variable indeces
      if( paraValue.size()!=3 && paraValue.size()!=4 && paraValue.size()!=6 ) {
	std::cerr << "WallGrowthStrainTRBScenterTriangulation::"
		  << "WallGrowthStrainTRBScenterTriangulation() "
		  << "Uses three parameters k_growth, strain_threshold and velocity threshold "
                  << "- or four parameters with the forth one equal to one " 
                  << "for the case of using two independent resting lengths for neighbours."
                  << "- or six parameters with the fifth one the K value for a concentration "
                  << "Hill function with n=2 and the sixth parameter v_max"
                  << std::endl;
	exit(EXIT_FAILURE);
      }
        
      if( indValue.size() !=3 || 
          indValue[0].size() !=1 || 
          indValue[1].size() !=1 || 
          (indValue[2].size() !=1 && indValue[2].size() !=2) ) {
	std::cerr << "WallGrowthStrainTRBScenterTriangulation::"
		  << "WallGrowthStrainTRBScenterTriangulation() "
                  << "wall length index is given in first level,"
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in second level, " 
		  << "and velocity_value_strore_index and potentially concentration "
                  << "index for domain dependent growth threshold in third level."
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowthStrainTRBScenterTriangulation");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "s_threshold";
      tmp[2] = "velocity_threshold";
      if(numParameter()>3)
        tmp[3] = "doubleFlag";
      if(numParameter()>4){
        tmp[4] = "k_hill";
        tmp[5] = "v_hill";
      }
      setParameterId( tmp );
    }
    
    void StrainTRBS::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
    }

    void StrainTRBS::
    update(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
           double h )
    {
      size_t dimension = vertexData[0].size();
      size_t numCells = T.numCell();
      size_t numWalls = T.numWall();
      size_t wallLengthIndex= variableIndex(0,0);
      size_t comIndex = variableIndex(1,0);
      size_t lengthInternalIndex = comIndex+dimension;
      size_t velocityStoreIndex = variableIndex(2,0);
      double strainThreshold=parameter(1);
      double velocityThreshold=parameter(2);
      size_t growthInd=30;

      static double growthtime=0;
      static double deltat=0;
      deltat +=h; 
            
      bool equil=true;
      for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){        
        if(cellData[cellIndex][velocityStoreIndex]>velocityThreshold)
          equil=false;
      }
      //std::cerr << equil << std::endl; //to check if it is ever updating
      //if(equil && deltat>0.0) {// If close to mechanical equilibrum        
      if(true) {  // If no check to be done
        growthtime+=h;  
        //std::cerr<<scaletmp<<std::endl;
        std::vector<std::vector<double> > mainWalls(numWalls);
        std::vector<std::vector<std::vector<double> > > internalWalls(numCells);
        
        for (size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex)
          mainWalls[wallIndex].resize(2);
        for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){
          size_t numCellWalls = T.cell(cellIndex).numWall();
          internalWalls[cellIndex].resize(numCellWalls);
          for (size_t cellWallIndex=0 ; cellWallIndex<numCellWalls ; ++cellWallIndex)
            internalWalls[cellIndex][cellWallIndex].resize(2);
        }
        
        for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
          //if(cellData[cellIndex][37]!=0) { // hypocotyl             
          
          size_t numCellWalls = T.cell(cellIndex).numWall(); 
          
          cellData[cellIndex][growthInd]=0;
          
          for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
            size_t wallIndexPlusOneMod = (wallIndex+1)%numCellWalls;
            //size_t v1 = com;
            size_t v2 = T.cell(cellIndex).vertex(wallIndex)->index();
            size_t v3 = T.cell(cellIndex).vertex(wallIndexPlusOneMod)->index();
            //size_t w1 = internal wallIndex
            size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();
            //size_t w3 = internal wallIndex+1
            
            DataMatrix position(3,vertexData[v2]);
            for (size_t d=0; d<dimension; ++d)
              position[0][d] = cellData[cellIndex][comIndex+d]; // com position
            //position[1] = vertexData[v2]; // given by initiation
            position[2] = vertexData[v3];
                        
            std::vector<double> restingLength(3);
            if(numParameter()>3 && parameter(3)==1){ // double resting length        
              restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1];
              restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod];      
              
              restingLength[1] = wallData[w2][wallLengthIndex]+
                cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex];
            }
            else{           // single resting length
              restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallIndex];
              restingLength[2] = cellData[cellIndex][lengthInternalIndex + wallIndexPlusOneMod];      
              restingLength[1] = wallData[w2][wallLengthIndex];
            } 
            
            // double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
            //                               (-restingLength[0]+restingLength[1]+restingLength[2])*
            //                               ( restingLength[0]-restingLength[1]+restingLength[2])*
            //                               ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
	    std::vector<double> length(3);
            length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                                   (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                                   (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );            
            length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
	    length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                                   (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                                   (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
            
            //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
            double CurrentAngle1=std::acos(  (length[0]*length[0]+
                                              length[1]*length[1]-
                                              length[2]*length[2])/
                                             (length[0]*length[1]*2)    );
            
            double Qa=std::cos(CurrentAngle1)*length[0];
            double Qc=std::sin(CurrentAngle1)*length[0];
            double Qb=length[1];
                        
            double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+
                                              restingLength[1]*restingLength[1]-
                                              restingLength[2]*restingLength[2])/
                                             (restingLength[0]*restingLength[1]*2)    );
            
            double Pa=std::cos(RestingAngle1)*restingLength[0];
            double Pc=std::sin(RestingAngle1)*restingLength[0];
            double Pb=restingLength[1];
            
            // shape vector matrix in resting shape in local coordinate system  = 
            // inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
            double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                              {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                              { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
            
            double positionLocal[3][2]={ {Qa , Qc}, 
                                         {0  , 0 },  
                                         {Qb , 0 }  };
            
            double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
            for ( int ii=0 ; ii<3 ; ++ii ) {
              DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
              DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
              DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
              DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
            }
            
            size_t growthVecInd=0;
            size_t growthVal1Ind=3;
            size_t growthVal2Ind=7;

            // size_t growthVecInd=32;
            // size_t growthVal1Ind=35;
            // size_t growthVal2Ind=15;

            double growthVal1=cellData[cellIndex][growthVal1Ind];
            double growthVal2=cellData[cellIndex][growthVal2Ind];


            double growthCurrGlob[3]={cellData[cellIndex][growthVecInd  ],
                                      cellData[cellIndex][growthVecInd+1],
                                      cellData[cellIndex][growthVecInd+2] };
            
            
            double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                                   (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                                   (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );
            
            double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                                   (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                                   (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
            
      
            double Xcurrent[3]=      
              { (position[2][0]-position[1][0])/tempA,
                (position[2][1]-position[1][1])/tempA,
                (position[2][2]-position[1][2])/tempA
              };
            
            double Bcurrent[3]=      
              { (position[0][0]-position[1][0])/tempB,
                (position[0][1]-position[1][1])/tempB,
                (position[0][2]-position[1][2])/tempB
              };
            
            double Zcurrent[3]=      
              { Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1],
                Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2],
                Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0]
              };
            
            tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
            Zcurrent[0]=Zcurrent[0]/tempA;
            Zcurrent[1]=Zcurrent[1]/tempA;
            Zcurrent[2]=Zcurrent[2]/tempA;
            
            double Ycurrent[3]=      
              { Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1],
                Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2],
                Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0]
              };
            
            double rotation[3][3]=
              { {Xcurrent[0] , Ycurrent[0] , Zcurrent[0] },
                {Xcurrent[1] , Ycurrent[1] , Zcurrent[1] },
                {Xcurrent[2] , Ycurrent[2] , Zcurrent[2] } };
      
            // rotating the growth vector from global coordinate system to the local one in the current shape
            double growthCurrLocal[3]=
              {  rotation[0][0]*growthCurrGlob[0]+
                 rotation[1][0]*growthCurrGlob[1]+
                 rotation[2][0]*growthCurrGlob[2],
                 
                 rotation[0][1]*growthCurrGlob[0]+
                 rotation[1][1]*growthCurrGlob[1]+
                 rotation[2][1]*growthCurrGlob[2],
                 
                 rotation[0][2]*growthCurrGlob[0]+
                 rotation[1][2]*growthCurrGlob[1]+
                 rotation[2][2]*growthCurrGlob[2]
              };

            // transform the growth vector from current local to the resting local
            double growthRestLocal[2]={ DeformGrad[1][1] *growthCurrLocal[0]-DeformGrad[0][1]*growthCurrLocal[1],
                                        -DeformGrad[1][0]*growthCurrLocal[0]+DeformGrad[0][0]*growthCurrLocal[1]};
            double tmpG=std::sqrt(growthRestLocal[0]*growthRestLocal[0]+growthRestLocal[1]*growthRestLocal[1]);
            growthRestLocal[0]/=tmpG;
            growthRestLocal[1]/=tmpG;
            
            std::vector<std::vector<double> > edgeRestLocal(3);
	    for (size_t d=0; d< 3; ++d)
              edgeRestLocal[d].resize(3);
                        
            edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
            edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
	    edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
            edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];            
            edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
            edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
            
            std:: vector<double> cosTet(3);
            std:: vector<double> sinTet(3);
            for (size_t j=0; j< 3; ++j){
              cosTet[j]=std::fabs((growthRestLocal[0]*edgeRestLocal[j][0]+
                                   growthRestLocal[1]*edgeRestLocal[j][1])/
                                  restingLength[j]);
              sinTet[j]=std::sqrt(std::fabs(1-cosTet[j]*cosTet[j]));
              
              if(cosTet[j]<0 || cosTet[j]>1)
                std::cerr<<"in growth.cc/strainTRBS cosTet is wrong: "
                         <<cosTet[j]<<" in cell here   ..."<<cellIndex<<std::endl;
            }
            
            std::vector<std::vector<double> > restingComp(3);
            for (size_t j=0; j< 3; ++j)
              restingComp[j].resize(2);
            
            for (size_t j=0; j< 3; ++j){
              restingComp[j][0]=restingLength[j]*cosTet[j];
              restingComp[j][1]=restingLength[j]*sinTet[j];
            }
            
            // for hypocotyl Siobhan
            // if(cellIndex==408 && wallIndex==0){
            //     // std::cout<<strainThreshold<<" "<<growthVal1<<" "<<growthVal2<<" "
            //     //          <<growthtime<<" "<<deltat<<" "
            //     //          <<vertexData[0][2]-vertexData[899][2]<<std::endl;
            //     std::cout<<strainThreshold<<" "<<growthVal1<<" "<<growthVal2<<" "
            //              <<(growthVal1-growthVal2)/growthVal1<<std::endl;
                            
            // }
            
            if (growthVal1>strainThreshold && growthVal2<strainThreshold){
	      cellData[cellIndex][growthInd]+=1;              
              double factor1=parameter(0)*h*(growthVal1-strainThreshold);
                                          
              for (size_t j=0; j< 3; ++j)
                restingComp[j][0]+=restingComp[j][0]*factor1;
            }
            if (growthVal1<strainThreshold && growthVal2>strainThreshold){
	      cellData[cellIndex][growthInd]+=1;              
              double factor2=parameter(0)*h*(growthVal2-strainThreshold);
              
              for (size_t j=0; j< 3; ++j)
                restingComp[j][1]+=restingComp[j][1]*factor2;
            }
	    if (growthVal1>strainThreshold && growthVal2>strainThreshold){
              cellData[cellIndex][growthInd]+=2;
              double factor1=parameter(0)*h*(growthVal1-strainThreshold);
              double factor2=parameter(0)*h*(growthVal2-strainThreshold);
              
              for (size_t j=0; j< 3; ++j){
                restingComp[j][0]+=restingComp[j][0]*factor1;
                restingComp[j][1]+=restingComp[j][1]*factor2;
              }
            }
            
            double internalTemp=std::sqrt(restingComp[0][0]*restingComp[0][0]+
                                          restingComp[0][1]*restingComp[0][1]);
            
            double externalTemp=std::sqrt(restingComp[1][0]*restingComp[1][0]+
                                          restingComp[1][1]*restingComp[1][1]);
            
            double internalTempPlusOne=std::sqrt(restingComp[2][0]*restingComp[2][0]+
                                                 restingComp[2][1]*restingComp[2][1]);
	    
            cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex]
              =externalTemp-wallData[w2][wallLengthIndex];
            cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1]=internalTemp;
            cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod]=internalTempPlusOne;
          } // walls

          // updating wall length
          //} // if(cellData[cellIndex][36]==0)// hypocotyl 
        } // cells
        //std::cerr<<growthtime<<"  "<<totalAreatmp<<std::endl;
      }      
    }
  
    //reserve begin
    // {

    //   size_t dimension = vertexData[0].size();
    //   size_t numCells = T.numCell();
    //   size_t numWalls = T.numWall();
    //   size_t wallLengthIndex= variableIndex(0,0);
    //   size_t comIndex = variableIndex(1,0);
    //   size_t lengthInternalIndex = comIndex+dimension;
    //   size_t velocityStoreIndex = variableIndex(2,0);
    //   double strainThreshold=parameter(1);
    //   double velocityThreshold=parameter(2);
    //   size_t growthInd=30;

    //   static double growthtime=0;
    //   static double deltat=0;
    //   deltat +=h; 

    //   //std::cerr<<"deltat  growth           "<<deltat<<std::endl;
      
    //   bool equil=true;
    //   for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
    //     if(cellData[cellIndex][velocityStoreIndex]>velocityThreshold)
    //       equil=false;
      

    //   // if(!equil){
    //   //   for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
    //   //     std::cerr<<"cell Ind "<<cellIndex
    //   //              <<",   strain Thresh "<<strainThreshold 
    //   //              <<",  cellStrain "<<cellData[cellIndex][11]
    //   //              <<",  velocity "<<cellData[cellIndex][22]
    //   //              <<std::endl;
    //   // }
      
      
    //   if(equil) {// if close to mechanical equilibrum
      
        
    //     growthtime+=h;  
    //     //std::cerr<<scaletmp<<std::endl;
    //     std::vector<std::vector<double> > mainWalls(numWalls);
    //     std::vector<std::vector<std::vector<double> > > internalWalls(numCells);
        
    //     for (size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex)
    //       mainWalls[wallIndex].resize(2);
    //     for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){
    //       size_t numCellWalls = T.cell(cellIndex).numWall();
    //       internalWalls[cellIndex].resize(numCellWalls);
    //       for (size_t cellWallIndex=0 ; cellWallIndex<numCellWalls ; ++cellWallIndex)
    //         internalWalls[cellIndex][cellWallIndex].resize(2);
    //     }
        
    //     for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
    //       //if(cellData[cellIndex][37]==-1 || cellData[cellIndex][37]==-2 ||cellData[cellIndex][37]==-3 )// hypocotyl
    //       if(cellData[cellIndex][36]==0)// hypocotyl
    //         {
    //       if(numParameter()>4){//domain dependent growth threshold
    //         double conc=cellData[cellIndex][variableIndex(2,1)];
    //         double nHill=2;
    //         double hill=std::pow(conc,nHill)/(std::pow(parameter(4),nHill)+std::pow(conc,nHill));
    //         strainThreshold=(1-hill)*parameter(1)+hill*parameter(5);
    //       }
    //       cellData[cellIndex][growthInd]=0;
    //       size_t numCellWalls = T.cell(cellIndex).numWall(); 
    //       // internalWallLenth[cellIndex].resize(numCellWalls);
    //       // externalWallLength[cellIndex].resize(numCellWalls);
          
    //       // for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
    //       //   size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();         
    //       //   size_t restAddIndex=0;
    //       //   if(numParameter()==3){         
    //       //     if (T.wall(w2).cell1()->index()==cellIndex)
    //       //       restAddIndex=1;
    //       //     else if (T.wall(w2).cell2()->index()==cellIndex)
    //       //       restAddIndex=2;
    //       //     else{
    //       //       std::cerr<<"something bad happened"<<std::endl;
    //       //       exit(-1);
    //       //     }
    //       //   }
                    
    //       //   std::cerr<<wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex]<<"  ";
    //       // }
    //       // std::cerr<<"      ";
    //       // for (size_t wallIndex=lengthInternalIndex; wallIndex<lengthInternalIndex+8; ++wallIndex) { 
          
    //       //   std::cerr<<cellData[0][wallIndex]<<"      ";
    //       // }          
          
    //       // std::cerr<<" h "<<h<<std::endl;
                              
    //       for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
    //         size_t wallIndexPlusOneMod = (wallIndex+1)%numCellWalls;
    //         //size_t v1 = com;
    //         size_t v2 = T.cell(cellIndex).vertex(wallIndex)->index();
    //         size_t v3 = T.cell(cellIndex).vertex(wallIndexPlusOneMod)->index();
    //         //size_t w1 = internal wallIndex
    //         size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();
    //         //size_t w3 = internal wallIndex+1
            
    //         DataMatrix position(3,vertexData[v2]);
    //         for (size_t d=0; d<dimension; ++d)
    //           position[0][d] = cellData[cellIndex][comIndex+d]; // com position
    //         //position[1] = vertexData[v2]; // given by initiation
    //         position[2] = vertexData[v3];
            
    //         // size_t restAddIndex=0;
    //         // if(numParameter()==4 && parameter(3)==1){         
    //         //   if (T.wall(w2).cell1()->index()==cellIndex)
    //         //     restAddIndex=1;
    //         //   else if (T.wall(w2).cell2()->index()==cellIndex)
    //         //     restAddIndex=2;
    //         //   else{
    //         //     std::cerr<<"something bad happened"<<std::endl;
    //         //     exit(-1);
    //         //   }
    //         // }
                        
    //         std::vector<double> restingLength(3);
    //         if(numParameter()>3 && parameter(3)==1){ // double resting length        
    //           restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1];
    //           restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod];      
              
    //           restingLength[1] = wallData[w2][wallLengthIndex]+
    //             cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex];
    //         }
    //         else{           // single resting length
    //           restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallIndex];
    //           restingLength[2] = cellData[cellIndex][lengthInternalIndex + wallIndexPlusOneMod];      
    //           restingLength[1] = wallData[w2][wallLengthIndex];
    //         } 
                                    
    //         double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
    //                                       (-restingLength[0]+restingLength[1]+restingLength[2])*
    //                                       ( restingLength[0]-restingLength[1]+restingLength[2])*
    //                                       ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
            
    //         std::vector<double> length(3);
    //         length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
    //                                (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
    //                                (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
            
    //         length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
            
    //         length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
    //                                (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
    //                                (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
                                    
    //         //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
    //         double CurrentAngle1=std::acos(  (length[0]*length[0]+
    //                                           length[1]*length[1]-
    //                                           length[2]*length[2])/
    //                                          (length[0]*length[1]*2)    );
            
    //         double Qa=std::cos(CurrentAngle1)*length[0];
    //         double Qc=std::sin(CurrentAngle1)*length[0];
    //         double Qb=length[1];
                        
    //         double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+
    //                                           restingLength[1]*restingLength[1]-
    //                                           restingLength[2]*restingLength[2])/
    //                                          (restingLength[0]*restingLength[1]*2)    );
            
    //         double Pa=std::cos(RestingAngle1)*restingLength[0];
    //         double Pc=std::sin(RestingAngle1)*restingLength[0];
    //         double Pb=restingLength[1];
            
    //         // shape vector matrix in resting shape in local coordinate system  = 
    //         // inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
    //         double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
    //                                           {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
    //                                           { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
            
    //         double positionLocal[3][2]={ {Qa , Qc}, 
    //                                      {0  , 0 },  
    //                                      {Qb , 0 }  };
            
    //         double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
    //         for ( int ii=0 ; ii<3 ; ++ii ) {
    //           DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
    //           DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
    //           DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
    //           DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
    //         }
            
    //         double Egreen[2][2];//E=0.5(C-I)
    //         Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
    //         Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
    //         Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
    //         Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);
                        
    //         double LeftCauchy[2][2]; // B=FFt
    //         LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
    //         LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
    //         LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
    //         LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];
                        
    //         double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
    //         double tempS=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
    //         StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/tempS));
    //         StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/tempS;
    //         StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/tempS;  
    //         StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/tempS));
            
    //         double detTrue=StrainAlmansi[0][0]*StrainAlmansi[1][1]-StrainAlmansi[0][1]*StrainAlmansi[1][0];
    //         double trTrue=StrainAlmansi[0][0]+StrainAlmansi[1][1];
            
    //         double strainValue1True=trTrue/2+std::sqrt(std::fabs(((trTrue*trTrue)/4)-detTrue));
    //         double strainValue2True=trTrue/2-std::sqrt(std::fabs(((trTrue*trTrue)/4)-detTrue));
            
    //         //std::cerr<<strainValue1True<<"  "<<strainValue2True<<"      ";
            
    //         double det=Egreen[0][0]*Egreen[1][1]-Egreen[0][1]*Egreen[1][0];
    //         double tr=Egreen[0][0]+Egreen[1][1];

    //         double strainValue1=tr/2+std::sqrt(std::fabs(((tr*tr)/4)-det));
    //         double strainValue2=tr/2-std::sqrt(std::fabs(((tr*tr)/4)-det));
                        
    //         //std::cerr<<"strain values"<<strainValue1True<<"  "<<strainValue2True<<std::endl;
    //         double strainRestLocal1[2]={0,0};
    //         double strainRestLocal2[2]={0,0};
    //                                                           // is this the best way?
    //         if (std::fabs(Egreen[1][0])>0.0000000000000001){
    //           strainRestLocal1[0]=strainValue1-Egreen[1][1];
    //           strainRestLocal1[1]=Egreen[1][0];
    //           strainRestLocal2[0]=strainValue2-Egreen[1][1];
    //           strainRestLocal2[1]=Egreen[1][0];
    //         }
    //         else if (std::fabs(Egreen[0][1])>0.0000000000000001){
    //           strainRestLocal1[0]=Egreen[0][1];
    //           strainRestLocal1[1]=strainValue1-Egreen[0][0];
    //           strainRestLocal2[0]=Egreen[0][1];
    //           strainRestLocal2[1]=strainValue2-Egreen[0][0];
    //         }
    //         else {
    //           strainRestLocal1[0]=1;
    //           strainRestLocal1[1]=0;
    //           strainRestLocal2[0]=0;
    //           strainRestLocal2[1]=1;
    //         }
            
    //         double tempAn=std::sqrt(strainRestLocal1[0]*strainRestLocal1[0]+
    //                                 strainRestLocal1[1]*strainRestLocal1[1]);
    //         if (tempAn !=0 ){
    //           strainRestLocal1[0]/=tempAn;
    //           strainRestLocal1[1]/=tempAn;
    //         }
    //         else
    //           std::cerr<<"in growth.cc/strainTRBS strange strain eigenvector at cell "
    //                    <<cellIndex<<"  wall  "<<wallIndex<<std::endl;
            
    //         tempAn=std::sqrt(strainRestLocal2[0]*strainRestLocal2[0]+
    //                          strainRestLocal2[1]*strainRestLocal2[1]);
    //         if (tempAn !=0 ){
    //           strainRestLocal2[0]/=tempAn;
    //           strainRestLocal2[1]/=tempAn;
    //         }
    //         else
    //           std::cerr<<"in growth.cc/strainTRBS strange strain eigenvector at cell "
    //                    <<cellIndex<<"  wall  "<<wallIndex<<std::endl;
            
    //         if (strainValue2>strainValue1){
    //           double temp=strainValue1;
    //           strainValue1=strainValue2;
    //           strainValue2=temp;
    //           temp=strainRestLocal1[0];
    //           strainRestLocal1[0]=strainRestLocal2[0];
    //           strainRestLocal2[0]=temp;        
    //           temp=strainRestLocal1[1];
    //           strainRestLocal1[1]=strainRestLocal2[1];
    //           strainRestLocal2[1]=temp;          
    //         }
    //         if (strainValue2True>strainValue1True){
    //           double temp=strainValue1True;
    //           strainValue1True=strainValue2True;
    //           strainValue2True=temp;
    //         }
            
    //         std::vector<std::vector<double> > edgeRestLocal(3);
            
    //         for (size_t d=0; d< 3; ++d)
    //           edgeRestLocal[d].resize(3);
                        
    //         edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
    //         edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
    //         //edgeRestLocal[0][2]= positionLocal[][2]-positionLocal[][2];
            
    //         edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
    //         edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];
    //         //edgeRestLocal[1][2]= positionLocal[][2]-positionLocal[][2];
            
    //         edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
    //         edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
    //         //edgeRestLocal[2][2]= positionLocal[][2]-positionLocal[][2];
                        
    //         std:: vector<double> cosTet(3);
    //         std:: vector<double> sinTet(3);
    //         for (size_t j=0; j< 3; ++j){
    //           cosTet[j]=std::fabs((strainRestLocal1[0]*edgeRestLocal[j][0]+
    //                                strainRestLocal1[1]*edgeRestLocal[j][1])/
    //                               restingLength[j]);
    //           sinTet[j]=std::sqrt(std::fabs(1-cosTet[j]*cosTet[j]));
              
    //           if(cosTet[j]<0 || cosTet[j]>1)
    //             std::cerr<<"in growth.cc/strainTRBS cosTet is wrong: "
    //                      <<cosTet[j]<<" in cell "<<cellIndex<<std::endl;
    //         }
                        
    //         std::vector<std::vector<double> > restingComp(3);
    //         for (size_t j=0; j< 3; ++j)
    //           restingComp[j].resize(2);
            
    //         for (size_t j=0; j< 3; ++j){
    //           restingComp[j][0]=restingLength[j]*cosTet[j];
    //           restingComp[j][1]=restingLength[j]*sinTet[j];
    //         }
    //         // std::cerr<<"cell Ind "<<cellIndex
    //         //          <<",   strain1 "<<strainValue1True 
    //         //          <<",   strain2  "<<strainValue2True 
    //         //          <<",   strain Thresh "<<strainThreshold 
    //         //          <<",  cellStrain "<<cellData[cellIndex][11]
    //         //          <<",  velocity "<<cellData[cellIndex][22]
    //         //          <<std::endl;
            
    //         if (strainValue1True>strainThreshold && strainValue2True<strainThreshold){
            
    //           cellData[cellIndex][growthInd]+=1;
    //           //std::cerr<<"growth 1"<<std::endl;  
    //            // double factor1=parameter(0)*h*(strainValue1True-strainThreshold)/(1-2*strainValue1True);
    //            // factor1=factor1+factor1*factor1;// first and second terms in taylor expansion
              
    //           double factor1=parameter(0)*h*(strainValue1True-strainThreshold);
              
    //           for (size_t j=0; j< 3; ++j)
    //             restingComp[j][0]+=restingComp[j][0]*factor1;
              
    //         }
            
    //         if (strainValue1True>strainThreshold && strainValue2True>strainThreshold){
    //           cellData[cellIndex][growthInd]+=2;
    //           //std::cerr<<"growth 2"<<std::endl;  
    //            // double 
    //            //   factor1=parameter(0)*h*(strainValue1True-strainThreshold)/(1-2*strainValue1True);
    //            // double 
    //            //   factor2=parameter(0)*h*(strainValue2True-strainThreshold)/(1-2*strainValue2True);
              
    //            // factor1=factor1+factor1*factor1;// first and second terms in taylor expansion
    //            // factor2=factor2+factor2*factor2;// first and second terms in taylor expansion
              
    //            double factor1=parameter(0)*h*(strainValue1True-strainThreshold);
    //            double factor2=parameter(0)*h*(strainValue2True-strainThreshold);
              
              
    //           for (size_t j=0; j< 3; ++j){
    //             restingComp[j][0]+=restingComp[j][0]*factor1;
    //             restingComp[j][1]+=restingComp[j][1]*factor2;
    //           }
    //         }
            
    //         double internalTemp=std::sqrt(restingComp[0][0]*restingComp[0][0]+
    //                                       restingComp[0][1]*restingComp[0][1]);
            
    //         double externalTemp=std::sqrt(restingComp[1][0]*restingComp[1][0]+
    //                                       restingComp[1][1]*restingComp[1][1]);
            
    //         double internalTempPlusOne=std::sqrt(restingComp[2][0]*restingComp[2][0]+
    //                                              restingComp[2][1]*restingComp[2][1]);
                                    
    //         if(numParameter()>3 && parameter(3)==1){ // double length
    //           cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex]
    //             =externalTemp-wallData[w2][wallLengthIndex];
    //           cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1]=internalTemp;
    //           cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod]=internalTempPlusOne;
    //         }
    //         else{                  // single length
              
    //           // WITH AREA AVERAGING
              
    //           size_t wallGlobalInd= T.cell(cellIndex).wall(wallIndex) ->index();
    //           //std::cerr<<" normalization area factor before "<<mainWalls[wallGlobalInd][0]<<std::endl;
    //           if (mainWalls[wallGlobalInd][0]==0){
    //             mainWalls[wallGlobalInd][0]=restingArea;
    //             mainWalls[wallGlobalInd][1]=restingArea*externalTemp;
    //             //std::cerr<<" normalization area factor middle "<<mainWalls[wallGlobalInd][0]<<std::endl;
    //           }
    //           else if (mainWalls[wallGlobalInd][0]!=0){
    //             mainWalls[wallGlobalInd][0]+=restingArea;
    //             mainWalls[wallGlobalInd][1]=
    //               (mainWalls[wallGlobalInd][1]+restingArea*externalTemp);
    //           }
              
    //           if (internalWalls[cellIndex][wallIndex][0]==0){
    //             internalWalls[cellIndex][wallIndex][0]=restingArea;
    //             internalWalls[cellIndex][wallIndex][1]=restingArea*internalTemp;
    //           }
    //           else if (internalWalls[cellIndex][wallIndex][0]!=0){
    //             internalWalls[cellIndex][wallIndex][0]+=restingArea;
    //             internalWalls[cellIndex][wallIndex][1]=
    //               (internalWalls[cellIndex][wallIndex][1]+restingArea*internalTemp);
    //           }         
              
    //           if (internalWalls[cellIndex][wallIndexPlusOneMod][0]==0){
    //             internalWalls[cellIndex][wallIndexPlusOneMod][0]=restingArea;
    //             internalWalls[cellIndex][wallIndexPlusOneMod][1]=restingArea*internalTempPlusOne;
    //           }
    //           else if (internalWalls[cellIndex][wallIndexPlusOneMod][0]!=0){
    //             internalWalls[cellIndex][wallIndexPlusOneMod][0]+=restingArea;
    //             internalWalls[cellIndex][wallIndexPlusOneMod][1]=
    //               (internalWalls[cellIndex][wallIndexPlusOneMod][1]+restingArea*internalTempPlusOne);
                
    //           }         
              
    //         }
    //       } // walls
             
    //       //std::cerr<<std::endl;
    //       // updating wall length
          
    //     } // cells
                        
    //     if(numParameter()==3 ||(numParameter()>3 && parameter(3)!=1)){ // single resting length
    //       for (size_t cellIndex=0; cellIndex< T.numCell(); cellIndex++)
    //         for (size_t wallIndex=0; wallIndex< T.cell(cellIndex).numWall(); wallIndex++){
    //           cellData[cellIndex][lengthInternalIndex + wallIndex]= 
    //             internalWalls[cellIndex][wallIndex][1]
    //             /internalWalls[cellIndex][wallIndex][0]; 
    //           size_t wallGlobalInd=T.cell(cellIndex).wall(wallIndex)->index();
    //           wallData[wallGlobalInd][wallLengthIndex]=
    //             mainWalls[wallGlobalInd][1]/mainWalls[wallGlobalInd][0];
    //         }
    //     }

    //     double totalCellStrain=0;
    //     double totalCellStress=0;
    //     double totalAreaRatio=0;
    //     double totalWallStrain=0;
    //     double totalWallStress=0;
    //     double totalCellEnergy=0;
    //     double totalWallEnergy=0;
                
    //     double totalEdge=0;
                
    //     for(size_t a=0; a<numWalls; ++a)
    //       totalEdge+=wallData[a][0];
                
    //     for(size_t a=0; a<numCells; ++a){
    //       totalCellStrain +=cellData[a][11];  
    //       totalCellStress +=cellData[a][7];  
    //       totalAreaRatio  +=cellData[a][19];  
    //       totalCellEnergy +=cellData[a][20];  
    //     }
    //     for(size_t a=0; a<numWalls; ++a)
    //       for(size_t a=0; a<numWalls; ++a){
    //         size_t v1 = T.wall(a).vertex1()->index();
    //         size_t v2 = T.wall(a).vertex2()->index();
    //         size_t dimension = vertexData[v1].size();
    //         double distance=0;
    //         for( size_t d=0 ; d<dimension ; d++ )
    //           distance += (vertexData[v1][d]-vertexData[v2][d])*
    //             (vertexData[v1][d]-vertexData[v2][d]);
    //         distance = std::sqrt(distance);
    //         totalWallStrain +=0.5*(1-(wallData[a][1]/distance)); 
    //         totalWallStress +=wallData[a][2];
    //         totalWallEnergy +=wallData[a][1]*wallData[a][2]*wallData[a][2]*0.5/10;// frce^2/spring constant  
    //       }
                                 
    //     // double totalAreatmp=0;
    //     // for(size_t a=0; a<numCells; ++a)
    //     //   totalAreatmp+=cellData[a][19];
    //     // static double totalArea1=totalAreatmp; 

    //     //rescaling
    //     //double scaletmp=std::sqrt(totalArea1/totalAreatmp);//0.995;      
    //     // for(size_t a=0; a<numWalls; ++a)
    //     //   wallData[a][0]*=scaletmp;
    //     // for(size_t a=0; a<numCells; ++a)
    //     //   for(size_t b=42; b<cellData[a].size(); ++b)
    //     //     cellData[a][b]*=scaletmp;
    //     //    // rescaling       
    //     // for(size_t a=0; a<numWalls; ++a)
    //     //   wallData[a][1]*=scaletmp;

    //     //std::cout<<"cells  "<<" "<<scaletmp<<"  "
    //       //totalCellStress<<" "<<totalCellEnergy<<" "<<totalAreaRatio<<" "
    //     //          <<totalWallStrain<<" "<<totalWallStress<<" "<<totalWallEnergy<<" "
    //     //         <<totalAreatmp<<"  "<<totalArea1<<" "<<totalEdge<<std::endl;

    //     //deltat=0; 
    //     //std::cout<<growthtime<<"  "<<totalAreatmp<<"  "<<deltat<<std::endl;
    //   }
 
    // }
    // reserve end


    StrainTRBSConcentrationHill::
    StrainTRBSConcentrationHill(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {

      //Do some checks on the parameters and variable indeces
      if( paraValue.size()!=7) {
	std::cerr << "WallGrowthStrainTRBScenterTriangulationConcentrationHill::"
		  << "WallGrowthStrainTRBScenterTriangulationConcentrationHill() "
		  << "Uses three parameters k_growth, strain_threshold and velocity threshold"
                  << std::endl;
	exit(EXIT_FAILURE);
      }
        
      if( indValue.size() !=4 || 
          indValue[0].size() !=1 || 
          indValue[1].size() !=1 || 
          (indValue[2].size() !=1 && indValue[3].size() !=1) ) {
	std::cerr << "WallGrowthStrainTRBScenterTriangulation::"
		  << "WallGrowthStrainTRBScenterTriangulation() "
                  << "wall length index is given in first level,"
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in second level, " 
		  << "and velocity_value_strore_index and concentration "
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowthStrainTRBScenterTriangulationConcentrationHill");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //

      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "k_growth_hill";
      tmp[2] = "k_hill";
      tmp[3] = "n_hill";
      tmp[4] = "s_threshold";
      tmp[5] = "velocity_threshold";
      tmp[6] = "doubleFlag";
      
      setParameterId( tmp );
    }
    
    void StrainTRBSConcentrationHill::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
    }

    void StrainTRBSConcentrationHill::
    update(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
           double h )
    {
      size_t dimension = vertexData[0].size();
      size_t numCells = T.numCell();
      size_t numWalls = T.numWall();
      size_t wallLengthIndex= variableIndex(0,0);
      size_t comIndex = variableIndex(1,0);
      size_t lengthInternalIndex = comIndex+dimension;
      size_t velocityStoreIndex = variableIndex(2,0);
      size_t concentrationIndex = variableIndex(3,0);
      double k_hill = parameter(2);
      double n_hill = parameter(3);
      double strainThreshold=parameter(4);
      double velocityThreshold=parameter(5);
      size_t growthInd=30;

      static double growthtime=0;
      static double deltat=0;
      deltat +=h; 
            
      bool equil=true;
      for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){        
        if(cellData[cellIndex][velocityStoreIndex]>velocityThreshold)
          equil=false;
      }
      //std::cerr << equil << std::endl; //to check if it is ever updating
      if(equil && deltat>0.0) {// If close to mechanical equilibrum        
	//if(true) {  // If no check to be done
        growthtime+=h;  
        //std::cerr<<scaletmp<<std::endl;
        std::vector<std::vector<double> > mainWalls(numWalls);
        std::vector<std::vector<std::vector<double> > > internalWalls(numCells);
        
        for (size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex)
          mainWalls[wallIndex].resize(2);
        for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){
          size_t numCellWalls = T.cell(cellIndex).numWall();
          internalWalls[cellIndex].resize(numCellWalls);
          for (size_t cellWallIndex=0 ; cellWallIndex<numCellWalls ; ++cellWallIndex)
            internalWalls[cellIndex][cellWallIndex].resize(2);
        }
        
        for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
          //if(cellData[cellIndex][37]!=0) { // hypocotyl             
          
          size_t numCellWalls = T.cell(cellIndex).numWall(); 
          
          cellData[cellIndex][growthInd]=0;
          
          for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
            size_t wallIndexPlusOneMod = (wallIndex+1)%numCellWalls;
            //size_t v1 = com;
            size_t v2 = T.cell(cellIndex).vertex(wallIndex)->index();
            size_t v3 = T.cell(cellIndex).vertex(wallIndexPlusOneMod)->index();
            //size_t w1 = internal wallIndex
            size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();
            //size_t w3 = internal wallIndex+1
            
            DataMatrix position(3,vertexData[v2]);
            for (size_t d=0; d<dimension; ++d)
              position[0][d] = cellData[cellIndex][comIndex+d]; // com position
            //position[1] = vertexData[v2]; // given by initiation
            position[2] = vertexData[v3];
                        
            std::vector<double> restingLength(3);
            if(parameter(6)==1){ // double resting length        
              restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1];
              restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod];      
              
              restingLength[1] = wallData[w2][wallLengthIndex]+
                cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex];
            }
            else{           // single resting length
              restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallIndex];
              restingLength[2] = cellData[cellIndex][lengthInternalIndex + wallIndexPlusOneMod];      
              restingLength[1] = wallData[w2][wallLengthIndex];
            } 
            
            // double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
            //                               (-restingLength[0]+restingLength[1]+restingLength[2])*
            //                               ( restingLength[0]-restingLength[1]+restingLength[2])*
            //                               ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
	    std::vector<double> length(3);
            length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                                   (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                                   (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );            
            length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
	    length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                                   (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                                   (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
            
            //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
            double CurrentAngle1=std::acos(  (length[0]*length[0]+
                                              length[1]*length[1]-
                                              length[2]*length[2])/
                                             (length[0]*length[1]*2)    );
            
            double Qa=std::cos(CurrentAngle1)*length[0];
            double Qc=std::sin(CurrentAngle1)*length[0];
            double Qb=length[1];
                        
            double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+
                                              restingLength[1]*restingLength[1]-
                                              restingLength[2]*restingLength[2])/
                                             (restingLength[0]*restingLength[1]*2)    );
            
            double Pa=std::cos(RestingAngle1)*restingLength[0];
            double Pc=std::sin(RestingAngle1)*restingLength[0];
            double Pb=restingLength[1];
            
            // shape vector matrix in resting shape in local coordinate system  = 
            // inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
            double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                              {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                              { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
            
            double positionLocal[3][2]={ {Qa , Qc}, 
                                         {0  , 0 },  
                                         {Qb , 0 }  };
            
            double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
            for ( int ii=0 ; ii<3 ; ++ii ) {
              DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
              DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
              DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
              DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
            }
            
            size_t growthVecInd=0;
            size_t growthVal1Ind=7;
            size_t growthVal2Ind=11;

            // size_t growthVecInd=32;
            // size_t growthVal1Ind=35;
            // size_t growthVal2Ind=15;

            double growthVal1=cellData[cellIndex][growthVal1Ind];
            double growthVal2=cellData[cellIndex][growthVal2Ind];


            double growthCurrGlob[3]={cellData[cellIndex][growthVecInd  ],
                                      cellData[cellIndex][growthVecInd+1],
                                      cellData[cellIndex][growthVecInd+2] };
            
            
            double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                                   (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                                   (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );
            
            double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                                   (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                                   (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
            
      
            double Xcurrent[3]=      
              { (position[2][0]-position[1][0])/tempA,
                (position[2][1]-position[1][1])/tempA,
                (position[2][2]-position[1][2])/tempA
              };
            
            double Bcurrent[3]=      
              { (position[0][0]-position[1][0])/tempB,
                (position[0][1]-position[1][1])/tempB,
                (position[0][2]-position[1][2])/tempB
              };
            
            double Zcurrent[3]=      
              { Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1],
                Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2],
                Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0]
              };
            
            tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
            Zcurrent[0]=Zcurrent[0]/tempA;
            Zcurrent[1]=Zcurrent[1]/tempA;
            Zcurrent[2]=Zcurrent[2]/tempA;
            
            double Ycurrent[3]=      
              { Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1],
                Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2],
                Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0]
              };
            
            double rotation[3][3]=
              { {Xcurrent[0] , Ycurrent[0] , Zcurrent[0] },
                {Xcurrent[1] , Ycurrent[1] , Zcurrent[1] },
                {Xcurrent[2] , Ycurrent[2] , Zcurrent[2] } };
      
            // rotating the growth vector from global coordinate system to the local one in the current shape
            double growthCurrLocal[3]=
              {  rotation[0][0]*growthCurrGlob[0]+
                 rotation[1][0]*growthCurrGlob[1]+
                 rotation[2][0]*growthCurrGlob[2],
                 
                 rotation[0][1]*growthCurrGlob[0]+
                 rotation[1][1]*growthCurrGlob[1]+
                 rotation[2][1]*growthCurrGlob[2],
                 
                 rotation[0][2]*growthCurrGlob[0]+
                 rotation[1][2]*growthCurrGlob[1]+
                 rotation[2][2]*growthCurrGlob[2]
              };

            // transform the growth vector from current local to the resting local
            double growthRestLocal[2]={ DeformGrad[1][1] *growthCurrLocal[0]-DeformGrad[0][1]*growthCurrLocal[1],
                                        -DeformGrad[1][0]*growthCurrLocal[0]+DeformGrad[0][0]*growthCurrLocal[1]};
            double tmpG=std::sqrt(growthRestLocal[0]*growthRestLocal[0]+growthRestLocal[1]*growthRestLocal[1]);
            growthRestLocal[0]/=tmpG;
            growthRestLocal[1]/=tmpG;
            
            std::vector<std::vector<double> > edgeRestLocal(3);
	    for (size_t d=0; d< 3; ++d)
              edgeRestLocal[d].resize(3);
                        
            edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
            edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
	    edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
            edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];            
            edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
            edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
            
            std:: vector<double> cosTet(3);
            std:: vector<double> sinTet(3);
            for (size_t j=0; j< 3; ++j){
              cosTet[j]=std::fabs((growthRestLocal[0]*edgeRestLocal[j][0]+
                                   growthRestLocal[1]*edgeRestLocal[j][1])/
                                  restingLength[j]);
              sinTet[j]=std::sqrt(std::fabs(1-cosTet[j]*cosTet[j]));
              
              if(cosTet[j]<0 || cosTet[j]>1)
                std::cerr<<"in growth.cc/strainTRBS cosTet is wrong: "
                         <<cosTet[j]<<" in cell here   ..."<<cellIndex<<std::endl;
            }
            
            std::vector<std::vector<double> > restingComp(3);
            for (size_t j=0; j< 3; ++j)
              restingComp[j].resize(2);
            
            for (size_t j=0; j< 3; ++j){
              restingComp[j][0]=restingLength[j]*cosTet[j];
              restingComp[j][1]=restingLength[j]*sinTet[j];
            }
            
            // for hypocotyl Siobhan
            // if(cellIndex==408 && wallIndex==0){
            //     // std::cout<<strainThreshold<<" "<<growthVal1<<" "<<growthVal2<<" "
            //     //          <<growthtime<<" "<<deltat<<" "
            //     //          <<vertexData[0][2]-vertexData[899][2]<<std::endl;
            //     std::cout<<strainThreshold<<" "<<growthVal1<<" "<<growthVal2<<" "
            //              <<(growthVal1-growthVal2)/growthVal1<<std::endl;
                            
            // }
            
            if (growthVal1>strainThreshold && growthVal2<strainThreshold){
	      cellData[cellIndex][growthInd]+=1;
	      double conc = cellData[cellIndex][concentrationIndex];
	      double hillFactor = std::pow(conc, n_hill)/(std::pow(k_hill, n_hill)+std::pow(conc, n_hill));
              double factor1=(parameter(0)+parameter(1)*hillFactor)*h*(growthVal1-strainThreshold);
                                          
              for (size_t j=0; j< 3; ++j)
                restingComp[j][0]+=restingComp[j][0]*factor1;
            }
            if (growthVal1<strainThreshold && growthVal2>strainThreshold){
	      cellData[cellIndex][growthInd]+=1;
	      double conc = cellData[cellIndex][concentrationIndex];
	      double hillFactor = std::pow(conc, n_hill)/(std::pow(k_hill, n_hill)+std::pow(conc, n_hill));
	      double factor2=(parameter(0)+parameter(1)*hillFactor)*h*(growthVal2-strainThreshold);
              
              for (size_t j=0; j< 3; ++j)
                restingComp[j][1]+=restingComp[j][1]*factor2;
            }
	    if (growthVal1>strainThreshold && growthVal2>strainThreshold){
              cellData[cellIndex][growthInd]+=2;
	      double conc = cellData[cellIndex][concentrationIndex];
	      double hillFactor = std::pow(conc, n_hill)/(std::pow(k_hill, n_hill)+std::pow(conc, n_hill));
	      double factor1=(parameter(0)+parameter(1)*hillFactor)*h*(growthVal1-strainThreshold);
	      double factor2=(parameter(0)+parameter(1)*hillFactor)*h*(growthVal2-strainThreshold);

              for (size_t j=0; j< 3; ++j){
                restingComp[j][0]+=restingComp[j][0]*factor1;
                restingComp[j][1]+=restingComp[j][1]*factor2;
              }
            }
            
            double internalTemp=std::sqrt(restingComp[0][0]*restingComp[0][0]+
                                          restingComp[0][1]*restingComp[0][1]);
            
            double externalTemp=std::sqrt(restingComp[1][0]*restingComp[1][0]+
                                          restingComp[1][1]*restingComp[1][1]);
            
            double internalTempPlusOne=std::sqrt(restingComp[2][0]*restingComp[2][0]+
                                                 restingComp[2][1]*restingComp[2][1]);
	    
            cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex]
              =externalTemp-wallData[w2][wallLengthIndex];
            cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1]=internalTemp;
            cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod]=internalTempPlusOne;
          } // walls

          // updating wall length
          //} // if(cellData[cellIndex][36]==0)// hypocotyl 
        } // cells
        //std::cerr<<growthtime<<"  "<<totalAreatmp<<std::endl;
      }      
    }
  

    VectorTRBS::
    VectorTRBS(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) {
      
      // Do some checks on the parameters and variable indeces
      //
      if (paraValue.size()!=4) {
	std::cerr << "WallGrowth::CenterTriangulation::VectorTRBS"
		  << "VectorTRBS() "
		  << "Uses four parameters k_growth, s_threshold, velocity_threshold "
                  << "and double_edge_flag = one " 
                  << "for the case of using two independent resting lengths for neigbours, "
		  << "zero for using a single edge."
                  << std::endl;
	exit(EXIT_FAILURE);
      }
      if (paraValue[3]!=1) {
	std::cerr << "WallGrowth::CenterTriangulation::VectorTRBS"
		  << "VectorTRBS() "
		  << "Fourth parameter needs curfrently be one, "
                  << "i.e. double_edge_flag = one " 
                  << "for the case of using two independent resting lengths for neigbours."
                  << std::endl;
	exit(EXIT_FAILURE);
      }
        
      if( indValue.size() !=4 || 
          indValue[0].size() !=1 || 
          indValue[1].size() !=1 || 
          indValue[2].size() !=1 ||
	  indValue[3].size() !=3 ) {
	std::cerr << "WallGrowth::CenterTriangulation::VectorTRBS"
		  << "VectorTRBS() "
                  << "wall length index is given in first level,"
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in second level. " 
		  << "In third level velocity_value_store_index is given (not used"
		  << " if velocity_threshold parameter >=100."
		  << "In 4th level the two vector magnitude indices and vector"
		  << " start index is given."
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowth::CenterTriangulation::VectorTRBS");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "s_threshold";
      tmp[2] = "velocity_threshold";
      tmp[3] = "doubleFlag";      
      setParameterId( tmp );
    }
    
    void VectorTRBS::
      derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
    }

    void VectorTRBS::
    update(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
           double h )
    {
      size_t dimension = vertexData[0].size();
      size_t numCells = T.numCell();
      size_t numWalls = T.numWall();
      size_t wallLengthIndex= variableIndex(0,0);
      size_t comIndex = variableIndex(1,0);
      size_t lengthInternalIndex = comIndex+dimension;
      size_t velocityStoreIndex = variableIndex(2,0);
      double strainThreshold=parameter(1);
      double velocityThreshold=parameter(2);
      
      size_t growthVecInd=variableIndex(3,2);
      size_t growthVal1Ind=variableIndex(3,0);;
      size_t growthVal2Ind=variableIndex(3,1);;
      
      static double growthtime=0;
      static double deltat=0;
      deltat +=h; 
      
      bool equil=true;
      if (velocityThreshold<=100) {
	for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){        
	  if(cellData[cellIndex][velocityStoreIndex]>velocityThreshold)
	    equil=false;
	}
      }

      if(equil) {//close to equilibrium
	
        growthtime+=h;  
	// ----------------------------------------
	// Temporary elements for storing triangular growth contributions for
	// external and internal edges (only used if single edge elements are used)
        //std::vector<std::vector<double> > mainWalls(numWalls);
        //std::vector<std::vector<std::vector<double> > > internalWalls(numCells);        
        //for (size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex)
	//mainWalls[wallIndex].resize(2);
        //for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){
	//size_t numCellWalls = T.cell(cellIndex).numWall();
	//internalWalls[cellIndex].resize(numCellWalls);
	//for (size_t cellWallIndex=0 ; cellWallIndex<numCellWalls ; ++cellWallIndex)
	//  internalWalls[cellIndex][cellWallIndex].resize(2);
        //}
	// ----------------------------------------
        
        for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
	  double growthVal1=cellData[cellIndex][growthVal1Ind];
	  double growthVal2=cellData[cellIndex][growthVal2Ind];

	  if (growthVal1>strainThreshold || growthVal2>strainThreshold) {// Growth n at least one direction 
	    size_t numCellWalls = T.cell(cellIndex).numWall(); 
	    
	    for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
	      
	      // Extract the information for the vertices/edges to be used for this triangle/plate
	      // ----------------------------------------
	      size_t wallIndexPlusOneMod = (wallIndex+1)%numCellWalls;
	      //size_t v1 = com;
	      size_t v2 = T.cell(cellIndex).vertex(wallIndex)->index();
	      size_t v3 = T.cell(cellIndex).vertex(wallIndexPlusOneMod)->index();
	      //size_t w1 = internal edge
	      size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();
	      //size_t w3 = internal edge+1
	      
	      DataMatrix position(3,vertexData[v2]);
	      for (size_t d=0; d<dimension; ++d)
		position[0][d] = cellData[cellIndex][comIndex+d]; // com position
	      //position[1] = vertexData[v2]; // given by initiation
	      position[2] = vertexData[v3];
	      
	      std::vector<double> restingLength(3);
	      if(numParameter()>3 && parameter(3)==1){ // double resting length        
		restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1];
		restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod];      
		restingLength[1] = wallData[w2][wallLengthIndex]+
		  cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex];
	      }
	      else{           // single resting length
		std::cerr << "WallGRowth::CenterTriangulation::VectorTRBS "
			  << "Current implementation requires double resting lengths."
			  << std::endl;
		exit(EXIT_FAILURE);		
		// ----------------------------------------
		//restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallIndex];
		//restingLength[2] = cellData[cellIndex][lengthInternalIndex + wallIndexPlusOneMod];      
		//restingLength[1] = wallData[w2][wallLengthIndex];
		// ----------------------------------------
	      } 
            
	      std::vector<double> length(3);
	      length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
				     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
				     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );            
	      length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
	      length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
				     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
				     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
	      
	      //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
	      double CurrentAngle1=std::acos( (length[0]*length[0]+
					       length[1]*length[1]-
					       length[2]*length[2])/
					      (length[0]*length[1]*2) );
            
	      double Qa=std::cos(CurrentAngle1)*length[0];
	      double Qc=std::sin(CurrentAngle1)*length[0];
	      double Qb=length[1];
	      
	      double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+
						restingLength[1]*restingLength[1]-
						restingLength[2]*restingLength[2])/
					       (restingLength[0]*restingLength[1]*2)    );
            
	      double Pa=std::cos(RestingAngle1)*restingLength[0];
	      double Pc=std::sin(RestingAngle1)*restingLength[0];
	      double Pb=restingLength[1];
            
	      // shape vector matrix in resting shape in local coordinate system  = 
	      // inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
	      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
						{-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
						{ 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
            
	      double positionLocal[3][2]={ {Qa , Qc}, 
					   {0  , 0 },  
					   {Qb , 0 }  };
            
	      double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
	      for ( int ii=0 ; ii<3 ; ++ii ) {
		DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
		DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
		DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
		DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
	      }
            
	      double growthCurrGlob[3]={cellData[cellIndex][growthVecInd  ],
					cellData[cellIndex][growthVecInd+1],
					cellData[cellIndex][growthVecInd+2] };
	      
	      double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
				     (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
				     (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );
	      
	      double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
				     (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
				     (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
            
	      double Xcurrent[3]=      
		{ (position[2][0]-position[1][0])/tempA,
		  (position[2][1]-position[1][1])/tempA,
		  (position[2][2]-position[1][2])/tempA
		};
	      
	      double Bcurrent[3]=      
		{ (position[0][0]-position[1][0])/tempB,
		  (position[0][1]-position[1][1])/tempB,
		  (position[0][2]-position[1][2])/tempB
		};
	      
	      double Zcurrent[3]=      
		{ Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1],
		  Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2],
		  Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0]
		};
            
	      tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
	      Zcurrent[0]=Zcurrent[0]/tempA;
	      Zcurrent[1]=Zcurrent[1]/tempA;
	      Zcurrent[2]=Zcurrent[2]/tempA;
	      
	      double Ycurrent[3]=      
		{ Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1],
		  Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2],
		  Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0]
		};
	      
	      double rotation[3][3]=
		{ {Xcurrent[0] , Ycurrent[0] , Zcurrent[0] },
		  {Xcurrent[1] , Ycurrent[1] , Zcurrent[1] },
		  {Xcurrent[2] , Ycurrent[2] , Zcurrent[2] } };
	      
	      // rotating the growth vector from global coordinate system to the local one in the current shape
	      // ----------------------------------------
	      double growthCurrLocal[3]=
		{  rotation[0][0]*growthCurrGlob[0]+
		   rotation[1][0]*growthCurrGlob[1]+
		   rotation[2][0]*growthCurrGlob[2],
		   
		   rotation[0][1]*growthCurrGlob[0]+
		   rotation[1][1]*growthCurrGlob[1]+
		   rotation[2][1]*growthCurrGlob[2],
		   
		   rotation[0][2]*growthCurrGlob[0]+
		   rotation[1][2]*growthCurrGlob[1]+
		   rotation[2][2]*growthCurrGlob[2]
		};
	      
	      // transform the growth vector from current local to the resting local
	      double growthRestLocal[2]={DeformGrad[1][1]*growthCurrLocal[0]-DeformGrad[0][1]*growthCurrLocal[1],
					 -DeformGrad[1][0]*growthCurrLocal[0]+DeformGrad[0][0]*growthCurrLocal[1]};
	      double tmpG=std::sqrt(growthRestLocal[0]*growthRestLocal[0]+growthRestLocal[1]*growthRestLocal[1]);
	      growthRestLocal[0]/=tmpG;
	      growthRestLocal[1]/=tmpG;
	      
	      std::vector<std::vector<double> > edgeRestLocal(3);
	      for (size_t d=0; d< 3; ++d)
		edgeRestLocal[d].resize(3);
	      
	      edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
	      edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
	      edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
	      edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];            
	      edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
	      edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
	      
	      std:: vector<double> cosTet(3);
	      std:: vector<double> sinTet(3);
	      for (size_t j=0; j< 3; ++j){
		cosTet[j]=std::fabs((growthRestLocal[0]*edgeRestLocal[j][0]+
				     growthRestLocal[1]*edgeRestLocal[j][1])/
				    restingLength[j]);
		sinTet[j]=std::sqrt(std::fabs(1-cosTet[j]*cosTet[j]));
		
		if(cosTet[j]<0 || cosTet[j]>1) {
		  std::cerr<<"WallGrowth::CenterTriangulation::VectorTRBS cosTet is wrong, "
			   <<cosTet[j]<<", in cell " << cellIndex << " and wall/triangle "
			   << wallIndex << " and edge " << j << std::endl;
		  exit(EXIT_FAILURE);
		}
	      }
	      
	      std::vector<std::vector<double> > restingComp(3);
	      for (size_t j=0; j< 3; ++j)
		restingComp[j].resize(2);
	      
	      for (size_t j=0; j< 3; ++j){
		restingComp[j][0]=restingLength[j]*cosTet[j];
		restingComp[j][1]=restingLength[j]*sinTet[j];
	      }
	      
	      // Apply the growth update if applicable
	      if (growthVal1>strainThreshold && growthVal2<strainThreshold){
		//cellData[cellIndex][growthInd]+=1;              
		double factor1=parameter(0)*h*(growthVal1-strainThreshold);		
		for (size_t j=0; j< 3; ++j)
		  restingComp[j][0]+=restingComp[j][0]*factor1;
	      }
	      if (growthVal1<strainThreshold && growthVal2>strainThreshold){
		//cellData[cellIndex][growthInd]+=1;              
		double factor2=parameter(0)*h*(growthVal2-strainThreshold);		
		for (size_t j=0; j< 3; ++j)
		  restingComp[j][1]+=restingComp[j][1]*factor2;
	      }
	      if (growthVal1>strainThreshold && growthVal2>strainThreshold){
		//cellData[cellIndex][growthInd]+=2;
		double factor1=parameter(0)*h*(growthVal1-strainThreshold);
		double factor2=parameter(0)*h*(growthVal2-strainThreshold);		
		for (size_t j=0; j< 3; ++j){
		  restingComp[j][0]+=restingComp[j][0]*factor1;
		  restingComp[j][1]+=restingComp[j][1]*factor2;
		}
	      }
            
	      double internalTemp=std::sqrt(restingComp[0][0]*restingComp[0][0]+
					    restingComp[0][1]*restingComp[0][1]);
	      
	      double externalTemp=std::sqrt(restingComp[1][0]*restingComp[1][0]+
					    restingComp[1][1]*restingComp[1][1]);
	      
	      double internalTempPlusOne=std::sqrt(restingComp[2][0]*restingComp[2][0]+
						   restingComp[2][1]*restingComp[2][1]);
	      
	      cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex]
		=externalTemp-wallData[w2][wallLengthIndex];
	      cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1]=internalTemp;
	      cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod]=internalTempPlusOne;
	    } // walls
	  } //if (growthVal1>strainThreshold || growthVal2>strainThreshold)
        } // cells
      } // if (equil)      
    }// end update()

  } //CenterTriangulation namespace
  
  StressSpatial::
  StressSpatial(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=6 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "Uses six parameters k_growth, stress(stretch)_threshold "
		<< "K_hill n_Hill "
		<< "stretch_flag and linear_flag" << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[4] != 0.0 && paraValue[4] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "Two variable index is used (wall length,spatial coordinate) at first "
		<< "level, and force variable index at second."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressSpatial");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    Kpow_=std::pow(paraValue[2],paraValue[3]);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "stretch_flag";
    tmp[5] = "linear_flag";
    setParameterId( tmp );
  }
  
  void StressSpatial::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    
    // Prepare spatial factor
    size_t sI=variableIndex(0,1);
    assert (sI<vertexData[0].size());
    size_t numVertices = vertexData.size();
    double sMax= vertexData[0][sI];
    size_t maxI=0;
    for (size_t i=1; i<numVertices; ++i)
      if (vertexData[i][sI]>sMax) {
	sMax=vertexData[i][sI];
	maxI=i;
      }
    std::vector<double> maxPos(dimension);
    for (size_t d=0; d<dimension; ++d)
      maxPos[d] = vertexData[maxI][d];
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(4)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	// Calculate spatial factor
	double maxDistance = 0.0;
	for (size_t d=0; d<dimension; ++d) {
	  double pos = 0.5*(vertexData[v1][d]+vertexData[v2][d]);
	  maxDistance += (maxPos[d]-pos)*(maxPos[d]-pos);
	}
	maxDistance = std::sqrt(maxDistance);
	double spatialFactor = Kpow_/(Kpow_+std::pow(maxDistance,parameter(3)));
	
	double growthRate = parameter(0)*(stress - parameter(1))*spatialFactor;
		
	if (parameter(5))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  StressSpatialSingle::
  StressSpatialSingle(std::vector<double> &paraValue, 
		      std::vector< std::vector<size_t> > 
		      &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=6 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "Uses six parameters k_growth, stress(stretch)_threshold "
		<< "K_hill n_Hill "
		<< "stretch_flag and linear_flag" << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[4] != 0.0 && paraValue[4] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "Two variable index is used (wall length,spatial coordinate) at first "
		<< "level, and force variable index at second."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressSpatialSingle");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    Kpow_=std::pow(paraValue[2],paraValue[3]);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "stretch_flag";
    tmp[5] = "linear_flag";
    setParameterId( tmp );
  }
  
  void StressSpatialSingle::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    
    // Prepare spatial factor
    size_t sI=variableIndex(0,1);
    assert (sI<vertexData[0].size());
    size_t numVertices = vertexData.size();
    double sMax= vertexData[0][sI];
    size_t maxI=0;
    for (size_t i=1; i<numVertices; ++i)
      if (vertexData[i][sI]>sMax) {
	sMax=vertexData[i][sI];
	maxI=i;
      }
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(4)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	// Calculate spatial factor
	double maxDistance = sMax - 0.5*(vertexData[v1][sI]+vertexData[v2][sI]);;
	double spatialFactor = Kpow_/(Kpow_+std::pow(maxDistance,parameter(3)));
	
	double growthRate = parameter(0)*(stress - parameter(1))*spatialFactor;
	
	if (parameter(5))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  StressConcentrationHill::
  StressConcentrationHill(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=7 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "Uses seven parameters k_growthConst, k_growthHill, K_Hill, n_Hill,"
		<< " stretch_threshold stretch_flag and linear_flag" << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch/strain used)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[6] != 0.0 && paraValue[6] != 1.0 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth (exponential))." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "wall length index and concentration index at first "
		<< "level, and spring constant variable indices at second"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressConcentrationHill");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "stretch_flag";
    tmp[3] = "linear_flag";
    setParameterId( tmp );
  }

  void StressConcentrationHill::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t concIndex = variableIndex(0,1);
    double Kpow = std::pow(parameter(2),parameter(3));
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(5)) {
	     for (size_t k=0; k<numVariableIndex(1); ++k)
	       stress += wallData[i][variableIndex(1,k)];
      }
      else {
      	double distance=0.0;
      	for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
      	  distance += (vertexData[v1][d]-vertexData[v2][d])*
      	    (vertexData[v1][d]-vertexData[v2][d]);
      	  distance = std::sqrt(distance);
      	  stress = (distance-wallData[i][lengthIndex]) /
      	  wallData[i][lengthIndex];
      }
      if (stress > parameter(4)) {
      	// Get the Hill factor from the two cells
      	double hillFactor=0.0;
      	if (T.wall(i).cell1() != T.background()) {
      	  double concpow = std::pow(cellData[T.wall(i).cell1()->index()][concIndex],parameter(3));
      	  hillFactor += concpow/(Kpow+concpow); 
      	}
      	if (T.wall(i).cell2() != T.background()) {
      	  double concpow = std::pow(cellData[T.wall(i).cell2()->index()][concIndex],parameter(3));
      	  hillFactor += concpow/(Kpow+concpow); 
      	}
      	double growthRate = (parameter(0)+hillFactor*parameter(1))*(stress - parameter(4));
      	if (parameter(6))
      	  growthRate *= wallData[i][lengthIndex];
      	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  ConstantStressEpidermalAsymmetric::
  ConstantStressEpidermalAsymmetric(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 ) {
      std::cerr << "WallGrowth::ConstantStressEpidermalAsymmetric::"
		<< "ConstantStressEpidermalAsymmetric() "
		<< "Uses two parameters k_growth and frac_epi\n";
      exit(EXIT_FAILURE);
    }  
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "WallGrowth::ConstantStressEpidermalAsymmetric::"
		<< "ConstantStressEpidermalAsymmetric() "
		<< "One variable index is used.\n";
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::ConstantStressEpidermalAsymmetric");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[0] = "frac_epi";
    setParameterId( tmp );
  }
  
  void ConstantStressEpidermalAsymmetric::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double kGrowth = parameter(0);
      if( T.wall(i).cell1() == T.background() || 
	  T.wall(i).cell2() == T.background() )
	kGrowth *= parameter(1);
      double distance=0.0;
      for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	distance += (vertexData[v1][d]-vertexData[v2][d])*
	  (vertexData[v1][d]-vertexData[v2][d]);
      distance = std::sqrt(distance);
      if( distance>wallData[i][lengthIndex] )
	wallDerivs[i][lengthIndex] += kGrowth*
	  (distance-wallData[i][lengthIndex]);
    }
  }
  
  Force::Force(std::vector<double> &paraValue,
	       std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 2) {
      std::cerr << "WallGrowth::Force::Force() "
		<< "Uses two parameters: k_growth and Force_threshold" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "WallGrowth::Force::Force() "
		<< "Wall length index must be given in first level.\n"
		<< "Wall force index/indices must be given in second level.\n";
      exit(EXIT_FAILURE);
    }
    
    setId("VertexFromWallSpringExperimental");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_L";
    tmp[1] = "phi";
    
    setParameterId(tmp);
  }

  void Force::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
  {
    for (size_t i = 0; i < T.numWall(); ++i) {
      size_t vertex1Index = T.wall(i).vertex1()->index();
      size_t vertex2Index = T.wall(i).vertex2()->index();
      size_t dimensions = vertexData[vertex1Index].size();
      
      double distance = 0.0;
      for (size_t d = 0; d < dimensions; ++d) {
	distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
	  * (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
      }
      distance = std::sqrt(distance);
      
      double F = 0.0;
      for (size_t j = 0; j < numVariableIndex(1); ++j)
	F += wallData[T.wall(i).index()][variableIndex(1, j)];
      
      double arg = F - parameter(1);
      if (arg > 0)
	wallDerivs[T.wall(i).index()][variableIndex(0, 0)] 
	  += parameter(0) * arg; 
      //*wallData[T.wall(i).index()][variableIndex(0, 0)];
    }
  }
} // end namespace WallGrowth
