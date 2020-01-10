//
// Filename     : fiberModel.cc
// Description  : Classes describing updates of mechanical properties, e.g. Young's modulii given input
//                from elsewhere
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//

#include "baseReaction.h"
#include "fiberModel.h"
#include "tissue.h"

namespace FiberModel {

  General::General(std::vector<double> &paraValue,
			 std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 8 ) {
      std::cerr << "FiberModel::General::General() Uses eight parameters:" << std::endl
		<< "k_rate, equilibrium_threshold (only update if vertex moving slower), " << std::endl
		<< "linear-hill flag (0=linear,1=Hill,2=Hill_direct)," << std::endl
		<< "K_hill, n_hill (only used if Hill versions selected)," << std::endl
		<< "young_matrix, young_fiber, (material limits, should be same as for the mechanical (TRBS) model used)"
		<< std::endl << "initialization flag (0=no initiation/1=isotropic/2=anisotropic" << std::endl;
      exit(EXIT_FAILURE);
    }
  
    if (indValue.size() != 3 || indValue[0].size() != 1 ||
	(indValue[1].size() != 1 && indValue[1].size() != 2) ||
	indValue[2].size() != 1) {
      std::cerr << "FiberModel::General::General() " << std::endl
		<< "First level gives stress/strain anisotropy index to read." << std::endl
		<< "Second level gives Young_Longitudinal index to update (p_2=1,2) or Y_L and F_L (p_2=2) "
		<< "where the first one is read and the second updated." << std::endl
		<< "Third level gives index for velocity to be compared with p_1." << std::endl;
      exit(EXIT_FAILURE);
    }
  
    setId("FiberModel::General");
    setParameter(paraValue);
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_rate";
    tmp[1] = "velocitythreshold";
    tmp[2] = "linear-hill-flag";
    tmp[3] = "k_hill";
    tmp[4] = "n_hill";
    tmp[5] = "Y_matrix";
    tmp[6] = "Y_fiber";
    tmp[7] = "init_flag";
    //if (parameter(7)==3)
    //  tmp[8] = "Poisson";

    setParameterId(tmp);
  }

  void General::initiate(Tissue &T,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs) {
    size_t numCell=cellData.size();
    size_t AnisoIndex=variableIndex(0,0);
    size_t YoungLIndex=variableIndex(1,0);
    double Kh=parameter(3);
    double Nh=parameter(4);
    double youngMatrix=parameter(5);
    double youngFiber=parameter(6);

    if (parameter(7)==1){ // initialize with iso material
      for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { // initiating with 0 anisotropy and isotropic material
	cellData[cellIndex][AnisoIndex] = 0;
	cellData[cellIndex][YoungLIndex] = youngMatrix+0.5*youngFiber; // youngL = youngMatrix + 0.5*youngFiber;
	//std::cerr<< cellData[cellIndex][variableIndex(1,0)] << std::endl;
      }
    }
    if (parameter(7)==2){
      for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { // initiating with current anisotropy 
	double anisotropy=cellData[cellIndex][AnisoIndex];
	if ( parameter(2)==0 ) // linear
	  cellData[cellIndex][YoungLIndex]=youngMatrix+0.5*(1+anisotropy)* youngFiber;
      
	if ( parameter(2)==1 || parameter(2)==2 ) // Hill
	  cellData[cellIndex][YoungLIndex] =  youngMatrix+ 
	    0.5*(1+(std::pow(anisotropy,Nh) /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
					      +std::pow(anisotropy,Nh))))* youngFiber;
	//std::cerr<< cellData[cellIndex][variableIndex(1,0)] << std::endl;
      }
    }
  }

  void General::derivs(Tissue &T,
		       DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       DataMatrix &cellDerivs,
		       DataMatrix &wallDerivs,
		       DataMatrix &vertexDerivs ) {} //

  void General::update(Tissue &T,
		       DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData, 
		       double h) 
  {
    size_t numCell=cellData.size();
    size_t AnisoIndex=variableIndex(0,0);
    size_t YoungLIndex=variableIndex(1,0);
    size_t velocityIndex=variableIndex(2,0);
    double Kh=parameter(3);
    double Nh=parameter(4);
    double youngMatrix=parameter(5);
    double youngFiber=parameter(6);
    //double poisson=parameter(7);

    if (parameter(0)==0.0)
      return;

    for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
      double anisotropy=cellData[cellIndex][AnisoIndex];

      // --- Linear feedback rule --- 
      if ( parameter(2)==0 // linear
	   && cellData[cellIndex][velocityIndex] < parameter(1) 
	   //&& cellData[cellIndex][YoungLIndex] < youngMatrix+0.5*(1+anisotropy)* youngFiber
	   && cellData[cellIndex][YoungLIndex] < youngMatrix+youngFiber // not to exceed maximum when using absolute stress anisotropy 
	   ){ // linear Fiber model
	double deltaYoungL=youngMatrix+0.5*(1+anisotropy)* youngFiber
	  -cellData[cellIndex][YoungLIndex];
	cellData[cellIndex][YoungLIndex] += parameter(0)*h*deltaYoungL;
	//cellData[cellIndex][YoungLIndex] = youngMatrix+0.5*(1+anisotropy)* youngFiber;
      }

      if ( parameter(2)==1 // Hill
	   && cellData[cellIndex][velocityIndex] < parameter(1) 
	   // && cellData[cellIndex][YoungLIndex] < 
	   // youngMatrix+
	   // 0.5*(1+(std::pow(anisotropy,Nh)
	   //         /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))* youngFiber
	   && cellData[cellIndex][YoungLIndex] < youngMatrix+youngFiber // not to exceed maximum when using absolute stress anisotropy 
	   ){ // Hill-like Fiber model
      
	double deltaYoungL=youngMatrix+
	  0.5*(1+(std::pow(anisotropy,Nh)
		  /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))* youngFiber
	  -cellData[cellIndex][YoungLIndex];
	cellData[cellIndex][YoungLIndex]+= parameter(0)*h*deltaYoungL;
	// cellData[cellIndex][YoungLIndex] = youngMatrix+
	//   0.5*(1+(std::pow(anisotropy,Nh)/(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))* youngFiber; 
      }
      if ( parameter(2)==2 ){ // Hill-like Fiber model instantenous update for fiber deposition model
	size_t fiberLIndex=variableIndex(1,1);
	youngFiber= cellData[cellIndex][YoungLIndex];
      
	cellData[cellIndex][fiberLIndex] =
	  0.5*(1+(std::pow(anisotropy,Nh)
		  /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))* youngFiber;      
      }
    }
  }

  Deposition::Deposition(std::vector<double> &paraValue,
			 std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 7) {
      std::cerr << "FiberModel::Deposition::Deposition() " 
		<< "Uses seven parameters: k_rate, initial randomness(percent),"
		<< "velocity threshold, initiation flag(0/1), k_hill and n_hill" 
		<< "fiber strength" 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
  
    if (indValue.size() != 1 || indValue[0].size() != 6) {
      std::cerr << "FiberModel::Deposition::Deposition() " 
		<< "First index: anisotropy." 
		<< "second index:  misses stress." 
		<< "third index: cell area."
		<< "fourth index: fiber." 
		<< "fifth index: longitudinal young fiber"
		<< "fourth index: velocity"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
  
    setId("FiberModel::Deposition");
    setParameter(paraValue);
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_rate";
    tmp[1] = "randomness";
    tmp[2] = "velocity_threshold";
    tmp[3] = "init_flag";
    tmp[4] = "k_Hill";
    tmp[5] = "n_Hill";
    tmp[6] = "fiber";
  
    setParameterId(tmp);
  }

  void Deposition::initiate(Tissue &T,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    DataMatrix &cellDerivs,
			    DataMatrix &wallDerivs,
			    DataMatrix &vertexDerivs) {
    size_t numCells=cellData.size();
    size_t YoungFiberIndex=variableIndex(0,3);
    size_t YoungFiberLIndex=variableIndex(0,4);
  
    double randCoef=parameter(1);
    double youngFiber=parameter(6);
    size_t initFlag=parameter(3);  

    if (initFlag==1) // initialize with uniform material
      for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){ 
	cellData[cellIndex][YoungFiberIndex] = youngFiber; 
	cellData[cellIndex][YoungFiberLIndex] = youngFiber/2; 
      }
    if (initFlag==2 || initFlag==3){ // initialize with random material stiffness
      srand((unsigned)time(NULL));
      for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){ 
	cellData[cellIndex][YoungFiberIndex] = youngFiber*(1+randCoef*2*(0.5-((double) rand()/(RAND_MAX)))); 
	cellData[cellIndex][YoungFiberLIndex] = cellData[cellIndex][YoungFiberIndex]/2; 

	std::cerr<<cellData[cellIndex][YoungFiberIndex]<<"  "<<cellData[cellIndex][YoungFiberLIndex]<<std::endl;
      }
    }
    return;  
  }

  void Deposition::derivs(Tissue &T,
			  DataMatrix &cellData,
			  DataMatrix &wallData,
			  DataMatrix &vertexData,
			  DataMatrix &cellDerivs,
			  DataMatrix &wallDerivs,
			  DataMatrix &vertexDerivs ) {}
  

  void Deposition::update(Tissue &T,
			  DataMatrix &cellData,
			  DataMatrix &wallData,
			  DataMatrix &vertexData, 
			  double h) 
  { 
    size_t numCells=cellData.size();
    
    size_t AnisoIndex=variableIndex(0,0);
    size_t MstressIndex=variableIndex(0,1);
    size_t areaIndex=variableIndex(0,2);
    size_t YoungFiberIndex=variableIndex(0,3);
    size_t YoungFiberLIndex=variableIndex(0,4);
    size_t velocityIndex=variableIndex(0,5);

    double k_rate=parameter(0);
    //double randCoef=parameter(1);
    double vThresh=parameter(2);
    double Kh=parameter(4);
    double Nh=parameter(5);
    double YoungFiber=parameter(6);
    if (parameter(0)==0.0)
      return;
    
    double smax=14;//28;
    double smin=8;//;
    
    bool equil=true;
    for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
      if(cellData[cellIndex][velocityIndex]>vThresh)
	equil=false;
    
    if(equil){
      double totalStressArea=0;
      double totalArea=0;
      
      for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
	if(cellData[cellIndex][MstressIndex]<smax && cellData[cellIndex][MstressIndex]>smin){
	  totalStressArea+=(cellData[cellIndex][MstressIndex]-smin)
	    *cellData[cellIndex][areaIndex];
	}
	else if(cellData[cellIndex][MstressIndex]>smax){
	  totalStressArea+=(smax-smin)
	    *cellData[cellIndex][areaIndex];
	}
	totalArea+=cellData[cellIndex][areaIndex];
      }
      
      for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
	double targetFiber=0;
	if(cellData[cellIndex][MstressIndex]<smax && cellData[cellIndex][MstressIndex]>smin){
	  targetFiber=(YoungFiber*totalArea*(cellData[cellIndex][MstressIndex]-smin))
	    /totalStressArea;      
	}
	else if(cellData[cellIndex][MstressIndex]>smax ){
	  targetFiber=(YoungFiber*totalArea*(smax-smin))
	    /totalStressArea;      
	}
	else if(cellData[cellIndex][MstressIndex]<smin ){
	  targetFiber=0;      
	}
	cellData[cellIndex][YoungFiberIndex]+=k_rate*
	  (targetFiber
	   -cellData[cellIndex][YoungFiberIndex]);
	
	if (parameter(3)==3){      
	  double anisotropy= cellData[cellIndex][AnisoIndex];
	  cellData[cellIndex][YoungFiberLIndex]+=k_rate*
	    ((0.5*(1+(std::pow(anisotropy,Nh)
		      /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
			+std::pow(anisotropy,Nh))))* cellData[cellIndex][YoungFiberIndex])
	     -cellData[cellIndex][YoungFiberLIndex]);
	}
	else
	  cellData[cellIndex][YoungFiberLIndex]=cellData[cellIndex][YoungFiberIndex]/2;
      }
      // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
      //   totalStressArea+=cellData[cellIndex][MstressIndex]
      //     *cellData[cellIndex][areaIndex];
      //   totalArea+=cellData[cellIndex][areaIndex];
      // }
      // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
      //   cellData[cellIndex][YoungFiberIndex]+=k_rate*
      //     ((YoungFiber*totalArea*cellData[cellIndex][MstressIndex])
      //      /totalStressArea
      //      -cellData[cellIndex][YoungFiberIndex]);
      //   cellData[cellIndex][YoungFiberLIndex]=cellData[cellIndex][YoungFiberIndex]/2;
      //   // double anisotropy= cellData[cellIndex][AnisoIndex];
      //   // cellData[cellIndex][YoungFiberLIndex]+=k_rate*
      //   //   ((0.5*(1+(std::pow(anisotropy,Nh)
      //   //           /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
      //   //             +std::pow(anisotropy,Nh))))* cellData[cellIndex][YoungFiberIndex])
      //   //   -cellData[cellIndex][YoungFiberLIndex]);
      // }
      // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
      //   totalStressArea+=cellData[cellIndex][MstressIndex]*cellData[cellIndex][areaIndex];
      //   totalArea+=cellData[cellIndex][areaIndex];
      // }
      // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
      //   cellData[cellIndex][YoungFiberIndex]+=k_rate*
      //     ((YoungFiber*totalArea*cellData[cellIndex][MstressIndex])/totalStressArea
      //      -cellData[cellIndex][YoungFiberIndex]);
      //   cellData[cellIndex][YoungFiberLIndex]=cellData[cellIndex][YoungFiberIndex]/2;
      //   // double anisotropy= cellData[cellIndex][AnisoIndex];
      //   // cellData[cellIndex][YoungFiberLIndex]+=k_rate*
      //   //   ((0.5*(1+(std::pow(anisotropy,Nh)
      //   //           /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
      //   //             +std::pow(anisotropy,Nh))))* cellData[cellIndex][YoungFiberIndex])
      //   //   -cellData[cellIndex][YoungFiberLIndex]);
      // }  
    }
    // double tmp;
    // //std::cerr<<"total    "<<totalStressArea<<std::endl;
  
    // if(totalStressArea !=0){
    //   for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex){ 
    //     tmp=h*k_rate*cellData[cellIndex][YoungFiberIndex];
    //     cellData[cellIndex][YoungFiberIndex]-=tmp;
    //     totalFiber+=tmp*cellData[cellIndex][areaIndex];
    //   }
    
    
    //   for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
    //     double anisotropy= cellData[cellIndex][AnisoIndex];
    //     cellData[cellIndex][YoungFiberIndex]+=(totalFiber*cellData[cellIndex][MstressIndex])/totalStressArea;
    //     cellData[cellIndex][YoungFiberLIndex] =
    //       0.5*(1+(std::pow(anisotropy,Nh)
    //               /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
    //                 +std::pow(anisotropy,Nh))))* cellData[cellIndex][YoungFiberIndex];
    //   }
    // }

  }
} // namespace FiberModel
