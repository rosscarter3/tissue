// Filename     : hypocotyl3D.cc
// Description  : Classes describing special reactions for a 3D hypocotyl growth simulation
// Author(s)    : Behruz Bozorg
//              : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : August 2018
// Revision     : $Id:$
//
#include <fstream>
#include <iostream>
#include <sstream>

#include "hypocotyl3D.h"
#include "baseReaction.h"

namespace Hypocotyl3D {

  limitZdis::limitZdis(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 0) {
      std::cerr << "Hypocotyl3D::limitZdis::limitZdis() "
		<< "Uses no parameter, distributes the forces at the "
		<< "two ends of a cylinder correctly for treatment of boundaries."
		<< " Uses hard coded collection of vertex indices "
		<< "and geometries defined "
		<< " only for hypocotyl paper (Bou Daher et al (2018)."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (indValue.size() != 1 || indValue[0].size() != 2) {
      std::cerr << "Hypocotyl3D::limitZdis::limitZdis() "
		<< "one index level with 2 indices."
		<< " Uses hard coded collection of vertex indices "
		<< "and geometries defined "
		<< " only for hypocotyl paper (Bou Daher et al (2018)."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    
    setId("Hypocotyl3D::limitZdis");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
  
  void limitZdis::initiate(Tissue &T,
			   DataMatrix &cellData,
			   DataMatrix &wallData,
			   DataMatrix &vertexData,
			   DataMatrix &cellDerivs,
			   DataMatrix &wallDerivs,
			   DataMatrix &vertexDerivs)
  {
    size_t numCells = T.numCell();
    for (size_t cellIndex= 0; cellIndex< numCells; ++cellIndex)
      if(cellData[cellIndex][38]==-1 && (cellData[cellIndex][37]==-4 ||cellData[cellIndex][37]==0) ) // for hypocotyl
	{ size_t numWalls = T.cell(cellIndex).numWall();      
	  for (size_t wallindex=0; wallindex<numWalls; ++wallindex) { 
	    size_t vInd= T.cell(cellIndex).vertex(wallindex)->index();
	    if(vertexData[vInd][2]>-50)
	      topVertices.push_back(vInd);
	    if(vertexData[vInd][2]<-50)
	      bottomVertices.push_back(vInd);
	  }
	}
  }

  void limitZdis::derivs(Tissue &T,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs) 
  {
    size_t numTop=topVertices.size();
    size_t numBottom=bottomVertices.size();
    
    double tmpZforce=0;
    for (size_t i= 0; i< numTop; ++i)
      tmpZforce+=vertexDerivs[topVertices[i]][2];
    tmpZforce/=numTop;
    for (size_t i= 0; i< numTop; ++i)
      vertexDerivs[topVertices[i]][2]=tmpZforce;
    
    tmpZforce=0;
    for (size_t i= 0; i< numBottom; ++i)
      tmpZforce+=vertexDerivs[bottomVertices[i]][2];
    tmpZforce/=numBottom;
    for (size_t i= 0; i< numBottom; ++i)
      vertexDerivs[bottomVertices[i]][2]=tmpZforce;
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
		<< "or four parameters with the forth one equal to one " 
		<< "for the case of using two independent resting lengths for neigbohrs. "
		<< "or six parameters with the fifth one the K value for a concentration "
		<<" hill function with n=2 and the sixth parameter v_max"
		<< std::endl;
      exit(0);
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
		<< "and velocity_value_strore_index and optional concentration "
		<< "index for domain dependent growth threshold in third level."
		<< std::endl;
      exit(0);
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
    
  //reserve begin
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
    
    //std::cerr<<"deltat  growth           "<<deltat<<std::endl;
    
    bool equil=true;
    for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
      if(cellData[cellIndex][velocityStoreIndex]>velocityThreshold)
	equil=false;
    
    
    // if(!equil){
    //   for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
    //     std::cerr<<"cell Ind "<<cellIndex
    //              <<",   strain Thresh "<<strainThreshold 
    //              <<",  cellStrain "<<cellData[cellIndex][11]
    //              <<",  velocity "<<cellData[cellIndex][22]
    //              <<std::endl;
    // }
    
    
    if(equil) {// if close to mechanical equilibrum
      
      
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
      
      for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) 
	
	if(parameter(3)!=-1 || (parameter(3)==-1 && cellData[cellIndex][36]==0))// hypocotyl boundary should not grow
	  {
	    if(numParameter()>4){//domain dependent growth threshold
	      double conc=cellData[cellIndex][variableIndex(2,1)];
	      double nHill=2;
	      double hill=std::pow(conc,nHill)/(std::pow(parameter(4),nHill)+std::pow(conc,nHill));
	      strainThreshold=(1-hill)*parameter(1)+hill*parameter(5);
	    }
	    cellData[cellIndex][growthInd]=0;
	    size_t numCellWalls = T.cell(cellIndex).numWall(); 
	    // internalWallLenth[cellIndex].resize(numCellWalls);
	    // externalWallLength[cellIndex].resize(numCellWalls);
	    
	    
	    
	    // for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
	    //   size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();         
	    //   size_t restAddIndex=0;
	    //   if(numParameter()==3){         
	    //     if (T.wall(w2).cell1()->index()==cellIndex)
	    //       restAddIndex=1;
	    //     else if (T.wall(w2).cell2()->index()==cellIndex)
	    //       restAddIndex=2;
	    //     else{
	    //       std::cerr<<"something bad happened"<<std::endl;
	    //       exit(-1);
	    //     }
	    //   }
	    
	    
	    //   std::cerr<<wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex]<<"  ";
	    // }
	    // std::cerr<<"      ";
	    // for (size_t wallIndex=lengthInternalIndex; wallIndex<lengthInternalIndex+8; ++wallIndex) { 
	    
	    //   std::cerr<<cellData[0][wallIndex]<<"      ";
	    // }
	    
	    
	    // std::cerr<<" h "<<h<<std::endl;
	    
	    
	    
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
	      
	      // size_t restAddIndex=0;
	      // if(numParameter()==4 && parameter(3)==1){         
	      //   if (T.wall(w2).cell1()->index()==cellIndex)
	      //     restAddIndex=1;
	      //   else if (T.wall(w2).cell2()->index()==cellIndex)
	      //     restAddIndex=2;
	      //   else{
	      //     std::cerr<<"something bad happened"<<std::endl;
	      //     exit(-1);
	      //   }
	      // }
	      
	      
	      std::vector<double> restingLength(3);
	      if(numParameter()>3 && (parameter(3)==1 ||parameter(3)==-1)){ // double resting length        
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
	      
	      
	      
	      double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
					    (-restingLength[0]+restingLength[1]+restingLength[2])*
					    ( restingLength[0]-restingLength[1]+restingLength[2])*
					    ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
	      
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
	      
	      double Egreen[2][2];//E=0.5(C-I)
	      Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
	      Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
	      Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
	      Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);
	      
	      
	      double LeftCauchy[2][2]; // B=FFt
	      LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
	      LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
	      LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
	      LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];
	      
	      
	      double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
	      double tempS=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
	      StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/tempS));
	      StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/tempS;
	      StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/tempS;  
	      StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/tempS));
	      
	      double detTrue=StrainAlmansi[0][0]*StrainAlmansi[1][1]-StrainAlmansi[0][1]*StrainAlmansi[1][0];
	      double trTrue=StrainAlmansi[0][0]+StrainAlmansi[1][1];
	      
	      double strainValue1True=trTrue/2+std::sqrt(std::fabs(((trTrue*trTrue)/4)-detTrue));
	      double strainValue2True=trTrue/2-std::sqrt(std::fabs(((trTrue*trTrue)/4)-detTrue));
	      
	      //std::cerr<<strainValue1True<<"  "<<strainValue2True<<"      ";
	      
	      double det=Egreen[0][0]*Egreen[1][1]-Egreen[0][1]*Egreen[1][0];
	      double tr=Egreen[0][0]+Egreen[1][1];
	      
	      double strainValue1=tr/2+std::sqrt(std::fabs(((tr*tr)/4)-det));
	      double strainValue2=tr/2-std::sqrt(std::fabs(((tr*tr)/4)-det));
	      
	      
	      //std::cerr<<"strain values"<<strainValue1True<<"  "<<strainValue2True<<std::endl;
	      double strainRestLocal1[2]={0,0};
	      double strainRestLocal2[2]={0,0};
	      // is this the best way?
	      if (std::fabs(Egreen[1][0])>0.0000000000000001){
		strainRestLocal1[0]=strainValue1-Egreen[1][1];
		strainRestLocal1[1]=Egreen[1][0];
		strainRestLocal2[0]=strainValue2-Egreen[1][1];
		strainRestLocal2[1]=Egreen[1][0];
	      }
	      else if (std::fabs(Egreen[0][1])>0.0000000000000001){
		strainRestLocal1[0]=Egreen[0][1];
		strainRestLocal1[1]=strainValue1-Egreen[0][0];
		strainRestLocal2[0]=Egreen[0][1];
		strainRestLocal2[1]=strainValue2-Egreen[0][0];
	      }
	      else {
		strainRestLocal1[0]=1;
		strainRestLocal1[1]=0;
		strainRestLocal2[0]=0;
		strainRestLocal2[1]=1;
	      }
	      
	      double tempAn=std::sqrt(strainRestLocal1[0]*strainRestLocal1[0]+
				      strainRestLocal1[1]*strainRestLocal1[1]);
	      if (tempAn !=0 ){
		strainRestLocal1[0]/=tempAn;
		strainRestLocal1[1]/=tempAn;
	      }
	      else
		std::cerr<<"in growth.cc/strainTRBS strange strain eigenvector at cell "
			 <<cellIndex<<"  wall  "<<wallIndex<<std::endl;
	      
	      tempAn=std::sqrt(strainRestLocal2[0]*strainRestLocal2[0]+
			       strainRestLocal2[1]*strainRestLocal2[1]);
	      if (tempAn !=0 ){
		strainRestLocal2[0]/=tempAn;
		strainRestLocal2[1]/=tempAn;
	      }
	      else
		std::cerr<<"in growth.cc/strainTRBS strange strain eigenvector at cell "
			 <<cellIndex<<"  wall  "<<wallIndex<<std::endl;
	      
	      if (strainValue2>strainValue1){
		double temp=strainValue1;
		strainValue1=strainValue2;
		strainValue2=temp;
		temp=strainRestLocal1[0];
		strainRestLocal1[0]=strainRestLocal2[0];
		strainRestLocal2[0]=temp;        
		temp=strainRestLocal1[1];
		strainRestLocal1[1]=strainRestLocal2[1];
		strainRestLocal2[1]=temp;          
	      }
	      if (strainValue2True>strainValue1True){
		double temp=strainValue1True;
		strainValue1True=strainValue2True;
		strainValue2True=temp;
	      }
	      
	      std::vector<std::vector<double> > edgeRestLocal(3);
	      
	      for (size_t d=0; d< 3; ++d)
		edgeRestLocal[d].resize(3);
	      
	      
	      edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
	      edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
	      //edgeRestLocal[0][2]= positionLocal[][2]-positionLocal[][2];
	      
	      edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
	      edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];
	      //edgeRestLocal[1][2]= positionLocal[][2]-positionLocal[][2];
	      
	      edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
	      edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
	      //edgeRestLocal[2][2]= positionLocal[][2]-positionLocal[][2];
	      
	      
	      std:: vector<double> cosTet(3);
	      std:: vector<double> sinTet(3);
	      for (size_t j=0; j< 3; ++j){
		cosTet[j]=std::fabs((strainRestLocal1[0]*edgeRestLocal[j][0]+
				     strainRestLocal1[1]*edgeRestLocal[j][1])/
				    restingLength[j]);
		sinTet[j]=std::sqrt(std::fabs(1-cosTet[j]*cosTet[j]));
		
		if(cosTet[j]<0 || cosTet[j]>1)
		  std::cerr<<"in growth.cc/strainTRBS cosTet is wrong: "
			   <<cosTet[j]<<" in cell "<<cellIndex<<std::endl;
	      }
	      
	      
	      std::vector<std::vector<double> > restingComp(3);
	      for (size_t j=0; j< 3; ++j)
		restingComp[j].resize(2);
	      
	      for (size_t j=0; j< 3; ++j){
		restingComp[j][0]=restingLength[j]*cosTet[j];
		restingComp[j][1]=restingLength[j]*sinTet[j];
	      }
	      // std::cerr<<"cell Ind "<<cellIndex
	      //          <<",   strain1 "<<strainValue1True 
	      //          <<",   strain2  "<<strainValue2True 
	      //          <<",   strain Thresh "<<strainThreshold 
	      //          <<",  cellStrain "<<cellData[cellIndex][11]
	      //          <<",  velocity "<<cellData[cellIndex][22]
	      //          <<std::endl;
	      //std::cerr<<strainValue1True<<"  "<<strainValue2True <<"  "<<strainThreshold <<std::endl;              
	      if (strainValue1True>strainThreshold && strainValue2True<strainThreshold){
		
		cellData[cellIndex][growthInd]+=1;
		//std::cerr<<"growth 1"<<std::endl;  
		// double factor1=parameter(0)*h*(strainValue1True-strainThreshold)/(1-2*strainValue1True);
		// factor1=factor1+factor1*factor1;// first and second terms in taylor expansion
		
		double factor1=parameter(0)*h*(strainValue1True-strainThreshold);
		
		for (size_t j=0; j< 3; ++j)
		  restingComp[j][0]+=restingComp[j][0]*factor1;
		
	      }
	      
	      if (strainValue1True>strainThreshold && strainValue2True>strainThreshold){
		
		cellData[cellIndex][growthInd]+=2;
		//std::cerr<<"growth 2"<<std::endl;  
		// double 
		//   factor1=parameter(0)*h*(strainValue1True-strainThreshold)/(1-2*strainValue1True);
		// double 
		//   factor2=parameter(0)*h*(strainValue2True-strainThreshold)/(1-2*strainValue2True);
		
		// factor1=factor1+factor1*factor1;// first and second terms in taylor expansion
		// factor2=factor2+factor2*factor2;// first and second terms in taylor expansion
		
		double factor1=parameter(0)*h*(strainValue1True-strainThreshold);
		double factor2=parameter(0)*h*(strainValue2True-strainThreshold);
		
		
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
	      
	      
	      
	      if(numParameter()>3 && (parameter(3)==1 ||parameter(3)==-1)){ // double length
		cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallIndex]
		  =externalTemp-wallData[w2][wallLengthIndex];
		cellData[cellIndex][lengthInternalIndex + 2*wallIndex+1]=internalTemp;
		cellData[cellIndex][lengthInternalIndex + 2*wallIndexPlusOneMod]=internalTempPlusOne;
	      }
	      else{                  // single length
		
		// WITH AREA AVERAGING
		
		size_t wallGlobalInd= T.cell(cellIndex).wall(wallIndex) ->index();
		//std::cerr<<" normalization area factor before "<<mainWalls[wallGlobalInd][0]<<std::endl;
		if (mainWalls[wallGlobalInd][0]==0){
		  mainWalls[wallGlobalInd][0]=restingArea;
		  mainWalls[wallGlobalInd][1]=restingArea*externalTemp;
		  //std::cerr<<" normalization area factor middle "<<mainWalls[wallGlobalInd][0]<<std::endl;
		}
		else if (mainWalls[wallGlobalInd][0]!=0){
		  mainWalls[wallGlobalInd][0]+=restingArea;
		  mainWalls[wallGlobalInd][1]=
		    (mainWalls[wallGlobalInd][1]+restingArea*externalTemp);
		}
		
		if (internalWalls[cellIndex][wallIndex][0]==0){
		  internalWalls[cellIndex][wallIndex][0]=restingArea;
		  internalWalls[cellIndex][wallIndex][1]=restingArea*internalTemp;
		}
		else if (internalWalls[cellIndex][wallIndex][0]!=0){
		  internalWalls[cellIndex][wallIndex][0]+=restingArea;
		  internalWalls[cellIndex][wallIndex][1]=
		    (internalWalls[cellIndex][wallIndex][1]+restingArea*internalTemp);
		}         
		
		if (internalWalls[cellIndex][wallIndexPlusOneMod][0]==0){
		  internalWalls[cellIndex][wallIndexPlusOneMod][0]=restingArea;
		  internalWalls[cellIndex][wallIndexPlusOneMod][1]=restingArea*internalTempPlusOne;
		}
		else if (internalWalls[cellIndex][wallIndexPlusOneMod][0]!=0){
		  internalWalls[cellIndex][wallIndexPlusOneMod][0]+=restingArea;
		  internalWalls[cellIndex][wallIndexPlusOneMod][1]=
		    (internalWalls[cellIndex][wallIndexPlusOneMod][1]+restingArea*internalTempPlusOne);
		  
		}         
		
	      }
	    } // walls
	    
	    //std::cerr<<std::endl;
	    // updating wall length
	    
	  } // cells
      
      
      
      if(numParameter()==3 ||(numParameter()>3 && parameter(3)!=1 && parameter(3)!=-1)){ // single resting length
	for (size_t cellIndex=0; cellIndex< T.numCell(); cellIndex++)
	  for (size_t wallIndex=0; wallIndex< T.cell(cellIndex).numWall(); wallIndex++){
	    cellData[cellIndex][lengthInternalIndex + wallIndex]= 
	      internalWalls[cellIndex][wallIndex][1]
	      /internalWalls[cellIndex][wallIndex][0]; 
	    size_t wallGlobalInd=T.cell(cellIndex).wall(wallIndex)->index();
	    wallData[wallGlobalInd][wallLengthIndex]=
	      mainWalls[wallGlobalInd][1]/mainWalls[wallGlobalInd][0];
	  }
      }
      
      double totalCellStrain=0;
      double totalCellStress=0;
      double totalAreaRatio=0;
      double totalWallStrain=0;
      double totalWallStress=0;
      double totalCellEnergy=0;
      double totalWallEnergy=0;
      
      
      double totalEdge=0;
      
      
      for(size_t a=0; a<numWalls; ++a)
	totalEdge+=wallData[a][0];
      
      for(size_t a=0; a<numCells; ++a){
	totalCellStrain +=cellData[a][11];  
	totalCellStress +=cellData[a][7];  
	totalAreaRatio  +=cellData[a][19];  
	totalCellEnergy +=cellData[a][20];  
      }
      for(size_t a=0; a<numWalls; ++a)
	for(size_t a=0; a<numWalls; ++a){
	  size_t v1 = T.wall(a).vertex1()->index();
	  size_t v2 = T.wall(a).vertex2()->index();
	  size_t dimension = vertexData[v1].size();
	  double distance=0;
	  for( size_t d=0 ; d<dimension ; d++ )
	    distance += (vertexData[v1][d]-vertexData[v2][d])*
	      (vertexData[v1][d]-vertexData[v2][d]);
	  distance = std::sqrt(distance);
	  totalWallStrain +=0.5*(1-(wallData[a][1]/distance)); 
	  totalWallStress +=wallData[a][2];
	  totalWallEnergy +=wallData[a][1]*wallData[a][2]*wallData[a][2]*0.5/10;// frce^2/spring constant  
	}
      
      
      
      
      // double totalAreatmp=0;
      // for(size_t a=0; a<numCells; ++a)
      //   totalAreatmp+=cellData[a][19];
      // static double totalArea1=totalAreatmp; 
      
      
      
      
      
      //rescaling
      //double scaletmp=std::sqrt(totalArea1/totalAreatmp);//0.995;      
      // for(size_t a=0; a<numWalls; ++a)
      //   wallData[a][0]*=scaletmp;
      // for(size_t a=0; a<numCells; ++a)
      //   for(size_t b=42; b<cellData[a].size(); ++b)
      //     cellData[a][b]*=scaletmp;
      //    // rescaling       
      // for(size_t a=0; a<numWalls; ++a)
      //   wallData[a][1]*=scaletmp;
      
      
      //std::cout<<"cells  "<<" "<<scaletmp<<"  "
      //totalCellStress<<" "<<totalCellEnergy<<" "<<totalAreaRatio<<" "
      //          <<totalWallStrain<<" "<<totalWallStress<<" "<<totalWallEnergy<<" "
      //         <<totalAreatmp<<"  "<<totalArea1<<" "<<totalEdge<<std::endl;
      
      
      
      //deltat=0; 
      //std::cout<<growthtime<<"  "<<totalAreatmp<<"  "<<deltat<<std::endl;
    }
    
  }
  // reserve end
  
  
  VertexFromTRBScenterTriangulationMT::
  VertexFromTRBScenterTriangulationMT(std::vector<double> &paraValue, 
				      std::vector< std::vector<size_t> > 
				      &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=11 ) { 
      std::cerr << "Hypocotyl3D::VertexFromTRBScenterTriangulationMT::"
		<< "VertexFromTRBScenterTriangulationMT() "
		<< "Uses 11 parameters: "
		<< "0,1: young modulus(matrix and fibre) " 
		<< "2,3 : poisson ratio (longitudinal (MT) and transverse directions)"
		<< "4 : MF flag(0 constant material anisotropy ,1: material anisotropy via FiberModel " 
		<< "5 : neighbor weight for correlation " 
		<< "6 : max stress/strain  " 
		<< "7 : 2nd flag(0: plane strain, 1: plane stress) " 
		<< "8 : MT direction angle"
		<< "9 (MT update flag): 0:for no feedback or direct feedback by indices,"
		<< "                    1:for MT direction from 7th parameter TETA, 2:force to Stress,  "
		<< "                    3: force to Strain ,4:force to perp-strain "
		<< "10 : 1:independent resting length for the elements, otherwise: common length "
		<<" also see line 2835 where an ad-hoc thing is implemented for parameter(9)==1 "
		<< std::endl;
      
      exit(EXIT_FAILURE);
    }
    
    if( (indValue.size()!=2 && indValue.size()!=4) || 
	indValue[0].size()!=11 ||
	indValue[1].size()!=1 ||
	(indValue.size()==4 && (indValue[2].size()!=0 && indValue[2].size()!=1 && indValue[2].size()!=2
				&& indValue[2].size()!=3)) ||
	(indValue.size()==4 && (indValue[3].size()!=0 && indValue[3].size()!=1 && indValue[3].size()!=2 ))
	) { 
      std::cerr << "Hypocotyl3D::VertexFromTRBScenterTriangulationMT::"
		<< "VertexFromTRBScenterTriangulationMT() "
		<< "11 indices including (0)Wall length index and (1)MT direction initial "
		<< "index and (2)strain and (3)stress "
		<< "anisotropy indices and indices for storing (4)area ratio, (5)isotropic energy  "
		<< "and (6)anisotropic energy and (7)young_fiber and (8)Misses stress  and start indices  "
		<< "for storing (9)stress tensor(6 elements) and (10)normal vector to cell plane  " 
		<< "(3 elements) given in first level." 
		<< "Start of additional Cell variable indices (center(x,y,z) "
		<< "L_1,...,L_n, n=num vertex) is given in second level (typically at the end)." 
		<< "Optionally two additional levels can be given where the strain, perpendicular " 
		<< "to strain and  2nd strain can be stored in 3rd, stress and 2nd stress can be "
		<< "stored in 4th level directions/values(dx dy dz value) can be stored at given indices."
		<< "If no index given at 3rd level, strain will not be stored, if one index given strain "
		<< "will be stored and if two indices are given maximal and perpendicular strain will  "
		<< "be stored and if 3 indices are given 2nd strain direction values will be stored at 3rd index"
		<< "If no index given at 4th level, stress will not be stored, if one index given " 
		<< "stress will be stored and if two indices are given maximal and 2nd stress will be "
		<< "stored at 1st and 2nd index respectively "
		<< std::endl;
      exit(EXIT_FAILURE);
    }
  
    // Set the variable values
    setId("Hypocotyl3D::VertexFromTRBScenterTriangulationMT");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "Y_mod_M";   // Matrix Young modulus
    tmp[1] = "Y_mod_F";   // Fiber Young modulus
    tmp[2] = "P_ratio_L"; // Longitudinal Poisson ratio
    tmp[3] = "P_ratio_T"; // Transverse Poisson ratio
    tmp[4] = "MF flag";
    tmp[5] = "neigborweight";
    tmp[6] = "stressmax";
    tmp[7] = "Strain-Stress flag";
    tmp[8] = "TETA anisotropy";
    tmp[9] = "MT update flag";
    tmp[10] = "double length flag";
    
    setParameterId( tmp );
    
    if( parameter(4)!=-1 && (parameter(2)<0 || parameter(2)>=0.5 || parameter(3)<0 || parameter(3)>=0.5) ) {
      std::cerr << "Hypocotyl3D::VertexFromTRBScenterTriangulationMT::"
		<< "VertexFromTRBScenterTriangulationMT() "
		<< "poisson ratio must be 0 <= p < 0.5 " << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // if( parameter(4)!=0 && parameter(4)!=1 && parameter(4)!=2 ) {
    //   std::cerr << " VertexFromTRBScenterTriangulationMT::"
    // 	      << " VertexFromTRBScenterTriangulationMT() "
    // 	      << " 5th parameter must be 0 or 1  "
    //             << " 0: constant material anisotropy, 1: material anisotropy via FiberModel  2: for energy landscape(see the code, do not use until you are sure)" << std::endl;
    //   exit(0);
    // }
    
    
    if( parameter(7)!=0 && parameter(7)!=1 ) {
      std::cerr << "Hypocotyl3D::VertexFromTRBScenterTriangulationMT::"
		<< "VertexFromTRBScenterTriangulationMT() "
		<< "6th parameter must be 0 or 1(0:plane strain, 1:plane stress) " << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( parameter(9)!=0 && parameter(9)!=1 && parameter(9)!=2 && parameter(9)!=3 && parameter(9)!=4) {
      std::cerr << "Hypocotyl3D::VertexFromTRBScenterTriangulationMT::"
		<< "VertexFromTRBScenterTriangulationMT() "
		<< "10th parameter must be 0/1/2/3/4"
		<< "0: for no feedback or direct feedback by indices "
		<< "1: for MT direction from 7th parameter TETA "
		<< "2: force to Stress "
		<< "3: force to Strain "
		<< "4: force to perp-strain " << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  

void VertexFromTRBScenterTriangulationMT::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
 

  
  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t numVertices = T.numVertex();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t comIndex = variableIndex(1,0);
  size_t lengthInternalIndex = comIndex+dimension;

  double neighborweight=parameter(5);
  double strainEpcilon =0.000001;
  double stressEpcilon =0.000001;    
  // double TotalVolume=0;
  // double deltaVolume=0;
  // for(size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex){ // stimating volume for equilibrium 
  //   TotalVolume +=std::sqrt(vertexData[vertexIndex][0]*vertexData[vertexIndex][0] +
  //                           vertexData[vertexIndex][1]*vertexData[vertexIndex][1] +
  //                           vertexData[vertexIndex][2]*vertexData[vertexIndex][2] );
  // }
  // deltaVolume=TotalVolume-cellData[0][25];
  // cellData[0][25]=TotalVolume;
  // cellData[0][24]=deltaVolume;
 
  size_t MTindex           =variableIndex(0,1);	 
  size_t isoEnergyIndex    =variableIndex(0,5);	
  size_t anisoEnergyIndex  =variableIndex(0,6);	
  size_t youngLIndex       =variableIndex(0,7);	
   

  double youngMatrix= parameter(0);    
  double youngFiber = parameter(1); 
  double poissonL   = parameter(2);    
  double poissonT   = parameter(3);
  double TETA       = parameter(8);  
  
  // static double tTotal=0, tDiag=0;
  // static double tRest=0, tRest2=0, tRest3=0, tRest4=0;

  // clock_t cpuTime0 ,cpuTimef, cpuTime1 ,cpuTime2, cpuTime3 ,cpuTime4;
  // cpuTime0=clock();

  for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    size_t numWalls = T.cell(cellIndex).numWall();
    
  

    if(  T.cell(cellIndex).numVertex()!= numWalls ) {
     
      std::cerr << "VertexFromTRBScenterTriangulationMT::derivs() same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
    double youngL=1;
    double youngT=1;
    


    if( parameter(4)==1){  // material anisotropy via FiberModel
      youngL = cellData[cellIndex][youngLIndex]; 
      youngT = 2*youngMatrix+youngFiber-youngL; 
    }
    
    if( parameter(4)==0 ){ // constant anisotropic material
      youngL = youngMatrix+youngFiber;
      youngT = youngMatrix; 
    }

        
    if( parameter(4)==-1 ){ //hypocotyl Siobhan
      youngL = 50;
      youngT = 50; 
      poissonL=0.2;
      poissonT=0.2;

      if(cellData[cellIndex][37]==-1){ // epidermis
	youngL *=parameter(0);
	youngT *=parameter(0);
      }
      if(cellData[cellIndex][37]==-2){ // inner
	youngL *=parameter(1); //1, 6 // hoop
	youngT *=1;                   // axial
      }
      if(cellData[cellIndex][37]==-3){  // anti. axial
	youngL *=parameter(2); //0.1, 0.5, 1
	youngT *=parameter(2); //0.1, 0.5, 1
      }
      if(cellData[cellIndex][37]==-4){  // anti transverse
	youngL *=parameter(3); //1, 10
	youngT *=parameter(3); //1, 10
      }
    }
    

      
      
    if( parameter(4)==2){  // for varrying material anisotropy with constant overall stiffness for energy landscape
      youngL =youngFiber; 
      youngT =youngMatrix-youngL;  // here youngMatrix is total stiffness
    }
    if( parameter(4)==3){  // for varrying material anisotropy with constant overall stiffness for energy landscape       
      double totalElast=youngMatrix;
      double Maniso=youngFiber;
      youngT =((1-Maniso)/(2-Maniso))*totalElast; 
      youngL=totalElast-youngT;
      
      //youngL =youngMatrix+youngFiber; 
      // youngT =youngMatrix+std::sqrt(
      //                               ((youngFiber-totalElast)*4*totalElast)
      //                               /
      //                               ((2-3.1415)*3.1415)
      //                               );  // here youngMatrix is total stiffness
      //youngT =youngMatrix+totalElast*std::sqrt(
      //                                         -(std::log((youngFiber/totalElast)-0.55))
      //                                         /6
      //                                         );
      //youngT =youngMatrix+2*(totalElast-youngFiber)/(3.1415-2);
      
      // youngL=youngMatrix+youngFiber; 
      // youngT =youngMatrix+(youngFiber-slope*totalElast)/(1-2*slope); 
    }
    
    // if( parameter(4)==5){  // for domain dependent material properties - Hypocotyl bending 3d (internal tissue softer)    
      //   if (cellData[cellIndex][38]==-1){
      //     youngL =youngMatrix+youngFiber; 
      //     youngT =youngMatrix;  
          
      //   }
      //   else{
      //     youngL =0.2*(youngMatrix+youngFiber); 
      //     youngT =0.2*youngMatrix;  
      //   }

      //   if (cellData[cellIndex][5]==5){
      //     youngL =0.01*(youngMatrix+youngFiber); 
      //     youngT =0.01*youngMatrix;  
      //   }
      // }

      if( parameter(4)==5){// for pavement cell resolution addaptive stiffness
        youngFiber=cellData[cellIndex][anisoEnergyIndex];
        double fiberL=cellData[cellIndex][youngLIndex];
        youngL = youngMatrix+fiberL;
        youngT = youngMatrix+youngFiber-fiberL;
       
        // if(cellData[cellIndex][30]==1){
        //   youngL = 20;
        //   youngT = 20;
        
        // }
      }
      
      if( parameter(4)==6){  // material anisotropy via FiberModel and loosening adhoc based on auxin
        youngL = cellData[cellIndex][youngLIndex]; 
        youngT = 2*youngMatrix+youngFiber-youngL; 
        
        double Kconc=0.005;
        double Nc1=2;
        double conc=cellData[cellIndex][13];
        double frac=(1-std::pow(conc,Nc1)/(std::pow(Kconc,Nc1)+std::pow(conc,Nc1)));
        youngL*=0.25+0.75*frac;
        youngT*=0.25+0.75*frac;
        //if(cellIndex==0) std::cerr<<cellData[cellIndex][13]<<"  "<<youngL<<"  " <<youngT<<std::endl;
      }
      
      if( parameter(4)==7){  // material anisotropy via FiberModel and destroying fibers based on auxin
	youngL = cellData[cellIndex][youngLIndex]-youngMatrix;
	youngT = 2*youngMatrix+youngFiber-youngL-2*youngMatrix;  
        
        double Kconc=0.005;
        double Nc1=2;
        double conc=cellData[cellIndex][13];
        double frac=(1-std::pow(conc,Nc1)/(std::pow(Kconc,Nc1)+std::pow(conc,Nc1)));
        youngL*=frac;
        youngT*=frac;
        youngL+=youngMatrix;
        youngT+=youngMatrix;
        //cellData[cellIndex][youngLIndex]=youngL; 
	
      }
      if( parameter(4)==8){  // material anisotropy via FiberModel and destroying fibers and a fraction of matrix based on auxin
	youngL = cellData[cellIndex][youngLIndex]-youngMatrix;
	youngT = 2*youngMatrix+youngFiber-youngL-2*youngMatrix;  
        
        double Kconc=0.005;
        double Nc1=2;
        double conc=cellData[cellIndex][13];
        double frac=(1-std::pow(conc,Nc1)/(std::pow(Kconc,Nc1)+std::pow(conc,Nc1)));
        youngL*=frac;
        youngT*=frac;
        youngL+=youngMatrix*(0.2+0.8*frac);
        youngT+=youngMatrix*(0.2+0.8*frac);
        //cellData[cellIndex][youngLIndex]=youngL; 
	
      }
      

                     
      // If( parameter(4)==-2){  // for heterogeneous stiffness(adhoc)
      //   double hFactor=0;
      //   double Hthreshold=0.01;
      //   if (std::fabs(cellData[cellIndex][youngLIndex])>Hthreshold){
      //     hFactor=Hthreshold;
      //   }
      //   else{
      //     hFactor=std::fabs(cellData[cellIndex][youngLIndex]); // take the heterogeneous info from this index
      //   }
      //   hFactor *=std::fabs(parameter(4)); // factor for heterogeneity      
      
      //   youngL = youngMatrix+youngFiber-hFactor*(youngMatrix+youngFiber); 
      //   youngT = youngMatrix-hFactor*(youngMatrix); 
      // }
      
      // if (cellData[cellIndex][31]==10){
      //   youngL = 150;
      //   youngT = 10;
      // }
     
      
      // if( parameter(4)<0){  // for spatial elasticity (ad-hoc)
      //   double x0=90;
      //   double y0=120;
        
      //   double xx=cellData[cellIndex][comIndex]-x0;
      //   double yy=cellData[cellIndex][comIndex+1]-y0;
        
      //   youngL = youngMatrix+youngFiber*(1-std::exp(parameter(4)*(xx*xx+yy*yy)));
      //   youngT = youngMatrix; 
      //   cellData[cellIndex][13]=youngL;
      // }
      
      // ad-hoc for regional loosening
      // if ( std::sqrt(cellData[cellIndex][comIndex  ]*cellData[cellIndex][comIndex  ]
      //                +cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])<30)
      //   cellData[cellIndex][youngLIndex]=20;
      // else
      //   cellData[cellIndex][youngLIndex]=200;
      
      
      
      // ad-hoc for 3d marcus
      
      // if(cellData[cellIndex][25]==-3){  // sidewalls (anticlinals)
      //   youngL = youngMatrix;
      //   youngT = youngMatrix; 
      //   if(cellData[cellIndex][29]==1){
      
      //   }
      // }
      
      //youngL = (1-cellData[cellIndex][30]/0.6)*youngL;
      //youngT = (1-cellData[cellIndex][30]/0.6)*youngT;
      
      // double sEpi=10;
      
      // if(cellData[cellIndex][27]==-2){ //L1 epidermis
      //   youngL = sEpi*youngL;
      //   youngT = sEpi*youngT; 
      // }
      // if(cellData[cellIndex][30]==10){ // L1 or L2 patch
      //   youngL = 0.1*youngL;
      //   youngT = 0.1*youngT; 
      // }
      
      // if(cellData[cellIndex][30]==5){ // b cells
      //   youngL = 0.5*youngL;
      //   youngT = 0.5*youngT; 
      // }
      // if(cellData[cellIndex][30]==2){ // c cells
      //   youngL = 0.8*youngL;
      //   youngT = 0.8*youngT; 
      // }
      
      
      // if(cellData[cellIndex][25]==-3 && cellData[cellIndex][29]==1){  // L1 anticlinals
      //   youngL = sEpi*youngL;
      //   youngT = sEpi*youngT; 
      // }
      
      // if(cellData[cellIndex][30]==400){
      //   youngL = 0.1*youngL;
      //   youngT = 0.1*youngT;
      // }
      
      
      // ad-hoc for root hair cell files and ....
      // double polar=180*std::acos(cellData[cellIndex][comIndex]/
      //                            std::sqrt(cellData[cellIndex][comIndex  ]*cellData[cellIndex][comIndex  ]+
      //                                      cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])
      //                            );
      // if (std::fmod(std::floor(polar/12),2)==1){
      //   youngL =110 ;
      //   youngT =110 ; 
      
      // }
      
      // if (std::sqrt(cellData[cellIndex][comIndex]*cellData[cellIndex][comIndex]+
      //               (cellData[cellIndex][comIndex+1]-60)*(cellData[cellIndex][comIndex+1]-60))<35){
      //   youngL =5 ;
      //   youngT =5 ; 
      
      // }
      
      // if (std::sqrt(cellData[cellIndex][comIndex]*cellData[cellIndex][comIndex]+
      //                         cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])<10){
      //   youngL =100 ;
      //   youngT =100 ; 
      
      // }
      



   
    double lambdaL, mioL, lambdaT, mioT; // ,lambdaTmatrix, mioTmatrix;
    
    if (parameter(7)==0){      
      // Lame coefficients based on plane strain (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/((1+poissonL)*(1-2*poissonL));
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/((1+poissonT)*(1-2*poissonT));
      mioT=youngT/(2*(1+poissonT));

   
    } 
    else{      
      // Lame coefficients based on plane stress (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/(1-poissonL*poissonL);
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/(1-poissonT*poissonT);
      mioT=youngT/(2*(1+poissonT));
     
      // lambdaTmatrix=youngMatrixA*poissonT/(1-poissonT*poissonT);
      // mioTmatrix=youngMatrixA/(2*(1+poissonT));
    }
    
    // Lame coefficients based on delin. paper (for 2D 0<poisson<1)
    // double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
    // double mioL=youngL/(1+poissonL);
    // double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
    // double mioT=youngT/(1+poissonT);
    

    

    double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double StressCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double normalGlob[3]={0,0,0};
    double TotalCellRestingArea=0;
    double TotalCellArea=0;
    double EnergyIso=0;     
    //double EnergyIsoFiber=0;            
    double EnergyAniso=0;
    double strainZ=0;

    if ( parameter(9)==1 ) {  // aniso direction from TETA 
      
      cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);    
      cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
      cellData[cellIndex][variableIndex(0,1)+2]=0;
      if(parameter(10)!=100){
	if(cellIndex==0 || cellIndex==2 ){
	  cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);    
	  cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
	  cellData[cellIndex][variableIndex(0,1)+2]=0;
	}
	if(cellIndex==1 || cellIndex==3 ){
	  cellData[cellIndex][variableIndex(0,1)]=std::cos(parameter(10));    
	  cellData[cellIndex][variableIndex(0,1)+1]=std::sin(parameter(10));
	  cellData[cellIndex][variableIndex(0,1)+2]=0;
	}
      }
    }
    
  

    // Aniso vector in current shape in global coordinate system
    double  AnisoCurrGlob[3]=
      { cellData[cellIndex][variableIndex(0,1)],  
        cellData[cellIndex][variableIndex(0,1)+1],
        cellData[cellIndex][variableIndex(0,1)+2]
      };
    
    //cpuTime1=clock(); 
      // One triangle per 'vertex' in cyclic order
      for (size_t wallindex=0; wallindex<numWalls; ++wallindex) {
        
        
      size_t kPlusOneMod = (wallindex+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
      size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //size_t w1 = internal wallindex
      size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //size_t w3 = internal wallindex+1

      // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      //position[0][2] z for vertex 1 of the current element
      

      // Resting lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)

      // size_t restAddIndex;
      // if (T.wall(w2).cell1()->index()==cellIndex)
      //   restAddIndex=1;
      // else if (T.wall(w2).cell2()->index()==cellIndex)
      //   restAddIndex=2;
      // else{
      //   std::cerr<<"something bad happened"<<std::endl;
      //   exit(-1);
      // }

      double restingLength[3];
      if(parameter(10)==1) { // double resting length  
        restingLength[0] =cellData[cellIndex][lengthInternalIndex + 2*wallindex+1];
        //restingLength[1] = wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex];
        restingLength[1] =wallData[w2][wallLengthIndex]+
          cellData[cellIndex][lengthInternalIndex+2*numWalls+wallindex];
        restingLength[2] =cellData[cellIndex][lengthInternalIndex + 2*kPlusOneMod];
      }      
      else{                 // single resting length
        restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
        restingLength[1] = wallData[w2][wallLengthIndex];
        restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
      }
      
      // Lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
      double length[3]=
        { std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) ),
          
          T.wall(w2).lengthFromVertexPosition(vertexData),
          
          std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) )
        };
      
      
      // Area of the element (using Heron's formula)                                      
      double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                    (-restingLength[0]+restingLength[1]+restingLength[2])*
                                    ( restingLength[0]-restingLength[1]+restingLength[2])*
                                    ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
      
      //double currentArea=std::sqrt( ( length[0]+length[1]+length[2])*
      //                              (-length[0]+length[1]+length[2])*
      //                              ( length[0]-length[1]+length[2])*
      //                              ( length[0]+length[1]-length[2])  )*0.25;

      //cpuTime3=clock();
     
      
      //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
      //std::vector<double> Angle(3);
      // can be changed by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
      // Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]
      //                       -restingLength[1]*restingLength[1])/
      //                      (restingLength[0]*restingLength[2]*2)    );
      // Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]
      //                       -restingLength[2]*restingLength[2])/
      //                      (restingLength[0]*restingLength[1]*2)    );
      // Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]
      //                       -restingLength[0]*restingLength[0])/
      //                      (restingLength[1]*restingLength[2]*2)    );
      
      // //Tensile Stiffness
      // std::vector<double>  tensileStiffness(3);
      // //double tensileStiffness[3];
      // double temp = 1.0/(restingArea*16);                                      
      // std::vector<double> cotan(3);
      // cotan[0] = 1.0/std::tan(Angle[0]);
      // cotan[1] = 1.0/std::tan(Angle[1]);
      // cotan[2] = 1.0/std::tan(Angle[2]);


      double cosAngle[3]=
        { (  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]
              -restingLength[1]*restingLength[1])/
             (restingLength[0]*restingLength[2]*2)    ),

          (  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]
              -restingLength[2]*restingLength[2])/
             (restingLength[0]*restingLength[1]*2)    ),

          (  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]
              -restingLength[0]*restingLength[0])/
             (restingLength[1]*restingLength[2]*2)    )
        };
      
      double temp = 1.0/(restingArea*16);                                      
      double cotan[3]=
        { cosAngle[0]/std::sqrt(1-cosAngle[0]*cosAngle[0]),
          cosAngle[1]/std::sqrt(1-cosAngle[1]*cosAngle[1]),
          cosAngle[2]/std::sqrt(1-cosAngle[2]*cosAngle[2])
        };
      
      
      //the force is calculated based on Transverse coefficients
      //Longitudinal coefficients are considered in deltaF
      //Tensile Stiffness
      
      double tensileStiffness[3]=
        { (2*cotan[2]*cotan[2]*(lambdaT+2*mioT)+2*mioT)*temp,
          (2*cotan[0]*cotan[0]*(lambdaT+2*mioT)+2*mioT)*temp,
          (2*cotan[1]*cotan[1]*(lambdaT+2*mioT)+2*mioT)*temp
        };
      
      //Angular Stiffness
      double angularStiffness[3]=
        { (2*cotan[1]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp,
          (2*cotan[0]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp,
          (2*cotan[0]*cotan[1]*(lambdaT+2*mioT)-2*mioT)*temp
        };
      
      //Calculate biquadratic strains  
      double Delta[3]=
        { (length[0])*(length[0])-(restingLength[0])*(restingLength[0]),
          (length[1])*(length[1])-(restingLength[1])*(restingLength[1]),
          (length[2])*(length[2])-(restingLength[2])*(restingLength[2])
        };
      
      //Area of the element (using Heron's formula)                                      
      double Area=std::sqrt( ( length[0]+length[1]+length[2])*
                             (-length[0]+length[1]+length[2])*
                             ( length[0]-length[1]+length[2])*
                             ( length[0]+length[1]-length[2])  )*0.25;
      
      
      
      // calculating the angles between shape vectors and anisotropy direction in resting shape when anisotropy vector is provided in current shape
      
      //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
      // double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
      //                                  (length[0]*length[1]*2)    );
      
      // double Qa=std::cos(CurrentAngle1)*length[0];
      // double Qc=std::sin(CurrentAngle1)*length[0];
      // double Qb=length[1];
      
      double CurrentAngle1=(length[0]*length[0]
                            +length[1]*length[1]
                            -length[2]*length[2])/
        (length[0]*length[1]*2);
      
      double Qa=(CurrentAngle1)*length[0];
      double Qc=std::sqrt(1-CurrentAngle1*CurrentAngle1)*length[0];
      double Qb=length[1];
      
      // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
      
      
      double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
                                        {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
                                        { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
    
      
      // Local coordinates of the resting shape ( counterclockwise )
      // double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]
      //                                   +restingLength[1]*restingLength[1]
      //                                   -restingLength[2]*restingLength[2])
      //                                  /(restingLength[0]*restingLength[1]*2)    );
      
      
      // double Pa=std::cos(RestingAngle1)*restingLength[0];
      // double Pc=std::sin(RestingAngle1)*restingLength[0];
      // double Pb=restingLength[1];
      
      
      double RestingAngle1=  (restingLength[0]*restingLength[0]
                              +restingLength[1]*restingLength[1]
                              -restingLength[2]*restingLength[2])/
        (restingLength[0]*restingLength[1]*2);
      
      
      double Pa=(RestingAngle1)*restingLength[0];
      double Pc=std::sqrt(1-RestingAngle1*RestingAngle1)*restingLength[0];
      double Pb=restingLength[1];
      
      // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                        {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                        { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
      
      
      
      // //Strain tensor  (clockwise ordering of nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );
      
      // double Qa=std::cos(CurrentAngle2)*length[2];
      // double Qb=length[1];
      // double Qc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVectorCurrent[3][2]={ {  0   ,       1/Qc      }, 
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) },       
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) }  };
      
      
      //cpuTime4=clock();
      //tRest4 +=cpuTime4-cpuTime3;
      
      
      // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
      
      
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
  
            
   
      
      // rotating the anisotropy vector from global coordinate system to the local one in the current shape
      double AnisoCurrLocal[3]=
        {  rotation[0][0]*AnisoCurrGlob[0]+
           rotation[1][0]*AnisoCurrGlob[1]+
           rotation[2][0]*AnisoCurrGlob[2],
           
           rotation[0][1]*AnisoCurrGlob[0]+
           rotation[1][1]*AnisoCurrGlob[1]+
           rotation[2][1]*AnisoCurrGlob[2],
           
           rotation[0][2]*AnisoCurrGlob[0]+
           rotation[1][2]*AnisoCurrGlob[1]+
           rotation[2][2]*AnisoCurrGlob[2]
        };
      


      // this part is  replaced by def. gradient transformation  
      // // Center of Mass for current shape in local coordinate
      // double CMCurrentLocal[2]={(Qa+Qb)/3, Qc/3};
      // // Tip of the ansotropy vector drown from the center of mass in the current shape
      // double  ACurrLocal[2] = {CMCurrentLocal[0]+AnisoCurrLocal[0],CMCurrentLocal[1]+AnisoCurrLocal[1]};
      // // Baricentric Coordinates of tip of anisotropy vector in the current shape wihich is equivalent to the baricentric coordinate of the corresponding point in the resting shape
      // double ABari[3]=
      //   {  ShapeVectorCurrent[0][0]*ACurrLocal[0]+
      //      ShapeVectorCurrent[0][1]*ACurrLocal[1]+
      //      ShapeVectorCurrent[0][2],
      
      //      ShapeVectorCurrent[1][0]*ACurrLocal[0]+
      //      ShapeVectorCurrent[1][1]*ACurrLocal[1]+
      //      ShapeVectorCurrent[1][2],
      
      //      ShapeVectorCurrent[2][0]*ACurrLocal[0]+
      //      ShapeVectorCurrent[2][1]*ACurrLocal[1]+
      //      ShapeVectorCurrent[2][2]
      //   };
      // //std::cerr<< "cell "<< cellIndex<< " baricentric coor of tip of anisoVector in local current "<<ABari[0]<<"  "<<ABari[1]<<"  "<<ABari[2]<< std::endl;
      // // Local coordinates of tip of AnisoVector in the resting shape from multyplying ABari by position matrix in the local resting shape coordinates           
      // double ARestLocal[2]=
      //   { Pa*ABari[0]+Pb*ABari[2],
      //     Pc*ABari[0]
      //   };
      // //std::cerr<< "cell "<< cellIndex<< " tip of anisoVector in local rest "<<ARestLocal[0]<<"  "<<ARestLocal[1]<< std::endl;     
      // // Center of Mass for resting shape in local coordinate
      // double CMRestingLocal[2]={(Pa+Pb)/3, Pc/3};
      // // Aniso Vector in the resting shape in local coordinate system
      // double AnisoRestLocal[2]={ARestLocal[0]-CMRestingLocal[0],ARestLocal[1]-CMRestingLocal[1]};
      // //std::cerr<< "cell "<< cellIndex<< " anisoVector Rest "<<AnisoRestLocal[0]<<"  "<<AnisoRestLocal[1]<<"  "<<AnisoRestLocal[2]<< std::endl;    
      
      


      //Anisotropic Correction is based on difference between Lam Coefficients of Longitudinal and Transverse dirrections:
      // double deltaLam=AnisoMeasure*(lambdaL-lambdaT);
      // double deltaMio=AnisoMeasure*(mioL-mioT);
      double deltaLam=lambdaL-lambdaT;
      double deltaMio=mioL-mioT;  

      // double deltaLamIsoFiber=lambdaT-lambdaTmatrix;
      // double deltaMioIsoFiber=mioT-mioTmatrix;  
 
      //cpuTime3=clock();
    
     
      //Angles between anisotropy vector and shape vectors for calculating the terms like a.Di , 
      //teta(k) = acos((dot(Anisorest,Dk))/(norm(Anisorest)*norm(Dk))),
      // std::vector<double> teta(3);
      // teta[0] = std::acos(  (ShapeVectorResting[0][0]*AnisoRestLocal[0]
      //                        +ShapeVectorResting[0][1]*AnisoRestLocal[1])/
      //                       std::sqrt(ShapeVectorResting[0][0]*ShapeVectorResting[0][0]
      //                                 +ShapeVectorResting[0][1]*ShapeVectorResting[0][1]+0.00000001) );
      
      // teta[1] = std::acos(  (ShapeVectorResting[1][0]*AnisoRestLocal[0]
      //                        +ShapeVectorResting[1][1]*AnisoRestLocal[1])/
      //                       std::sqrt(ShapeVectorResting[1][0]*ShapeVectorResting[1][0]
      //                                 +ShapeVectorResting[1][1]*ShapeVectorResting[1][1]+0.00000001) );
      
      // teta[2] = std::acos(  (ShapeVectorResting[2][0]*AnisoRestLocal[0]
      //                        +ShapeVectorResting[2][1]*AnisoRestLocal[1])/
      //                       std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]
      //                                 +ShapeVectorResting[2][1]*ShapeVectorResting[2][1]+0.00000001) );
      //cpuTime4=clock();
      //tRest4 +=cpuTime4-cpuTime3;
     
     
    
       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
      // strain tensor in resting shape E=0.5(FtF-I)
      // trE
      // B=FFt
      // axa (direct product of aniso vector in resting shape)
      // atEa
      // E(axa) and (axa)E
      
      
     

      
      double positionLocal[3][2]={ {Qa , Qc}, 
                                   {0  , 0 },  
                                   {Qb , 0 }  };
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F=Sigma i Qi x Di   <---------------------------------------
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];        
      } 




      // transforming aniso vector in local coordinate from current to resting shape using inverse of deformation gradient
      double AnisoRestLocal[2]={ DeformGrad[1][1] *AnisoCurrLocal[0]-DeformGrad[0][1]*AnisoCurrLocal[1],
                                 -DeformGrad[1][0]*AnisoCurrLocal[0]+DeformGrad[0][0]*AnisoCurrLocal[1]};
                
      double AnisoMeasure=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]);
      // std::cerr<< "cell "<< cellIndex<<" AnisoMeasure "<< AnisoMeasure<<std::endl;
      
      // choosing a random normalized dirrection for anisotropy if AnisoVector is close to perpendicular to the cell plane
      if ( AnisoMeasure<0.0001) {
        double randomAngle=((rand()%360)*2*3.14159265)/360;
        
        AnisoRestLocal[0]=std::cos(randomAngle);
        AnisoRestLocal[1]=std::sin(randomAngle);
      }  
      else {// normalizing the anisoVector if it is not random
        AnisoRestLocal[0]=AnisoRestLocal[0]/AnisoMeasure;
        AnisoRestLocal[1]=AnisoRestLocal[1]/AnisoMeasure;        
      }
      
      double directAniso[2][2]={{AnisoRestLocal[0]*AnisoRestLocal[0],AnisoRestLocal[0]*AnisoRestLocal[1]},
                                {AnisoRestLocal[1]*AnisoRestLocal[0],AnisoRestLocal[1]*AnisoRestLocal[1]}};



      double LeftCauchy[2][2]= // B=FFt
        { { DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1],
            DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1]
          },
          { DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1],
            DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1]
          }
        };
      
      double Egreen[2][2]= //E=0.5(C-I)
        { { 0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1),
            0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1])   
          },
          { 0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]), 
            0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1)
          } 
        };
     
      double trE=Egreen[0][0]+Egreen[1][1];
      
      double E2[2][2]= // used for energy calculation only
        { { Egreen[0][0]*Egreen[0][0]+Egreen[0][1]*Egreen[1][0],
            Egreen[0][0]*Egreen[0][1]+Egreen[0][1]*Egreen[1][1]
          },
          { Egreen[1][0]*Egreen[0][0]+Egreen[1][1]*Egreen[1][0],
            Egreen[1][0]*Egreen[0][1]+Egreen[1][1]*Egreen[1][1]
          }
        };

      double I2=E2[0][0]+E2[1][1]; //trE2 used for energy calculation only

      
      double I5=AnisoRestLocal[0]*AnisoRestLocal[0]*E2[0][0]   //atE2a used for energy calculation only
        +AnisoRestLocal[0]*AnisoRestLocal[1]*(E2[0][1]+E2[1][0])
        +AnisoRestLocal[1]*AnisoRestLocal[1]*E2[1][1];
      
    
     
      temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
      double StrainAlmansi[2][2]= // e=0.5(1-B^-1)  True strain tensor
        { { 0.5*(1-(LeftCauchy[1][1]/temp)) , 0.5*LeftCauchy[0][1]/temp },
          { 0.5*LeftCauchy[1][0]/temp       , 0.5*(1-(LeftCauchy[0][0]/temp))}
        };
      
      if(StrainAlmansi[0][0] != StrainAlmansi[0][0] ||
         StrainAlmansi[1][1] != StrainAlmansi[1][1] ||
         StrainAlmansi[0][1] != StrainAlmansi[0][1] ||
         StrainAlmansi[1][0] != StrainAlmansi[1][0] ) 
        std::cerr<<std::endl<<" strain is wrong "<<StrainAlmansi[0][0]<<StrainAlmansi[0][0]
                 <<StrainAlmansi[0][0] << StrainAlmansi[0][0] <<"   Q "<<restingLength[0]<<" "<<restingLength[1]<<" "<<restingLength[2]<<"    P "<<Pa<<" "<<Pb<<" "<<Pc<<std::endl;
      double atEa=AnisoRestLocal[0]*AnisoRestLocal[0]*Egreen[0][0]
        +AnisoRestLocal[0]*AnisoRestLocal[1]*(Egreen[0][1]+Egreen[1][0])
        +AnisoRestLocal[1]*AnisoRestLocal[1]*Egreen[1][1];

      double I4=atEa;


    

  
      double Eaa[2][2]=
        { { Egreen[0][0]*directAniso[0][0]+Egreen[0][1]*directAniso[1][0],
            Egreen[0][0]*directAniso[0][1]+Egreen[0][1]*directAniso[1][1]
          },
          { Eaa[1][0]= Egreen[1][0]*directAniso[0][0]+Egreen[1][1]*directAniso[1][0],             
            Eaa[1][1]= Egreen[1][0]*directAniso[0][1]+Egreen[1][1]*directAniso[1][1]
          }
        };        

      double aaE[2][2]=
        { { directAniso[0][0]*Egreen[0][0]+directAniso[0][1]*Egreen[1][0],        
            directAniso[0][0]*Egreen[0][1]+directAniso[0][1]*Egreen[1][1]
          },        
          { directAniso[1][0]*Egreen[0][0]+directAniso[1][1]*Egreen[1][0],              
            directAniso[1][0]*Egreen[0][1]+directAniso[1][1]*Egreen[1][1]        
          }
        };

      double B2[2][2]=// LeftCauchy^2
        { { LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0],
            LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1]
          },
          { LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0],
            LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1]
          }
        };
      
      double areaFactor=restingArea/Area; // 1/detF
      //double areaFactor=Area/restingArea; // detF
      
      double Sigma[2][2]= // true stress tensor (isotropic term) 
        { { areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][0]+(mioT)*B2[0][0]),
            areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][1]+(mioT)*B2[0][1])
          },      
          { areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][0]+(mioT)*B2[1][0]),      
            areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][1]+(mioT)*B2[1][1])
          }
        };
      

      // ENERGY MODELS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


      // double deltaS[2][2]; // based on Delin. paper
      // deltaS[0][0]=deltaLam*(trE*directAniso[0][0]+atEa)+(2*deltaMio)*(Eaa[0][0]+aaE[0][0])-(deltaLam+2*deltaMio)*atEa*directAniso[0][0];
      // deltaS[1][0]=deltaLam*(trE*directAniso[1][0]     )+(2*deltaMio)*(Eaa[1][0]+aaE[1][0])-(deltaLam+2*deltaMio)*atEa*directAniso[1][0];
      // deltaS[0][1]=deltaLam*(trE*directAniso[0][1]     )+(2*deltaMio)*(Eaa[0][1]+aaE[0][1])-(deltaLam+2*deltaMio)*atEa*directAniso[0][1];
      // deltaS[1][1]=deltaLam*(trE*directAniso[1][1]+atEa)+(2*deltaMio)*(Eaa[1][1]+aaE[1][1])-(deltaLam+2*deltaMio)*atEa*directAniso[1][1];

      double deltaS[2][2]= // based on  equipartitioning
        { { (deltaLam/2)*(trE*directAniso[0][0]+atEa)+(deltaMio)*(Eaa[0][0]+aaE[0][0]),
            (deltaLam/2)*(trE*directAniso[0][1]     )+(deltaMio)*(Eaa[0][1]+aaE[0][1])
          },
          { (deltaLam/2)*(trE*directAniso[1][0]     )+(deltaMio)*(Eaa[1][0]+aaE[1][0]),      
            (deltaLam/2)*(trE*directAniso[1][1]+atEa)+(deltaMio)*(Eaa[1][1]+aaE[1][1])
          }
        };
      
      // double deltaMio1=(youngL-youngT)/2;
      // double deltaS[2][2]; // based on  ... 
      // deltaS[0][0]=2*(deltaMio1)* atEa * directAniso[0][0];
      // deltaS[1][0]=2*(deltaMio1)* atEa * directAniso[1][0];
      // deltaS[0][1]=2*(deltaMio1)* atEa * directAniso[0][1];
      // deltaS[1][1]=2*(deltaMio1)* atEa * directAniso[1][1];




      strainZ +=restingArea*(1-poissonT*((2*lambdaT*trE+2*mioT*trE)+deltaS[0][0]+deltaS[1][1])/youngT);
      //std::cerr<< "cell "<< cellIndex<< " thickness :  " << strainZ << std::endl;
      
      //<<<<<<<<<<<<<<<<<<<isotropic force from stress tensor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      // double ss[2][2];//lambda(trE)I+2mioE
      // ss[0][0]=lambdaT*trE+2*mioT*Egreen[0][0];
      // ss[0][1]=            2*mioT*Egreen[0][1];
      // ss[1][0]=            2*mioT*Egreen[1][0];
      // ss[1][1]=lambdaT*trE+2*mioT*Egreen[1][1];

       

      // double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      // TPK[0][0]=restingArea*(DeformGrad[0][0]*ss[0][0]+DeformGrad[0][1]*ss[1][0]);
      // TPK[1][0]=restingArea*(DeformGrad[1][0]*ss[0][0]+DeformGrad[1][1]*ss[1][0]);
      // TPK[0][1]=restingArea*(DeformGrad[0][0]*ss[0][1]+DeformGrad[0][1]*ss[1][1]);
      // TPK[1][1]=restingArea*(DeformGrad[1][0]*ss[0][1]+DeformGrad[1][1]*ss[1][1]);

      // //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      // //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      // double deltaFTPKlocal[2][2];
      // deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      // deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      // double deltaFTPK[2][2]; 
      // deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      TPK[0][0]=restingArea*(DeformGrad[0][0]*deltaS[0][0]+DeformGrad[0][1]*deltaS[1][0]);
      TPK[1][0]=restingArea*(DeformGrad[1][0]*deltaS[0][0]+DeformGrad[1][1]*deltaS[1][0]);
      TPK[0][1]=restingArea*(DeformGrad[0][0]*deltaS[0][1]+DeformGrad[0][1]*deltaS[1][1]);
      TPK[1][1]=restingArea*(DeformGrad[1][0]*deltaS[0][1]+DeformGrad[1][1]*deltaS[1][1]);

      //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      double deltaFTPKlocal[3][2];
      deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      deltaFTPKlocal[1][0]= TPK[0][0]*ShapeVectorResting[1][0]+TPK[0][1]*ShapeVectorResting[1][1];
      deltaFTPKlocal[1][1]= TPK[1][0]*ShapeVectorResting[1][0]+TPK[1][1]*ShapeVectorResting[1][1];

      deltaFTPKlocal[2][0]= TPK[0][0]*ShapeVectorResting[2][0]+TPK[0][1]*ShapeVectorResting[2][1];
      deltaFTPKlocal[2][1]= TPK[1][0]*ShapeVectorResting[2][0]+TPK[1][1]*ShapeVectorResting[2][1];


      double deltaFTPK[3][3]; 
      deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];

      deltaFTPK[1][0]= rotation[0][0]*deltaFTPKlocal[1][0]+rotation[0][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][1]= rotation[1][0]*deltaFTPKlocal[1][0]+rotation[1][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][2]= rotation[2][0]*deltaFTPKlocal[1][0]+rotation[2][1]*deltaFTPKlocal[1][1];
      
      deltaFTPK[2][0]= rotation[0][0]*deltaFTPKlocal[2][0]+rotation[0][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][1]= rotation[1][0]*deltaFTPKlocal[2][0]+rotation[1][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][2]= rotation[2][0]*deltaFTPKlocal[2][0]+rotation[2][1]*deltaFTPKlocal[2][1];


      // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



      


      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
      //......................Rstrain here........ ( or clockwise ordering of nodes/edges)
          

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  


      // for (int r=0 ; r<2 ; r++) 
      //   for (int s=0 ; s<2 ; s++) 
      //     std::cerr<<" almansi  "<<StrainAlmansi[r][s]<<std::endl;      


   
    //---- Anisotropic Correction Force-------------------------------
       double deltaF[3][3];




        for ( int i=0 ; i<3 ; ++i )  // from stress tensor(equipartitioning energy)
          for ( int j=0 ; j<3 ; ++j )
            deltaF[i][j]=(-deltaFTPK[i][j]);
      
        double  I1=trE;

        // energy
        EnergyIso +=( (lambdaT/2)*I1*I1 + mioT*I2 )*restingArea; //<<<<<<<<<<<<<<<<<<<<<
      
        EnergyAniso +=( (deltaLam/2)*I4*I1 + deltaMio*I5 )*restingArea; //<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        
        
        //Forces of vertices   
        double Force[3][3]={{0,0,0},{0,0,0},{0,0,0}};                                          
        


        Force[0][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][0]-position[0][0])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][0]-position[0][0])
          + deltaF[0][0]; 
        Force[0][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][1]-position[0][1])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][1]-position[0][1])
          + deltaF[0][1];  
        Force[0][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][2]-position[0][2])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][2]-position[0][2])
          + deltaF[0][2]; 

        Force[1][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][0]-position[1][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][0]-position[1][0])
          + deltaF[1][0];  
        Force[1][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][1]-position[1][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][1]-position[1][1])
          + deltaF[1][1];  
        Force[1][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][2]-position[1][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][2]-position[1][2])
          + deltaF[1][2];  
        
        Force[2][0]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][0]-position[2][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][0]-position[2][0])
          + deltaF[2][0];  
        Force[2][1]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][1]-position[2][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][1]-position[2][1])
          + deltaF[2][1];  
        Force[2][2]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][2]-position[2][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][2]-position[2][2])
          + deltaF[2][2];
        
        
        bool isSliver=false;
        if (std::abs(1-cosAngle[0])<0.0001 || std::abs(1-cosAngle[1])<0.0001 ||  std::abs(1-cosAngle[2])<0.0001) {
          std::cerr<<" sliver in "<<cellIndex;
          double tmp1=cellDerivs[cellIndex][comIndex  ]+vertexDerivs[v2][0]+vertexDerivs[v3][0],
            tmp2=cellDerivs[cellIndex][comIndex+1]+vertexDerivs[v2][1]+vertexDerivs[v3][1],
            tmp3=cellDerivs[cellIndex][comIndex+2]+vertexDerivs[v2][2]+vertexDerivs[v3][2];

          cellDerivs[cellIndex][comIndex  ] +=tmp1;
          cellDerivs[cellIndex][comIndex+1] +=tmp2;
          cellDerivs[cellIndex][comIndex+2] +=tmp3;
          
          vertexDerivs[v2][0] += tmp1;
          vertexDerivs[v2][1] += tmp2;
          vertexDerivs[v2][2] += tmp3;
          
          vertexDerivs[v3][0] += tmp1;
          vertexDerivs[v3][1] += tmp2;
          vertexDerivs[v3][2] += tmp3;

          isSliver=true;
          std::cerr<<" there is a sliver!!!"<< " in cell "<< cellIndex<<" and wall  "<<wallindex<<std::endl;
        }  
        


        
        // adding TRBSMT forces to the total vertexDerives
        if (! isSliver) {
          cellDerivs[cellIndex][comIndex  ] += Force[0][0];
          cellDerivs[cellIndex][comIndex+1] += Force[0][1];
          cellDerivs[cellIndex][comIndex+2] += Force[0][2];
          
          vertexDerivs[v2][0] += Force[1][0];
          vertexDerivs[v2][1] += Force[1][1];
          vertexDerivs[v2][2] += Force[1][2];
          
          vertexDerivs[v3][0] += Force[2][0];
          vertexDerivs[v3][1] += Force[2][1];
          vertexDerivs[v3][2] += Force[2][2];
        }
     
              
    }
    
  
      //cpuTime2=clock();
      //tRest+=cpuTime2-cpuTime1;


      //cpuTime2=clock();
      //tRest+=cpuTime2-cpuTime1;
   
 
  }
  
  //cpuTimef=clock();
  //tTotal+=cpuTimef-cpuTime0;
    
}



void VertexFromTRBScenterTriangulationMT::
update(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData, 
       double h) {
  
  timeC+=h;

  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t numVertices = T.numVertex();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t comIndex = variableIndex(1,0);
  size_t lengthInternalIndex = comIndex+dimension;

  double neighborweight=parameter(5);
  double strainEpcilon =0.000001;
  double stressEpcilon =0.000001;    
  // double TotalVolume=0;
  // double deltaVolume=0;
  // for(size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex){ // stimating volume for equilibrium 
  //   TotalVolume +=std::sqrt(vertexData[vertexIndex][0]*vertexData[vertexIndex][0] +
  //                           vertexData[vertexIndex][1]*vertexData[vertexIndex][1] +
  //                           vertexData[vertexIndex][2]*vertexData[vertexIndex][2] );
  // }
  // deltaVolume=TotalVolume-cellData[0][25];
  // cellData[0][25]=TotalVolume;
  // cellData[0][24]=deltaVolume;
 
  size_t MTindex           =variableIndex(0,1);	 
  size_t strainAnIndex     =variableIndex(0,2);	
  size_t stressAnIndex     =variableIndex(0,3);	
  size_t areaRatioIndex    =variableIndex(0,4);	
  size_t isoEnergyIndex    =variableIndex(0,5);	
  size_t anisoEnergyIndex  =variableIndex(0,6);	
  size_t youngLIndex       =variableIndex(0,7);	
  size_t MisesStressIndex  =variableIndex(0,8);	
  size_t stressTensorIndex =variableIndex(0,9);	
  size_t normalVectorIndex =variableIndex(0,10);

  double youngMatrix= parameter(0);    
  double youngFiber = parameter(1); 
  double poissonL   = parameter(2);    
  double poissonT   = parameter(3);
  double TETA       = parameter(8);  
  
 

  //std::cerr<<" here is update "<<std::endl;
  //std::cout<<"begin:"<<std::endl;
  for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    size_t numWalls = T.cell(cellIndex).numWall();
    
  

    if(  T.cell(cellIndex).numVertex()!= numWalls ) {
     
      std::cerr << "VertexFromTRBScenterTriangulationMT::update() same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
  
 
   
 double youngL=1;
    double youngT=1;
    


    if( parameter(4)==1){  // material anisotropy via FiberModel
      youngL = cellData[cellIndex][youngLIndex]; 
      youngT = 2*youngMatrix+youngFiber-youngL; 
    }
    
    if( parameter(4)==0 ){ // constant anisotropic material
      youngL = youngMatrix+youngFiber;
      youngT = youngMatrix; 
    }

        
    if( parameter(4)==-1 ){ //hypocotyl Siobhan
      youngL = 50;
      youngT = 50; 
      poissonL=0.2;
      poissonT=0.2;

      if(cellData[cellIndex][37]==-1){ // epidermis
	youngL *=parameter(0);
	youngT *=parameter(0);
      }
      if(cellData[cellIndex][37]==-2){ // inner
	youngL *=parameter(1); //1, 6 // hoop
	youngT *=1;                   // axial
      }
      if(cellData[cellIndex][37]==-3){  // anti. axial
	youngL *=parameter(2); //0.1, 0.5, 1
	youngT *=parameter(2); //0.1, 0.5, 1
      }
      if(cellData[cellIndex][37]==-4){  // anti transverse
	youngL *=parameter(3); //1, 10
	youngT *=parameter(3); //1, 10
      }
    }
    

      
      
    if( parameter(4)==2){  // for varrying material anisotropy with constant overall stiffness for energy landscape
      youngL =youngFiber; 
      youngT =youngMatrix-youngL;  // here youngMatrix is total stiffness
    }
    if( parameter(4)==3){  // for varrying material anisotropy with constant overall stiffness for energy landscape       
      double totalElast=youngMatrix;
      double Maniso=youngFiber;
      youngT =((1-Maniso)/(2-Maniso))*totalElast; 
      youngL=totalElast-youngT;
      
      //youngL =youngMatrix+youngFiber; 
      // youngT =youngMatrix+std::sqrt(
      //                               ((youngFiber-totalElast)*4*totalElast)
      //                               /
      //                               ((2-3.1415)*3.1415)
      //                               );  // here youngMatrix is total stiffness
      //youngT =youngMatrix+totalElast*std::sqrt(
      //                                         -(std::log((youngFiber/totalElast)-0.55))
      //                                         /6
      //                                         );
      //youngT =youngMatrix+2*(totalElast-youngFiber)/(3.1415-2);
      
      // youngL=youngMatrix+youngFiber; 
      // youngT =youngMatrix+(youngFiber-slope*totalElast)/(1-2*slope); 
    }
    
    // if( parameter(4)==5){  // for domain dependent material properties - Hypocotyl bending 3d (internal tissue softer)    
      //   if (cellData[cellIndex][38]==-1){
      //     youngL =youngMatrix+youngFiber; 
      //     youngT =youngMatrix;  
          
      //   }
      //   else{
      //     youngL =0.2*(youngMatrix+youngFiber); 
      //     youngT =0.2*youngMatrix;  
      //   }

      //   if (cellData[cellIndex][5]==5){
      //     youngL =0.01*(youngMatrix+youngFiber); 
      //     youngT =0.01*youngMatrix;  
      //   }
      // }

      if( parameter(4)==5){// for pavement cell resolution addaptive stiffness
        youngFiber=cellData[cellIndex][anisoEnergyIndex];
        double fiberL=cellData[cellIndex][youngLIndex];
        youngL = youngMatrix+fiberL;
        youngT = youngMatrix+youngFiber-fiberL;
       
        // if(cellData[cellIndex][30]==1){
        //   youngL = 20;
        //   youngT = 20;
        
        // }
      }
      
      if( parameter(4)==6){  // material anisotropy via FiberModel and loosening adhoc based on auxin
        youngL = cellData[cellIndex][youngLIndex]; 
        youngT = 2*youngMatrix+youngFiber-youngL; 
        
        double Kconc=0.005;
        double Nc1=2;
        double conc=cellData[cellIndex][13];
        double frac=(1-std::pow(conc,Nc1)/(std::pow(Kconc,Nc1)+std::pow(conc,Nc1)));
        youngL*=0.25+0.75*frac;
        youngT*=0.25+0.75*frac;
        //if(cellIndex==0) std::cerr<<cellData[cellIndex][13]<<"  "<<youngL<<"  " <<youngT<<std::endl;
      }
      
      if( parameter(4)==7){  // material anisotropy via FiberModel and destroying fibers based on auxin
	youngL = cellData[cellIndex][youngLIndex]-youngMatrix;
	youngT = 2*youngMatrix+youngFiber-youngL-2*youngMatrix;  
        
        double Kconc=0.005;
        double Nc1=2;
        double conc=cellData[cellIndex][13];
        double frac=(1-std::pow(conc,Nc1)/(std::pow(Kconc,Nc1)+std::pow(conc,Nc1)));
        youngL*=frac;
        youngT*=frac;
        youngL+=youngMatrix;
        youngT+=youngMatrix;
        //cellData[cellIndex][youngLIndex]=youngL; 
	
      }
      if( parameter(4)==8){  // material anisotropy via FiberModel and destroying fibers and a fraction of matrix based on auxin
	youngL = cellData[cellIndex][youngLIndex]-youngMatrix;
	youngT = 2*youngMatrix+youngFiber-youngL-2*youngMatrix;  
        
        double Kconc=0.005;
        double Nc1=2;
        double conc=cellData[cellIndex][13];
        double frac=(1-std::pow(conc,Nc1)/(std::pow(Kconc,Nc1)+std::pow(conc,Nc1)));
        youngL*=frac;
        youngT*=frac;
        youngL+=youngMatrix*(0.2+0.8*frac);
        youngT+=youngMatrix*(0.2+0.8*frac);
        //cellData[cellIndex][youngLIndex]=youngL; 
	
      }
      

                     
      // If( parameter(4)==-2){  // for heterogeneous stiffness(adhoc)
      //   double hFactor=0;
      //   double Hthreshold=0.01;
      //   if (std::fabs(cellData[cellIndex][youngLIndex])>Hthreshold){
      //     hFactor=Hthreshold;
      //   }
      //   else{
      //     hFactor=std::fabs(cellData[cellIndex][youngLIndex]); // take the heterogeneous info from this index
      //   }
      //   hFactor *=std::fabs(parameter(4)); // factor for heterogeneity      
      
      //   youngL = youngMatrix+youngFiber-hFactor*(youngMatrix+youngFiber); 
      //   youngT = youngMatrix-hFactor*(youngMatrix); 
      // }
      
      // if (cellData[cellIndex][31]==10){
      //   youngL = 150;
      //   youngT = 10;
      // }
     
      
      // if( parameter(4)<0){  // for spatial elasticity (ad-hoc)
      //   double x0=90;
      //   double y0=120;
        
      //   double xx=cellData[cellIndex][comIndex]-x0;
      //   double yy=cellData[cellIndex][comIndex+1]-y0;
        
      //   youngL = youngMatrix+youngFiber*(1-std::exp(parameter(4)*(xx*xx+yy*yy)));
      //   youngT = youngMatrix; 
      //   cellData[cellIndex][13]=youngL;
      // }
      
      // ad-hoc for regional loosening
      // if ( std::sqrt(cellData[cellIndex][comIndex  ]*cellData[cellIndex][comIndex  ]
      //                +cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])<30)
      //   cellData[cellIndex][youngLIndex]=20;
      // else
      //   cellData[cellIndex][youngLIndex]=200;
      
      
      
      // ad-hoc for 3d marcus
      
      // if(cellData[cellIndex][25]==-3){  // sidewalls (anticlinals)
      //   youngL = youngMatrix;
      //   youngT = youngMatrix; 
      //   if(cellData[cellIndex][29]==1){
      
      //   }
      // }
      
      //youngL = (1-cellData[cellIndex][30]/0.6)*youngL;
      //youngT = (1-cellData[cellIndex][30]/0.6)*youngT;
      
      // double sEpi=10;
      
      // if(cellData[cellIndex][27]==-2){ //L1 epidermis
      //   youngL = sEpi*youngL;
      //   youngT = sEpi*youngT; 
      // }
      // if(cellData[cellIndex][30]==10){ // L1 or L2 patch
      //   youngL = 0.1*youngL;
      //   youngT = 0.1*youngT; 
      // }
      
      // if(cellData[cellIndex][30]==5){ // b cells
      //   youngL = 0.5*youngL;
      //   youngT = 0.5*youngT; 
      // }
      // if(cellData[cellIndex][30]==2){ // c cells
      //   youngL = 0.8*youngL;
      //   youngT = 0.8*youngT; 
      // }
      
      
      // if(cellData[cellIndex][25]==-3 && cellData[cellIndex][29]==1){  // L1 anticlinals
      //   youngL = sEpi*youngL;
      //   youngT = sEpi*youngT; 
      // }
      
      // if(cellData[cellIndex][30]==400){
      //   youngL = 0.1*youngL;
      //   youngT = 0.1*youngT;
      // }
      
      
      // ad-hoc for root hair cell files and ....
      // double polar=180*std::acos(cellData[cellIndex][comIndex]/
      //                            std::sqrt(cellData[cellIndex][comIndex  ]*cellData[cellIndex][comIndex  ]+
      //                                      cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])
      //                            );
      // if (std::fmod(std::floor(polar/12),2)==1){
      //   youngL =110 ;
      //   youngT =110 ; 
      
      // }
      
      // if (std::sqrt(cellData[cellIndex][comIndex]*cellData[cellIndex][comIndex]+
      //               (cellData[cellIndex][comIndex+1]-60)*(cellData[cellIndex][comIndex+1]-60))<35){
      //   youngL =5 ;
      //   youngT =5 ; 
      
      // }
      
      // if (std::sqrt(cellData[cellIndex][comIndex]*cellData[cellIndex][comIndex]+
      //                         cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])<10){
      //   youngL =100 ;
      //   youngT =100 ; 
      
      // }
      

  

    double lambdaL, mioL, lambdaT, mioT; // ,lambdaTmatrix, mioTmatrix;
    
    if (parameter(7)==0){      
      // Lame coefficients based on plane strain (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/((1+poissonL)*(1-2*poissonL));
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/((1+poissonT)*(1-2*poissonT));
      mioT=youngT/(2*(1+poissonT));

   
    } 
    else{      
      // Lame coefficients based on plane stress (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/(1-poissonL*poissonL);
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/(1-poissonT*poissonT);
      mioT=youngT/(2*(1+poissonT));
     
      // lambdaTmatrix=youngMatrixA*poissonT/(1-poissonT*poissonT);
      // mioTmatrix=youngMatrixA/(2*(1+poissonT));
    }
    

    //std::cerr<<"youngL  "<<youngL<<", youngT"<<youngT<<std::endl;

    // Lame coefficients based on delin. paper (for 2D 0<poisson<1)
    // double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
    // double mioL=youngL/(1+poissonL);
    // double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
    // double mioT=youngT/(1+poissonT);
    

    

    double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double StressCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double normalGlob[3]={0,0,0};
    double TotalCellRestingArea=0;
    double TotalCellArea=0;
    double EnergyIso=0;     
    //double EnergyIsoFiber=0;            
    double EnergyAniso=0;
    double strainZ=0;

    if ( parameter(9)==1 ) {  // aniso direction from TETA 
      
      cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);    
      cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
      cellData[cellIndex][variableIndex(0,1)+2]=0;
      // if(parameter(10)!=100){
      //   if(cellIndex==0 || cellIndex==2 ){
      //     cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);    
      //     cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
      //     cellData[cellIndex][variableIndex(0,1)+2]=0;
      //   }
      //   if(cellIndex==1 || cellIndex==3 ){
      //     cellData[cellIndex][variableIndex(0,1)]=std::cos(parameter(10));    
      //     cellData[cellIndex][variableIndex(0,1)+1]=std::sin(parameter(10));
      //     cellData[cellIndex][variableIndex(0,1)+2]=0;
      //   }
      // }
    }
    
  

    
  
      // One triangle per 'vertex' in cyclic order
      for (size_t wallindex=0; wallindex<numWalls; ++wallindex) {
        
      size_t kPlusOneMod = (wallindex+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
      size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //size_t w1 = internal wallindex
      size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //size_t w3 = internal wallindex+1

      // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      //position[0][2] z for vertex 1 of the current element
      

      // Resting lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)

      // size_t restAddIndex;
      // if (T.wall(w2).cell1()->index()==cellIndex)
      //   restAddIndex=1;
      // else if (T.wall(w2).cell2()->index()==cellIndex)
      //   restAddIndex=2;
      // else{
      //   std::cerr<<"something bad happened"<<std::endl;
      //   exit(-1);
      // }

      double restingLength[3];
      if(parameter(10)==1) { // double resting length  
        restingLength[0] =cellData[cellIndex][lengthInternalIndex + 2*wallindex+1];
        //restingLength[1] = wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex];
        restingLength[1] =wallData[w2][wallLengthIndex]+
          cellData[cellIndex][lengthInternalIndex+2*numWalls+wallindex];
        restingLength[2] =cellData[cellIndex][lengthInternalIndex + 2*kPlusOneMod];
        //cellData[cellIndex][31]=cellIndex;        
        // if(cellIndex==408 && wallindex==0) 
        //   std::cout<<restingLength[0]<<" "<<restingLength[1]<<" "<<restingLength[2]<<std::endl;
      }      
      else{                 // single resting length
        restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
        restingLength[1] = wallData[w2][wallLengthIndex];
        restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
      }
      
      // Lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
      double length[3]=
        { std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) ),
          
          T.wall(w2).lengthFromVertexPosition(vertexData),
          
          std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) )
        };
      
      
      // Area of the element (using Heron's formula)                                      
      double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                    (-restingLength[0]+restingLength[1]+restingLength[2])*
                                    ( restingLength[0]-restingLength[1]+restingLength[2])*
                                    ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
      
      //double currentArea=std::sqrt( ( length[0]+length[1]+length[2])*
      //                              (-length[0]+length[1]+length[2])*
      //                              ( length[0]-length[1]+length[2])*
      //                              ( length[0]+length[1]-length[2])  )*0.25;

 
     
   

      
      double temp = 1.0/(restingArea*16);                                      
      
      //Area of the element (using Heron's formula)                                      
      double Area=std::sqrt( ( length[0]+length[1]+length[2])*
                             (-length[0]+length[1]+length[2])*
                             ( length[0]-length[1]+length[2])*
                             ( length[0]+length[1]-length[2])  )*0.25;
      
      
      
      // calculating the angles between shape vectors and anisotropy direction in resting shape when anisotropy vector is provided in current shape
      
      //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
      // double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
      //                                  (length[0]*length[1]*2)    );
      
      // double Qa=std::cos(CurrentAngle1)*length[0];
      // double Qc=std::sin(CurrentAngle1)*length[0];
      // double Qb=length[1];
      
      double CurrentAngle1=(length[0]*length[0]
                            +length[1]*length[1]
                            -length[2]*length[2])/
        (length[0]*length[1]*2);
      
      double Qa=(CurrentAngle1)*length[0];
      double Qc=std::sqrt(1-CurrentAngle1*CurrentAngle1)*length[0];
      double Qb=length[1];
      
      // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
      
      
      double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
                                        {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
                                        { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
    
      
      // Local coordinates of the resting shape ( counterclockwise )
      // double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]
      //                                   +restingLength[1]*restingLength[1]
      //                                   -restingLength[2]*restingLength[2])
      //                                  /(restingLength[0]*restingLength[1]*2)    );
      
      
      // double Pa=std::cos(RestingAngle1)*restingLength[0];
      // double Pc=std::sin(RestingAngle1)*restingLength[0];
      // double Pb=restingLength[1];
      
      
      double RestingAngle1=  (restingLength[0]*restingLength[0]
                              +restingLength[1]*restingLength[1]
                              -restingLength[2]*restingLength[2])/
        (restingLength[0]*restingLength[1]*2);
      
      
      double Pa=(RestingAngle1)*restingLength[0];
      double Pc=std::sqrt(1-RestingAngle1*RestingAngle1)*restingLength[0];
      double Pb=restingLength[1];
      
      // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                        {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                        { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
      
      
      
      // //Strain tensor  (clockwise ordering of nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );
      
      // double Qa=std::cos(CurrentAngle2)*length[2];
      // double Qb=length[1];
      // double Qc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVectorCurrent[3][2]={ {  0   ,       1/Qc      }, 
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) },       
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) }  };
      
      
  
      
      
      // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
      
      
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
  
            
   
      // Aniso vector in current shape in global coordinate system
      double  AnisoCurrGlob[3]=
        { cellData[cellIndex][variableIndex(0,1)],  
          cellData[cellIndex][variableIndex(0,1)+1],
          cellData[cellIndex][variableIndex(0,1)+2]
        };
      
      // rotating the anisotropy vector from global coordinate system to the local one in the current shape
      double AnisoCurrLocal[3]=
        {  rotation[0][0]*AnisoCurrGlob[0]+
           rotation[1][0]*AnisoCurrGlob[1]+
           rotation[2][0]*AnisoCurrGlob[2],
           
           rotation[0][1]*AnisoCurrGlob[0]+
           rotation[1][1]*AnisoCurrGlob[1]+
           rotation[2][1]*AnisoCurrGlob[2],
           
           rotation[0][2]*AnisoCurrGlob[0]+
           rotation[1][2]*AnisoCurrGlob[1]+
           rotation[2][2]*AnisoCurrGlob[2]
        };
      
      //Anisotropic Correction is based on difference between Lam Coefficients of Longitudinal and Transverse dirrections:
      // double deltaLam=AnisoMeasure*(lambdaL-lambdaT);
      // double deltaMio=AnisoMeasure*(mioL-mioT);
      double deltaLam=lambdaL-lambdaT;
      double deltaMio=mioL-mioT;  

      // double deltaLamIsoFiber=lambdaT-lambdaTmatrix;
      // double deltaMioIsoFiber=mioT-mioTmatrix;  
 

    
     
     
     
    
       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
      // strain tensor in resting shape E=0.5(FtF-I)
      // trE
      // B=FFt
      // axa (direct product of aniso vector in resting shape)
      // atEa
      // E(axa) and (axa)E

      
     

     
      
      double positionLocal[3][2]={ {Qa , Qc}, 
                                   {0  , 0 },  
                                   {Qb , 0 }  };
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F=Sigma i Qi x Di   <---------------------------------------
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];        
      } 
    


      // transforming aniso vector in local coordinate from current to resting shape using inverse of deformation gradient
      double AnisoRestLocal[2]={ DeformGrad[1][1] *AnisoCurrLocal[0]-DeformGrad[0][1]*AnisoCurrLocal[1],
                                 -DeformGrad[1][0]*AnisoCurrLocal[0]+DeformGrad[0][0]*AnisoCurrLocal[1]};
                
      double AnisoMeasure=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]);
      // std::cerr<< "cell "<< cellIndex<<" AnisoMeasure "<< AnisoMeasure<<std::endl;
      
      // choosing a random normalized dirrection for anisotropy if AnisoVector is close to perpendicular to the cell plane
      if ( AnisoMeasure<0.0001) {
        double randomAngle=((rand()%360)*2*3.14159265)/360;
        
        AnisoRestLocal[0]=std::cos(randomAngle);
        AnisoRestLocal[1]=std::sin(randomAngle);
      }  
      else {// normalizing the anisoVector if it is not random
        AnisoRestLocal[0]=AnisoRestLocal[0]/AnisoMeasure;
        AnisoRestLocal[1]=AnisoRestLocal[1]/AnisoMeasure;        
      }
      
      double directAniso[2][2]={{AnisoRestLocal[0]*AnisoRestLocal[0],AnisoRestLocal[0]*AnisoRestLocal[1]},
                                {AnisoRestLocal[1]*AnisoRestLocal[0],AnisoRestLocal[1]*AnisoRestLocal[1]}};



      double LeftCauchy[2][2]= // B=FFt
        { { DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1],
            DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1]
          },
          { DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1],
            DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1]
          }
        };
      
      double Egreen[2][2]= //E=0.5(C-I)
        { { 0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1),
            0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1])   
          },
          { 0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]), 
            0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1)
          } 
        };
      //std::cout<<std::abs(Egreen[0][0]+Egreen[1][1])/2<<" "<<std::abs(Egreen[0][1])<<std::endl;

      double E2[2][2]= // used for energy calculation only
        { { Egreen[0][0]*Egreen[0][0]+Egreen[0][1]*Egreen[1][0],
            Egreen[0][0]*Egreen[0][1]+Egreen[0][1]*Egreen[1][1]
          },
          { Egreen[1][0]*Egreen[0][0]+Egreen[1][1]*Egreen[1][0],
            Egreen[1][0]*Egreen[0][1]+Egreen[1][1]*Egreen[1][1]
          }
        };

      double trE=Egreen[0][0]+Egreen[1][1];

      double I2=E2[0][0]+E2[1][1]; //trE2 used for energy calculation only

      
      double I5=AnisoRestLocal[0]*AnisoRestLocal[0]*E2[0][0]   //atE2a used for energy calculation only
        +AnisoRestLocal[0]*AnisoRestLocal[1]*(E2[0][1]+E2[1][0])
        +AnisoRestLocal[1]*AnisoRestLocal[1]*E2[1][1];
      
    
     
      temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
      double StrainAlmansi[2][2]= // e=0.5(1-B^-1)  True strain tensor
        { { 0.5*(1-(LeftCauchy[1][1]/temp)) , 0.5*LeftCauchy[0][1]/temp },
          { 0.5*LeftCauchy[1][0]/temp       , 0.5*(1-(LeftCauchy[0][0]/temp))}
        };
      
      if(StrainAlmansi[0][0] != StrainAlmansi[0][0] ||
         StrainAlmansi[1][1] != StrainAlmansi[1][1] ||
         StrainAlmansi[0][1] != StrainAlmansi[0][1] ||
         StrainAlmansi[1][0] != StrainAlmansi[1][0] ) 
        std::cerr<<std::endl<<" strain is wrong "<<StrainAlmansi[0][0]<<StrainAlmansi[0][0]
                 <<StrainAlmansi[0][0] << StrainAlmansi[0][0] <<"   Q "<<restingLength[0]<<" "<<restingLength[1]<<" "<<restingLength[2]<<"    P "<<Pa<<" "<<Pb<<" "<<Pc<<std::endl;
      double atEa=AnisoRestLocal[0]*AnisoRestLocal[0]*Egreen[0][0]
        +AnisoRestLocal[0]*AnisoRestLocal[1]*(Egreen[0][1]+Egreen[1][0])
        +AnisoRestLocal[1]*AnisoRestLocal[1]*Egreen[1][1];

      double I4=atEa;

      
    

  
      double Eaa[2][2]=
        { { Egreen[0][0]*directAniso[0][0]+Egreen[0][1]*directAniso[1][0],
            Egreen[0][0]*directAniso[0][1]+Egreen[0][1]*directAniso[1][1]
          },
          { Eaa[1][0]= Egreen[1][0]*directAniso[0][0]+Egreen[1][1]*directAniso[1][0],             
            Eaa[1][1]= Egreen[1][0]*directAniso[0][1]+Egreen[1][1]*directAniso[1][1]
          }
        };        

      double aaE[2][2]=
        { { directAniso[0][0]*Egreen[0][0]+directAniso[0][1]*Egreen[1][0],        
            directAniso[0][0]*Egreen[0][1]+directAniso[0][1]*Egreen[1][1]
          },        
          { directAniso[1][0]*Egreen[0][0]+directAniso[1][1]*Egreen[1][0],              
            directAniso[1][0]*Egreen[0][1]+directAniso[1][1]*Egreen[1][1]        
          }
        };

      double B2[2][2]=// LeftCauchy^2
        { { LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0],
            LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1]
          },
          { LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0],
            LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1]
          }
        };
      
      double areaFactor=restingArea/Area; // 1/detF
      //double areaFactor=Area/restingArea; // detF
      
      double Sigma[2][2]= // true stress tensor (isotropic term) 
        { { areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][0]+(mioT)*B2[0][0]),
            areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][1]+(mioT)*B2[0][1])
          },      
          { areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][0]+(mioT)*B2[1][0]),      
            areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][1]+(mioT)*B2[1][1])
          }
        };
      

      // ENERGY MODELS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


      // double deltaS[2][2]; // based on Delin. paper
      // deltaS[0][0]=deltaLam*(trE*directAniso[0][0]+atEa)+(2*deltaMio)*(Eaa[0][0]+aaE[0][0])-(deltaLam+2*deltaMio)*atEa*directAniso[0][0];
      // deltaS[1][0]=deltaLam*(trE*directAniso[1][0]     )+(2*deltaMio)*(Eaa[1][0]+aaE[1][0])-(deltaLam+2*deltaMio)*atEa*directAniso[1][0];
      // deltaS[0][1]=deltaLam*(trE*directAniso[0][1]     )+(2*deltaMio)*(Eaa[0][1]+aaE[0][1])-(deltaLam+2*deltaMio)*atEa*directAniso[0][1];
      // deltaS[1][1]=deltaLam*(trE*directAniso[1][1]+atEa)+(2*deltaMio)*(Eaa[1][1]+aaE[1][1])-(deltaLam+2*deltaMio)*atEa*directAniso[1][1];

      double deltaS[2][2]= // based on  equipartitioning
        { { (deltaLam/2)*(trE*directAniso[0][0]+atEa)+(deltaMio)*(Eaa[0][0]+aaE[0][0]),
            (deltaLam/2)*(trE*directAniso[0][1]     )+(deltaMio)*(Eaa[0][1]+aaE[0][1])
          },
          { (deltaLam/2)*(trE*directAniso[1][0]     )+(deltaMio)*(Eaa[1][0]+aaE[1][0]),      
            (deltaLam/2)*(trE*directAniso[1][1]+atEa)+(deltaMio)*(Eaa[1][1]+aaE[1][1])
          }
        };
      
      // double deltaMio1=(youngL-youngT)/2;
      // double deltaS[2][2]; // based on  ... 
      // deltaS[0][0]=2*(deltaMio1)* atEa * directAniso[0][0];
      // deltaS[1][0]=2*(deltaMio1)* atEa * directAniso[1][0];
      // deltaS[0][1]=2*(deltaMio1)* atEa * directAniso[0][1];
      // deltaS[1][1]=2*(deltaMio1)* atEa * directAniso[1][1];




      strainZ +=restingArea*(1-poissonT*((2*lambdaT*trE+2*mioT*trE)+deltaS[0][0]+deltaS[1][1])/youngT);
      //std::cerr<< "cell "<< cellIndex<< " thickness :  " << strainZ << std::endl;
      
      //<<<<<<<<<<<<<<<<<<<isotropic force from stress tensor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      // double ss[2][2];//lambda(trE)I+2mioE
      // ss[0][0]=lambdaT*trE+2*mioT*Egreen[0][0];
      // ss[0][1]=            2*mioT*Egreen[0][1];
      // ss[1][0]=            2*mioT*Egreen[1][0];
      // ss[1][1]=lambdaT*trE+2*mioT*Egreen[1][1];

       

      // double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      // TPK[0][0]=restingArea*(DeformGrad[0][0]*ss[0][0]+DeformGrad[0][1]*ss[1][0]);
      // TPK[1][0]=restingArea*(DeformGrad[1][0]*ss[0][0]+DeformGrad[1][1]*ss[1][0]);
      // TPK[0][1]=restingArea*(DeformGrad[0][0]*ss[0][1]+DeformGrad[0][1]*ss[1][1]);
      // TPK[1][1]=restingArea*(DeformGrad[1][0]*ss[0][1]+DeformGrad[1][1]*ss[1][1]);

      // //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      // //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      // double deltaFTPKlocal[2][2];
      // deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      // deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      // double deltaFTPK[2][2]; 
      // deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      TPK[0][0]=restingArea*(DeformGrad[0][0]*deltaS[0][0]+DeformGrad[0][1]*deltaS[1][0]);
      TPK[1][0]=restingArea*(DeformGrad[1][0]*deltaS[0][0]+DeformGrad[1][1]*deltaS[1][0]);
      TPK[0][1]=restingArea*(DeformGrad[0][0]*deltaS[0][1]+DeformGrad[0][1]*deltaS[1][1]);
      TPK[1][1]=restingArea*(DeformGrad[1][0]*deltaS[0][1]+DeformGrad[1][1]*deltaS[1][1]);

      //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      double deltaFTPKlocal[3][2];
      deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      deltaFTPKlocal[1][0]= TPK[0][0]*ShapeVectorResting[1][0]+TPK[0][1]*ShapeVectorResting[1][1];
      deltaFTPKlocal[1][1]= TPK[1][0]*ShapeVectorResting[1][0]+TPK[1][1]*ShapeVectorResting[1][1];

      deltaFTPKlocal[2][0]= TPK[0][0]*ShapeVectorResting[2][0]+TPK[0][1]*ShapeVectorResting[2][1];
      deltaFTPKlocal[2][1]= TPK[1][0]*ShapeVectorResting[2][0]+TPK[1][1]*ShapeVectorResting[2][1];


      double deltaFTPK[3][3]; 
      deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];

      deltaFTPK[1][0]= rotation[0][0]*deltaFTPKlocal[1][0]+rotation[0][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][1]= rotation[1][0]*deltaFTPKlocal[1][0]+rotation[1][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][2]= rotation[2][0]*deltaFTPKlocal[1][0]+rotation[2][1]*deltaFTPKlocal[1][1];
      
      deltaFTPK[2][0]= rotation[0][0]*deltaFTPKlocal[2][0]+rotation[0][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][1]= rotation[1][0]*deltaFTPKlocal[2][0]+rotation[1][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][2]= rotation[2][0]*deltaFTPKlocal[2][0]+rotation[2][1]*deltaFTPKlocal[2][1];


      // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



      double deltaSFt[2][2];
      deltaSFt[0][0]=deltaS[0][0]*DeformGrad[0][0]+deltaS[0][1]*DeformGrad[0][1];
      deltaSFt[1][0]=deltaS[1][0]*DeformGrad[0][0]+deltaS[1][1]*DeformGrad[0][1];
      deltaSFt[0][1]=deltaS[0][0]*DeformGrad[1][0]+deltaS[0][1]*DeformGrad[1][1];
      deltaSFt[1][1]=deltaS[1][0]*DeformGrad[1][0]+deltaS[1][1]*DeformGrad[1][1];
      
      double deltaSigma[2][2];// true stress tensor (anisotropic correction term)deltaLambda and deltaMio (Longitudinal-Transverse)
      deltaSigma[0][0]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][0]+DeformGrad[0][1]*deltaSFt[1][0]);
      deltaSigma[1][0]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][0]+DeformGrad[1][1]*deltaSFt[1][0]);
      deltaSigma[0][1]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][1]+DeformGrad[0][1]*deltaSFt[1][1]);
      deltaSigma[1][1]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][1]+DeformGrad[1][1]*deltaSFt[1][1]);

      double sigmafactor=1;
      double deltasigmafactor=1;

      double StressTensor[3][3];
      StressTensor[0][0]=sigmafactor*Sigma[0][0]+deltasigmafactor*deltaSigma[0][0];
      StressTensor[1][0]=sigmafactor*Sigma[1][0]+deltasigmafactor*deltaSigma[1][0];
      StressTensor[0][1]=sigmafactor*Sigma[0][1]+deltasigmafactor*deltaSigma[0][1];
      StressTensor[1][1]=sigmafactor*Sigma[1][1]+deltasigmafactor*deltaSigma[1][1];


      


      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
      //......................Rstrain here........ ( or clockwise ordering of nodes/edges)
          

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  


      // for (int r=0 ; r<2 ; r++) 
      //   for (int s=0 ; s<2 ; s++) 
      //     std::cerr<<" almansi  "<<StrainAlmansi[r][s]<<std::endl;      

      double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
      StrainTensor[0][0]=StrainAlmansi[0][0];
      StrainTensor[1][0]=StrainAlmansi[1][0];
      StrainTensor[0][1]=StrainAlmansi[0][1];
      StrainTensor[1][1]=StrainAlmansi[1][1];


      StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StrainTensor[1][2]=0;
      StrainTensor[2][2]=0;
      StrainTensor[2][0]=0;
      StrainTensor[2][1]=0;
      
      StressTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StressTensor[1][2]=0;
      StressTensor[2][2]=0;
      StressTensor[2][0]=0;
      StressTensor[2][1]=0;

      //rotation matrix to  global coordinate system based on counterclockwise ordering;   rotation[3][3] calculated above  

     
      
      // rotating strain tensor to the global coordinate system
      double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s] += rotation[r][w]*StrainTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StrainTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StrainTensor[r][s] += tempR[r][w]*rotation[s][w]; 
        }

      // rotating stress tensor to the global coordinate system
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          tempR[r][s]=0;
        
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s] += rotation[r][w]*StressTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StressTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StressTensor[r][s] += tempR[r][w]*rotation[s][w]; 
        }
      
 


           
     
      // accumulating strain and stress tensors and normal to cell plane vector to be averaged later
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StrainCellGlobal[r][s] += Area*StrainTensor[r][s];

     

      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StressCellGlobal[r][s] += Area*StressTensor[r][s];
      
      for (int r=0 ; r<3 ; r++) 
        normalGlob[r] += Area*Zcurrent[r];

      TotalCellRestingArea=TotalCellRestingArea+restingArea;
      TotalCellArea=TotalCellArea+Area;
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSORS (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   
    //---- Anisotropic Correction Force-------------------------------
       double deltaF[3][3];




        for ( int i=0 ; i<3 ; ++i )  // from stress tensor(equipartitioning energy)
          for ( int j=0 ; j<3 ; ++j )
            deltaF[i][j]=(-deltaFTPK[i][j]);
      
        double  I1=trE;

        // energy
        EnergyIso +=( (lambdaT/2)*I1*I1 + mioT*I2 )*restingArea; //<<<<<<<<<<<<<<<<<<<<<
      
        EnergyAniso +=( (deltaLam/2)*I4*I1 + deltaMio*I5 )*restingArea; //<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        
        
              
    }
    

     for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellArea; 

        
    
    
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StressCellGlobal[r][s]= StressCellGlobal[r][s]/TotalCellArea; 

    for (int r=0 ; r<3 ; r++)   
      normalGlob[r]/= TotalCellArea;
    
    double temp=std::sqrt(normalGlob[0]*normalGlob[0]+normalGlob[1]*normalGlob[1]+normalGlob[2]*normalGlob[2]);
    
    if(temp>0)
      for (int r=0 ; r<3 ; r++)   
        normalGlob[r]/=temp;
   
    double areaRatio=TotalCellArea/TotalCellRestingArea; 
    
    strainZ=strainZ/TotalCellRestingArea; 
    // if( parameter(4)==5){
    //   cellData[cellIndex][areaRatioIndex]=TotalCellRestingArea;
    // }  
      
    


    // eigenvalue/eigenvectors of averaged STRAIN and STRESS tensors in global coordinate system. (Jacobi method)
    if(neighborweight>0){// in this case stress calculations are done in a seperate loop
      cellData[cellIndex][stressTensorIndex  ]=StressCellGlobal[0][0];
      cellData[cellIndex][stressTensorIndex+1]=StressCellGlobal[1][1];
      cellData[cellIndex][stressTensorIndex+2]=StressCellGlobal[2][2];
      cellData[cellIndex][stressTensorIndex+3]=StressCellGlobal[0][1];
      cellData[cellIndex][stressTensorIndex+4]=StressCellGlobal[0][2];
      cellData[cellIndex][stressTensorIndex+5]=StressCellGlobal[1][2];
      
      cellData[cellIndex][normalVectorIndex  ]=normalGlob[0];
      cellData[cellIndex][normalVectorIndex+1]=normalGlob[1];
      cellData[cellIndex][normalVectorIndex+2]=normalGlob[2];
    }

    cellData[cellIndex][normalVectorIndex  ]=normalGlob[0];
    cellData[cellIndex][normalVectorIndex+1]=normalGlob[1];
    cellData[cellIndex][normalVectorIndex+2]=normalGlob[2];

    if(neighborweight==0){
      
      // stress component along MT direction 
      // cellData[cellIndex][MTstressIndex ] =
      //   cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex  ] *StressCellGlobal[0][0]  +
      //   cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex+1] *StressCellGlobal[0][1]  +
      //   cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex+2] *StressCellGlobal[0][2]  +
      //   cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex  ] *StressCellGlobal[1][0]  +
      //   cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex+1] *StressCellGlobal[1][1]  +
      //   cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex+2] *StressCellGlobal[1][2]  +
      //   cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex  ] *StressCellGlobal[2][0]  +
      //   cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex+1] *StressCellGlobal[2][1]  +
      //   cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex+2] *StressCellGlobal[2][2]       ;
      
      
      
      // eigenvalue/eigenvectors of  stress tensor in global coordinate system. (Jacobi method)
      
      double pi=3.14159265;
      
      int I,J;    
      double pivot=1;
      double tanRotAngle,Si,Co;
      double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      pivot=1;
      
      while (pivot>stressEpcilon) {
        pivot=std::fabs(StressCellGlobal[1][0]);
        I=1;
        J=0;
        if (std::fabs(StressCellGlobal[2][0])>pivot) {
          pivot=std::fabs(StressCellGlobal[2][0]);
          I=2;
          J=0;
        }
        if (std::fabs(StressCellGlobal[2][1])>pivot) {
          pivot=std::fabs(StressCellGlobal[2][1]);
          I=2;
          J=1;
        }
        if (std::fabs(StressCellGlobal[I][I]-StressCellGlobal[J][J])<stressEpcilon) {
          
          //RotAngle=pi/4;
          Si=0.70710678118;
          Co=0.70710678118;
        }            
        else {
          tanRotAngle=(2*StressCellGlobal[I][J])/(StressCellGlobal[J][J]-StressCellGlobal[I][I]);
          tanRotAngle=1/std::sqrt(1+tanRotAngle*tanRotAngle);
          if(tanRotAngle<1)
            Si=0.70710678118*std::sqrt(1-tanRotAngle);
          else
            Si=0.70710678118*std::sqrt(tanRotAngle-1);
          Co=0.70710678118*std::sqrt(1+tanRotAngle);
        }
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;
        
        double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStress[r][s]=tempStress[r][s]+StressCellGlobal[r][w]*tempRot[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StressCellGlobal[r][s]=0;
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StressCellGlobal[r][s]=StressCellGlobal[r][s]+tempRot[w][r]*tempStress[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStress[r][s]=eigenVectorStress[r][s];
        
        for (size_t ii=0; ii<3; ++ii) {
          for (size_t jj=0; jj<3; ++jj) {
            eigenVectorStress[ii][jj] = 0.0;
          }
        }
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
      }
      
      // normalizing eigenvectors (remove if not necessary)  
      double temp=std::sqrt(eigenVectorStress[0][0]*eigenVectorStress[0][0] +
                            eigenVectorStress[1][0]*eigenVectorStress[1][0] +
                            eigenVectorStress[2][0]*eigenVectorStress[2][0] );
      if(temp>0){
        eigenVectorStress[0][0]/=temp;
        eigenVectorStress[1][0]/=temp;
        eigenVectorStress[2][0]/=temp;
      }
      temp=std::sqrt(eigenVectorStress[0][1]*eigenVectorStress[0][1] +
                     eigenVectorStress[1][1]*eigenVectorStress[1][1] +
                     eigenVectorStress[2][1]*eigenVectorStress[2][1] );
      if(temp>0){
        eigenVectorStress[0][1]/=temp;
        eigenVectorStress[1][1]/=temp;
        eigenVectorStress[2][1]/=temp;
      }
      temp=std::sqrt(eigenVectorStress[0][2]*eigenVectorStress[0][2] +
                     eigenVectorStress[1][2]*eigenVectorStress[1][2] +
                     eigenVectorStress[2][2]*eigenVectorStress[2][2] );
      if(temp>0){
        eigenVectorStress[0][2]/=temp;
        eigenVectorStress[1][2]/=temp;
        eigenVectorStress[2][2]/=temp;
      } 
      
      // maximal stress direction
      double maximalStressValue=StressCellGlobal[0][0];
      int Istress=0;
      if ( std::fabs(StressCellGlobal[1][1]) > std::fabs(maximalStressValue) ) 
        {
          maximalStressValue=StressCellGlobal[1][1];
          Istress=1;
        }
      if ( std::fabs(StressCellGlobal[2][2]) > std::fabs(maximalStressValue) ) 
        {
          maximalStressValue=StressCellGlobal[2][2];
          Istress=2;
        }
      
      
      // 2nd maximalstress direction/value
      double maximalStressValue2, maximalStressValue3;
      int Istress2,Istress3;
      if (Istress==0) {
        Istress2=1;
        Istress3=2;
      }
      if (Istress==1) {
        Istress2=2;
        Istress3=0;
      }
      if (Istress==2) {
        Istress2=0;
        Istress3=1;
      }
      if(std::fabs(StressCellGlobal[Istress3][Istress3]) > std::fabs(StressCellGlobal[Istress2][Istress2]) ) {
        temp=Istress2;
        Istress2=Istress3;
        Istress3=temp;
      }
      maximalStressValue2=StressCellGlobal[Istress2][Istress2];
      maximalStressValue3=StressCellGlobal[Istress3][Istress3];
      
      
      if(std::fabs(normalGlob[0]*eigenVectorStress[0][Istress]+
                   normalGlob[1]*eigenVectorStress[1][Istress]+
                   normalGlob[2]*eigenVectorStress[2][Istress]) > .7) {
        // std::cerr << "max stress normal to the cell plane "<<cellIndex<<std::endl;
        Istress=Istress2; 
        Istress2=Istress3;
        maximalStressValue=maximalStressValue2;
        maximalStressValue2=maximalStressValue3; 
      }
      
      // von mises stress for principal plane stress
      double misesStress=std::sqrt( maximalStressValue    * maximalStressValue  
                                    + maximalStressValue2 * maximalStressValue2
                                    - maximalStressValue  * maximalStressValue2);

      cellData[cellIndex][MisesStressIndex ]=misesStress;
      // storing a measure for stress anisotropy in cell vector
      if (std::fabs(maximalStressValue)<  0.000001) cellData[cellIndex][stressAnIndex]=0;
      
      if (std::fabs(maximalStressValue)>= 0.000001) {
        if(parameter(6)==0)
          cellData[cellIndex][stressAnIndex]=1-std::fabs(maximalStressValue2/maximalStressValue);
        else
          cellData[cellIndex][stressAnIndex]=(1-std::fabs(maximalStressValue2/maximalStressValue))*(maximalStressValue/parameter(6));
      }    
      
      
      
      
      
      if (numVariableIndexLevel()==4 &&(numVariableIndex(3)==1 || numVariableIndex(3)==2 ) ) { // storing maximal stress
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
            cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
          }
        if (dimension==3)
          { 
            cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
            cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
            cellData[cellIndex][variableIndex(3,0)+2]=eigenVectorStress[2][Istress];
            
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
          }
      }
      
      
      if (numVariableIndexLevel()==4 && numVariableIndex(3)==2 ) { // storing 2nd maximal stress
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(3,1)]  =eigenVectorStress[0][Istress2];
            cellData[cellIndex][variableIndex(3,1)+1]=eigenVectorStress[1][Istress2];
            cellData[cellIndex][variableIndex(3,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
          }
        if (dimension==3)
          {
            cellData[cellIndex][variableIndex(3,1)]  =eigenVectorStress[0][Istress2];
            cellData[cellIndex][variableIndex(3,1)+1]=eigenVectorStress[1][Istress2];
            cellData[cellIndex][variableIndex(3,1)+2]=eigenVectorStress[2][Istress2];
            cellData[cellIndex][variableIndex(3,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
          }
      }
      
      // if( parameter(4)==5){// misses stress is stored instead of iso energy
      //   cellData[cellIndex][isoEnergyIndex]=std::sqrt(maximalStressValue*maximalStressValue+
      //                                                  maximalStressValue2*maximalStressValue2);      
      // }  
      
      
    
    }// end of neighborweight==0
  
   

    // STRAIN:
    

    double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    double pivot=1;
    double pi=3.1415927;
    int I,J;
    double tanRotAngle,Si,Co;
    


    pivot=std::fabs(StrainCellGlobal[1][0]);
    I=1;
    J=0;
    if (std::fabs(StrainCellGlobal[2][0])>pivot) {
      pivot=std::fabs(StrainCellGlobal[2][0]);
      I=2;
      J=0;
    }
    if (std::fabs(StrainCellGlobal[2][1])>pivot) {
      pivot=std::fabs(StrainCellGlobal[2][1]);
      I=2;
      J=1;
    }

    
    

    while (pivot> strainEpcilon) {  
      if (std::fabs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])< strainEpcilon ) {
        //RotAngle=pi/4;
        Si=0.70710678118;
        Co=0.70710678118;
      }            
      else {
        tanRotAngle=(2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]);
        tanRotAngle=1/std::sqrt(1+tanRotAngle*tanRotAngle);
        if(tanRotAngle<1)
          Si=0.70710678118*std::sqrt(1-tanRotAngle);
        else
          Si=0.70710678118*std::sqrt(tanRotAngle-1);
        Co=0.70710678118*std::sqrt(1+tanRotAngle);
      }
      double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      tempRot[I][I]=Co;
      tempRot[J][J]=Co;
      tempRot[I][J]=Si;
      tempRot[J][I]=-Si;
      double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
      
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          StrainCellGlobal[r][s]=0;
      
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
      
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          tempStrain[r][s]=eigenVectorStrain[r][s];
      
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          eigenVectorStrain[r][s]=0;
      
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
      
      
      pivot=std::fabs(StrainCellGlobal[1][0]);
      I=1;
      J=0;
      if (std::fabs(StrainCellGlobal[2][0])>pivot) {
        pivot=std::fabs(StrainCellGlobal[2][0]);
        I=2;
        J=0;
      }
      if (std::fabs(StrainCellGlobal[2][1])>pivot) {
        pivot=std::fabs(StrainCellGlobal[2][1]);
        I=2;
        J=1;
      }
    }
    

    

    
    // normalizing eigenvectors (remove if not necessary)  
    temp=std::sqrt(eigenVectorStrain[0][0]*eigenVectorStrain[0][0] +
                   eigenVectorStrain[1][0]*eigenVectorStrain[1][0] +
                   eigenVectorStrain[2][0]*eigenVectorStrain[2][0] );
    if(temp>0){
      eigenVectorStrain[0][0]/=temp;
      eigenVectorStrain[1][0]/=temp;
      eigenVectorStrain[2][0]/=temp;
    }
    temp=std::sqrt(eigenVectorStrain[0][1]*eigenVectorStrain[0][1] +
                   eigenVectorStrain[1][1]*eigenVectorStrain[1][1] +
                   eigenVectorStrain[2][1]*eigenVectorStrain[2][1] );
    if(temp>0){
      eigenVectorStrain[0][1]/=temp;
      eigenVectorStrain[1][1]/=temp;
      eigenVectorStrain[2][1]/=temp;
    }
    temp=std::sqrt(eigenVectorStrain[0][2]*eigenVectorStrain[0][2] +
                   eigenVectorStrain[1][2]*eigenVectorStrain[1][2] +
                   eigenVectorStrain[2][2]*eigenVectorStrain[2][2] );
    if(temp>0){
      eigenVectorStrain[0][2]/=temp;
      eigenVectorStrain[1][2]/=temp;
      eigenVectorStrain[2][2]/=temp;
    }
 
      // maximal strain direction
      double maximalStrainValue=StrainCellGlobal[0][0];
      int Istrain=0;
      if (std::fabs(StrainCellGlobal[1][1])>std::fabs(maximalStrainValue)) 
        {
          maximalStrainValue=StrainCellGlobal[1][1];
          Istrain=1;
        }
      if (std::fabs(StrainCellGlobal[2][2])>std::fabs(maximalStrainValue)) 
        {
          maximalStrainValue=StrainCellGlobal[2][2];
          Istrain=2;
        }
   

      // 2nd maximalstrain direction/value
      double maximalStrainValue2,maximalStrainValue3;
      int Istrain2,Istrain3;
      if (Istrain==0) {
        Istrain2=1;
        Istrain3=2;
      }
      if (Istrain==1) {
        Istrain2=0;
        Istrain3=2;
      }
      if (Istrain==2) {
        Istrain2=0;
        Istrain3=1;
      }
      if(std::fabs(StrainCellGlobal[Istrain3][Istrain3]) > std::fabs(StrainCellGlobal[Istrain2][Istrain2] )) {
        temp=Istrain2;
        Istrain2=Istrain3;
        Istrain3=temp;
      }
      maximalStrainValue2=StrainCellGlobal[Istrain2][Istrain2];
      maximalStrainValue3=StrainCellGlobal[Istrain3][Istrain3];
      
      
      if(std::fabs(normalGlob[0]*eigenVectorStrain[0][Istrain]+   //<<<<<<<<<<<<<<<<<<
                  normalGlob[1]*eigenVectorStrain[1][Istrain]+
                  normalGlob[2]*eigenVectorStrain[2][Istrain]) > 1) {

	std::cerr <<cellIndex;
        std::cerr << "max strain vector is out of cell plane in the cell " <<cellIndex <<std::endl;
     
      }
      
      
      
      if(maximalStrainValue != maximalStrainValue){
        
        std::cerr<<"mechanicalTRBS maximal strain here "<<maximalStrainValue<<std::endl;
        std::cerr<<"in the cell "<<cellIndex<<std::endl;
        exit(-1);
      }


     


      // ad-hoc to check resting length estimation for strain maintenance
      // if (timeC>100 && !lengthout ){

      //   for (size_t wallindex=0; wallindex<numWalls; ++wallindex) {
          
      //     size_t kPlusOneMod = (wallindex+1)%numWalls;
      //     size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
      //     size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //     //size_t w1 = internal wallindex
      //     size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //     //size_t w3 = internal wallindex+1
          
      //     // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      //     DataMatrix position(3,vertexData[v2]);
      //     for (size_t d=0; d<dimension; ++d)
      //       position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //     //position[1] = vertexData[v2]; // given by initiation
      //     position[2] = vertexData[v3];
  
      //     std::vector<double> L(3),Le(3);


          
      //     L[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      //     L[1] = wallData[w2][wallLengthIndex];
      //     L[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];  
          
      //     std::vector<std::vector<double> > wallVec(3);
          
      //      std::vector<double> length(3);
      //      length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
      //                             (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
      //                             (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
           
      //      length[1] = std::sqrt( (position[2][0]-position[1][0])*(position[2][0]-position[1][0]) +
      //                             (position[2][1]-position[1][1])*(position[2][1]-position[1][1]) +
      //                             (position[2][2]-position[1][2])*(position[2][2]-position[1][2]) );
           
      //      length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
      //                             (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
      //                             (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
           
      //     for (size_t d=0; d<3; ++d) 
      //       wallVec[0].push_back(position[0][d]-position[1][d]);          
      //     for (size_t d=0; d<3; ++d) 
      //       wallVec[1].push_back(position[1][d]-position[2][d]);          
      //     for (size_t d=0; d<3; ++d) 
      //       wallVec[2].push_back(position[2][d]-position[0][d]);          
          
      //     for (size_t w=0; w<3; ++w) 
      //       Le[w]=std::sqrt( 
      //                       ( wallVec[w][0]*eigenVectorStrain[0][Istrain]
      //                         +wallVec[w][1]*eigenVectorStrain[1][Istrain]
      //                         +wallVec[w][2]*eigenVectorStrain[2][Istrain])
      //                       *( wallVec[w][0]*eigenVectorStrain[0][Istrain]
      //                          +wallVec[w][1]*eigenVectorStrain[1][Istrain]
      //                          +wallVec[w][2]*eigenVectorStrain[2][Istrain])
      //                       *(1-2*maximalStrainValue)
                            
      //                       +( wallVec[w][0]*eigenVectorStrain[0][Istrain2]
      //                          +wallVec[w][1]*eigenVectorStrain[1][Istrain2]
      //                          +wallVec[w][2]*eigenVectorStrain[2][Istrain2])
      //                       *( wallVec[w][0]*eigenVectorStrain[0][Istrain2]
      //                          +wallVec[w][1]*eigenVectorStrain[1][Istrain2]
      //                          +wallVec[w][2]*eigenVectorStrain[2][Istrain2])
      //                       *(1-2*maximalStrainValue2)
                            
      //                       +( wallVec[w][0]*eigenVectorStrain[0][Istrain3]
      //                          +wallVec[w][1]*eigenVectorStrain[1][Istrain3]
      //                          +wallVec[w][2]*eigenVectorStrain[2][Istrain3])
      //                       *( wallVec[w][0]*eigenVectorStrain[0][Istrain3]
      //                          +wallVec[w][1]*eigenVectorStrain[1][Istrain3]
      //                          +wallVec[w][2]*eigenVectorStrain[2][Istrain3])
      //                       *(1-2*maximalStrainValue3)
      //                        );
          
      //     std::cerr<<length[0]<<"  "<<L[0]<<"  "<<Le[0]<<"  "
      //              <<length[0]-L[0]<<"  "<<length[0]-Le[0]<<"  "
      //              <<(length[0]-L[0])/length[0]<<"  "<<(length[0]-Le[0])/length[0]<<"  "
      //              <<(length[0]*length[0]-L[0]*L[0])/(length[0]*length[0])<<"  "
      //              <<(length[0]*length[0]-Le[0]*Le[0])/(length[0]*length[0])<<"  "
      //              <<std::endl;
      //     std::cerr<<length[1]<<"  "<<L[1]<<"  "<<Le[1]<<"  "
      //              <<length[1]-L[1]<<"  "<<length[1]-Le[1]<<"  "
      //              <<(length[1]-L[1])/length[1]<<"  "<<(length[1]-Le[1])/length[1]<<"  "
      //              <<(length[1]*length[1]-L[1]*L[1])/(length[1]*length[1])<<"  "
      //              <<(length[1]*length[1]-Le[1]*Le[1])/(length[1]*length[1])<<"  "
      //              <<std::endl;
      //     std::cerr<<length[2]<<"  "<<L[2]<<"  "<<Le[2]<<"  "
      //              <<length[2]-L[2]<<"  "<<length[2]-Le[2]<<"  "
      //              <<(length[2]-L[2])/length[2]<<"  "<<(length[2]-Le[2])/length[2]<<"  "
      //              <<(length[2]*length[2]-L[2]*L[2])/(length[2]*length[2])<<"  "
      //              <<(length[2]*length[2]-Le[2]*Le[2])/(length[2]*length[2])<<"  "
      //              <<std::endl;

      //   }
        
      // }
      // END of ad-hoc to check resting length estimation for strain maintenance


     
      double growthStrain=maximalStrainValue2;
      if ( cellData[cellIndex][MTindex  ]*eigenVectorStrain[0][Istrain] +
	   cellData[cellIndex][MTindex+1]*eigenVectorStrain[1][Istrain] +
	   cellData[cellIndex][MTindex+2]*eigenVectorStrain[2][Istrain] < 0.01  ){
      growthStrain=maximalStrainValue;
      }

      // normal to the cell plane in global coordinate is averaged Zcurrent[], vector product gives the perpendicular strain direction
      double PerpStrain[3];
      PerpStrain[0]=normalGlob[1]*eigenVectorStrain[2][Istrain]-normalGlob[2]*eigenVectorStrain[1][Istrain];
      PerpStrain[1]=normalGlob[2]*eigenVectorStrain[0][Istrain]-normalGlob[0]*eigenVectorStrain[2][Istrain];
      PerpStrain[2]=normalGlob[0]*eigenVectorStrain[1][Istrain]-normalGlob[1]*eigenVectorStrain[0][Istrain];
      temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     
      
       if(std::fabs(temp)<0.0001){ // if maximal strain is normal to the cell plane storing strain direction instead as it should not be used
         PerpStrain[0]=eigenVectorStrain[0][Istrain];
         PerpStrain[1]=eigenVectorStrain[1][Istrain];
         PerpStrain[2]=eigenVectorStrain[2][Istrain];
       }

       
    // storing a measure for strain anisotropy in cell vector
    if (std::fabs(maximalStrainValue) <  0.000001) cellData[cellIndex][strainAnIndex]=0;
    if (std::fabs(maximalStrainValue) >= 0.000001) cellData[cellIndex][strainAnIndex]=1-std::fabs(maximalStrainValue2/maximalStrainValue);// relative
    //if (std::fabs(maximalStrainValue) >= 0.000001) cellData[cellIndex][strainAnIndex]=(maximalStrainValue-maximalStrainValue2)/0.08;// absolute

    // storing strain/stress direction/value in cellData
    if (numVariableIndexLevel()==4 && (numVariableIndex(2)==1 || numVariableIndex(2)==2 || numVariableIndex(2)==3) ) {// storing maximal strain
      if (dimension==2)
        {
          cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
          cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
          cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
        }
      if (dimension==3)
        {
          cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
          cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
          cellData[cellIndex][variableIndex(2,0)+2]=eigenVectorStrain[2][Istrain];
          cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
        }
    }
    //std::cerr<< maximalStrainValue<< std::endl;
    if (numVariableIndexLevel()==4 && numVariableIndex(2)==3) {//storing 2nd maximal strain
      if (dimension==2)
        {
          cellData[cellIndex][variableIndex(2,2)]  =eigenVectorStrain[0][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+1]=eigenVectorStrain[1][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
        }
      if (dimension==3)
        {
          cellData[cellIndex][variableIndex(2,2)]  =eigenVectorStrain[0][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+1]=eigenVectorStrain[1][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+2]=eigenVectorStrain[2][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
        }
    }
    
  
    
    if (numVariableIndexLevel()==4 && ( numVariableIndex(2)==2 || numVariableIndex(2)==3) ) {//storing perpendicular to maximal strain
      if (dimension==2)
        {
          cellData[cellIndex][variableIndex(2,1)]  =PerpStrain[0];
          cellData[cellIndex][variableIndex(2,1)+1]=PerpStrain[1];
          cellData[cellIndex][variableIndex(2,1)+3]=maximalStrainValue2;  //maximal Strain Value is stored after its eigenvector
          
        }
      if (dimension==3)
        { 
          cellData[cellIndex][variableIndex(2,1)]  =PerpStrain[0];
          cellData[cellIndex][variableIndex(2,1)+1]=PerpStrain[1];
          cellData[cellIndex][variableIndex(2,1)+2]=PerpStrain[2];
          cellData[cellIndex][variableIndex(2,1)+3]=maximalStrainValue2;//growthStrain;  //growth Strain Value is stored after its eigenvector
          
        }
    }
    
    
    
    // normalizing MT vector
    temp=std::sqrt(cellData[cellIndex][MTindex  ]*cellData[cellIndex][MTindex  ] +
                   cellData[cellIndex][MTindex+1]*cellData[cellIndex][MTindex+1] +
                   cellData[cellIndex][MTindex+2]*cellData[cellIndex][MTindex+2]   );
    if(temp>0){      
      cellData[cellIndex][MTindex  ]=cellData[cellIndex][MTindex  ]/temp;
      cellData[cellIndex][MTindex+1]=cellData[cellIndex][MTindex+1]/temp;
      cellData[cellIndex][MTindex+2]=cellData[cellIndex][MTindex+2]/temp;
    }


 
    //<<<<<<<<<<<<<<<<<<<<<< angles between  vectors and circumferential direction <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   
    //cellData[cellIndex][15]= cellData[cellIndex][comIndex+2]; // z coordinate of central vertex of the cell // 15 --> Z coordinate
        
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
    
    //cellData[cellIndex][areaRatioIndex  ]=areaRatio;
    cellData[cellIndex][areaRatioIndex  ]=TotalCellRestingArea;
    //cellData[cellIndex][23]=areaRatio;TotalCellRestingArea;
    

    // if(TotalCellRestingArea<180  )
    //   cellData[cellIndex][37]=-4;
    // if(TotalCellRestingArea<400  && TotalCellRestingArea>390  )
    //   cellData[cellIndex][37]=-3;
    // if(TotalCellRestingArea<620  && TotalCellRestingArea>610  )
    //   cellData[cellIndex][37]=-2;
    // if(TotalCellRestingArea<780  && TotalCellRestingArea>760  )
    //   cellData[cellIndex][37]=-1;
        




    cellData[cellIndex][isoEnergyIndex  ]=EnergyIso/TotalCellRestingArea;
    cellData[cellIndex][anisoEnergyIndex]=EnergyAniso/TotalCellRestingArea;
    
    
    
    
    //cellData[cellIndex][isoEnergyIndex  ]= TotalCellRestingArea;
    //cellData[cellIndex][anisoEnergyIndex]=EnergyAniso/TotalCellRestingArea;
    
    
    
  }
  if(timeC>100)
    lengthout=true;
  
  
  

  // for overal energy calculation  
  // double totalEnergyIso=0;
  // double totalEnergyAniso=0;
  //  for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
  //   totalEnergyIso   +=cellData[cellIndex][isoEnergyIndex];
  //   totalEnergyAniso +=cellData[cellIndex][anisoEnergyIndex];
  //  }




  if(neighborweight>0)
    for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
      
      const size_t numWalls = T.cell(cellIndex).numWall();
      
      
      
      Cell *  cell1=&(T.cell(cellIndex));
      
      //std::vector<int> neighbor(numWalls);
      std::vector<size_t> neighbor(numWalls);
      for   ( size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex ){
        neighbor[wallIndex]=(cell1->cellNeighbor(wallIndex))->index();
      }
      
      // std::cerr<<" cell   "<<cellIndex << "cell neighbors    ";
      // for   ( size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex ){
      //   std::cerr<< "  "<<neighbor[wallIndex]<<" " ;
      // }
      
      
      
      double normalGlob[3]={0,0,0};
      
      normalGlob[0]= cellData[cellIndex][normalVectorIndex  ];
      normalGlob[1]= cellData[cellIndex][normalVectorIndex+1];
      normalGlob[2]= cellData[cellIndex][normalVectorIndex+2]; 
      
      double StressTensor[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      
      int counter=0;
      for (size_t nn=0 ; nn<numWalls ; nn++){
        if (neighbor[nn]<numCells && neighbor[nn]>-1){
          
          StressTensor[0][0]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex  ];
          StressTensor[1][1]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+1];
          StressTensor[2][2]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+2];
          StressTensor[0][1]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+3];
          StressTensor[2][0]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+4];
          StressTensor[1][2]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+5];
          counter+=1;
        }
      }
      if(counter !=0){
        StressTensor[0][0]/=counter;
        StressTensor[1][1]/=counter;
        StressTensor[2][2]/=counter;
        StressTensor[0][1]/=counter;
        StressTensor[2][0]/=counter;
        StressTensor[1][2]/=counter;
      }
      
      StressTensor[0][0]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex  ];
      StressTensor[1][1]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+1];
      StressTensor[2][2]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+2];
      StressTensor[0][1]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+3];
      StressTensor[2][0]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+4];
      StressTensor[1][2]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+5];
      
      
      
      StressTensor[0][2]=StressTensor[2][0];
      StressTensor[1][0]=StressTensor[0][1];
      StressTensor[2][1]=StressTensor[1][2];
      
      
      // stress component along MT direction 
      // cellData[cellIndex][MTstressIndex ] =
      //   cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex  ] *StressTensor[0][0]  +
      //   cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex+1] *StressTensor[0][1]  +
      //   cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex+2] *StressTensor[0][2]  +
      //   cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex  ] *StressTensor[1][0]  +
      //   cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex+1] *StressTensor[1][1]  +
      //   cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex+2] *StressTensor[1][2]  +
      //   cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex  ] *StressTensor[2][0]  +
      //   cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex+1] *StressTensor[2][1]  +
      //   cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex+2] *StressTensor[2][2]       ;
      
      
      
      // eigenvalue/eigenvectors of  stress tensor in global coordinate system. (Jacobi method)
      
      double pi=3.14159265;
      
      int I,J;    
      double pivot=1;
      double tanRotAngle,Si,Co;
      double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      pivot=1;
      
      while (pivot>stressEpcilon) {
        pivot=std::fabs(StressTensor[1][0]);
        I=1;
        J=0;
        if (std::fabs(StressTensor[2][0])>pivot) {
          pivot=std::fabs(StressTensor[2][0]);
          I=2;
          J=0;
        }
        if (std::fabs(StressTensor[2][1])>pivot) {
          pivot=std::fabs(StressTensor[2][1]);
          I=2;
          J=1;
        }
        if (std::fabs(StressTensor[I][I]-StressTensor[J][J])<stressEpcilon) {
          
          //RotAngle=pi/4;
          Si=0.70710678118;
          Co=0.70710678118;
        }            
        else {
          tanRotAngle=(2*StressTensor[I][J])/(StressTensor[J][J]-StressTensor[I][I]);
          tanRotAngle=1/std::sqrt(1+tanRotAngle*tanRotAngle);
          if(tanRotAngle<1)
            Si=0.70710678118*std::sqrt(1-tanRotAngle);
          else
            Si=0.70710678118*std::sqrt(tanRotAngle-1);
          Co=0.70710678118*std::sqrt(1+tanRotAngle);
        }
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;
        
        double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStress[r][s]=tempStress[r][s]+StressTensor[r][w]*tempRot[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StressTensor[r][s]=0;
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StressTensor[r][s]=StressTensor[r][s]+tempRot[w][r]*tempStress[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStress[r][s]=eigenVectorStress[r][s];
        
        for (size_t ii=0; ii<3; ++ii) {
          for (size_t jj=0; jj<3; ++jj) {
            eigenVectorStress[ii][jj] = 0.0;
          }
        }
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
      }
      
      // normalizing eigenvectors (remove if not necessary)  
      double temp=std::sqrt(eigenVectorStress[0][0]*eigenVectorStress[0][0] +
                            eigenVectorStress[1][0]*eigenVectorStress[1][0] +
                            eigenVectorStress[2][0]*eigenVectorStress[2][0] );
      if(temp>0){
        eigenVectorStress[0][0]/=temp;
        eigenVectorStress[1][0]/=temp;
        eigenVectorStress[2][0]/=temp;
      }
      temp=std::sqrt(eigenVectorStress[0][1]*eigenVectorStress[0][1] +
                     eigenVectorStress[1][1]*eigenVectorStress[1][1] +
                     eigenVectorStress[2][1]*eigenVectorStress[2][1] );
      if(temp>0){
        eigenVectorStress[0][1]/=temp;
        eigenVectorStress[1][1]/=temp;
        eigenVectorStress[2][1]/=temp;
      }
      temp=std::sqrt(eigenVectorStress[0][2]*eigenVectorStress[0][2] +
                     eigenVectorStress[1][2]*eigenVectorStress[1][2] +
                     eigenVectorStress[2][2]*eigenVectorStress[2][2] );
      if(temp>0){
        eigenVectorStress[0][2]/=temp;
        eigenVectorStress[1][2]/=temp;
        eigenVectorStress[2][2]/=temp;
      } 
      
      // maximal stress direction
      double maximalStressValue=StressTensor[0][0];
      int Istress=0;
      if ( StressTensor[1][1] > maximalStressValue ) 
        {
          maximalStressValue=StressTensor[1][1];
          Istress=1;
        }
      if ( StressTensor[2][2] > maximalStressValue ) 
        {
          maximalStressValue=StressTensor[2][2];
          Istress=2;
        }
      
      
      // 2nd maximalstress direction/value
      double maximalStressValue2, maximalStressValue3;
      int Istress2,Istress3;
      if (Istress==0) {
        Istress2=1;
        Istress3=2;
      }
      if (Istress==1) {
        Istress2=2;
        Istress3=0;
      }
      if (Istress==2) {
        Istress2=0;
        Istress3=1;
      }
      if(StressTensor[Istress3][Istress3] > StressTensor[Istress2][Istress2] ) {
        temp=Istress2;
        Istress2=Istress3;
        Istress3=temp;
      }
      maximalStressValue2=StressTensor[Istress2][Istress2];
      maximalStressValue3=StressTensor[Istress3][Istress3];
      
      
      if(std::fabs(normalGlob[0]*eigenVectorStress[0][Istress]+
                   normalGlob[1]*eigenVectorStress[1][Istress]+
                   normalGlob[2]*eigenVectorStress[2][Istress]) > .7) {
        // std::cerr << "max stress normal to the cell plane "<<cellIndex<<std::endl;
        Istress=Istress2; 
        Istress2=Istress3;
        maximalStressValue=maximalStressValue2;
        maximalStressValue2=maximalStressValue3; 
      }
      
      
      
      
      // storing a measure for stress anisotropy in cell vector
      if (std::fabs(maximalStressValue)<  0.000001) cellData[cellIndex][stressAnIndex]=0;
      
      if (std::fabs(maximalStressValue)>= 0.000001) {
        if(parameter(6)==0)
          cellData[cellIndex][stressAnIndex]=1-std::fabs(maximalStressValue2/maximalStressValue);
        else
          cellData[cellIndex][stressAnIndex]=(1-std::fabs(maximalStressValue2/maximalStressValue))*(maximalStressValue/parameter(6));
      }    
      
      
      
      
      
      if (numVariableIndexLevel()==4 &&(numVariableIndex(3)==1 || numVariableIndex(3)==2 ) ) { // storing maximal stress
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
            cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
          }
        if (dimension==3)
          { 
            cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
            cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
            cellData[cellIndex][variableIndex(3,0)+2]=eigenVectorStress[2][Istress];
            
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
          }
      }
      
      
      if (numVariableIndexLevel()==4 && numVariableIndex(3)==2 ) { // storing 2nd maximal stress
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(3,1)]  =eigenVectorStress[0][Istress2];
            cellData[cellIndex][variableIndex(3,1)+1]=eigenVectorStress[1][Istress2];
            cellData[cellIndex][variableIndex(3,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
          }
        if (dimension==3)
          {
            cellData[cellIndex][variableIndex(3,1)]  =eigenVectorStress[0][Istress2];
            cellData[cellIndex][variableIndex(3,1)+1]=eigenVectorStress[1][Istress2];
            cellData[cellIndex][variableIndex(3,1)+2]=eigenVectorStress[2][Istress2];
            cellData[cellIndex][variableIndex(3,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
          }
      }
      
      
     
      
    }
  

  
 
}



void VertexFromTRBScenterTriangulationMT::
printState(Tissue *T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData, 
	   std::ostream &os)
{ // VTK style
  size_t tbsInd=28;  // 1, -1, 0
  size_t L123Ind=29; // 1, 2, 3
  size_t dimension = vertexData[0].size();
  
  size_t numCells = T->numCell();

  static size_t index=0;

  std::stringstream name;
  name << "tmp/VTKL1Top" << index  <<".vtu";
  std::ofstream myfile;
  myfile.open (name.str());
  //  myfile.open ("tmp/VTK_IntWalls.vtu");//("VTK_IntWalls%06d.vtu",index);
  size_t pointCounter=0; 
  size_t cellCounter=0; 




  //////////////////////L1t/////////////////////////////////
  std::vector<std::vector<double> > connectivityL1t;
  std::vector<double> cellIndicesL1t;
  std::vector<std::vector<double> > positionsL1t;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==1 && cellData[cellIndex][L123Ind]==1){
      cellCounter++;
      connectivityL1t.resize(cellCounter);
      cellIndicesL1t.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL1t[cellCounter-1].resize(numCellVertices);
      positionsL1t.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL1t[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL1t[pointCounter-numCellVertices+intc].resize(3);
        positionsL1t[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL1t[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL1t[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL1t[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  size_t tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL1t[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL1t[d]][0]<<" "
           <<cellData[cellIndicesL1t[d]][1]<<" "
           <<cellData[cellIndicesL1t[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL1t[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  size_t cellCompNum=cellData[0].size();
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL1t[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////


//////////////////////L2t/////////////////////////////////

  name.str("");
  name << "tmp/VTKL2Top" << index  <<".vtu";
  myfile.open (name.str());
  
  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL2t;
  std::vector<double> cellIndicesL2t;
  std::vector<std::vector<double> > positionsL2t;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==1 && cellData[cellIndex][L123Ind]==2){
      cellCounter++;
      connectivityL2t.resize(cellCounter);
      cellIndicesL2t.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL2t[cellCounter-1].resize(numCellVertices);
      positionsL2t.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL2t[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL2t[pointCounter-numCellVertices+intc].resize(3);
        positionsL2t[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL2t[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL2t[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL2t[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL2t[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL2t[d]][0]<<" "
           <<cellData[cellIndicesL2t[d]][1]<<" "
           <<cellData[cellIndicesL2t[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL2t[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL2t[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////


//////////////////////L3t////////////////////////////////
  name.str("");
  name << "tmp/VTKL3Top" << index  <<".vtu";
  myfile.open (name.str());

  pointCounter=0; 
  cellCounter=0; 



  std::vector<std::vector<double> > connectivityL3t;
  std::vector<double> cellIndicesL3t;
  std::vector<std::vector<double> > positionsL3t;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==1 && cellData[cellIndex][L123Ind]==3){
      cellCounter++;
      connectivityL3t.resize(cellCounter);
      cellIndicesL3t.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL3t[cellCounter-1].resize(numCellVertices);
      positionsL3t.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL3t[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL3t[pointCounter-numCellVertices+intc].resize(3);
        positionsL3t[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL3t[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL3t[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL3t[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL3t[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL3t[d]][0]<<" "
           <<cellData[cellIndicesL3t[d]][1]<<" "
           <<cellData[cellIndicesL3t[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL3t[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL3t[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////


//////////////////////L1s/////////////////////////////////
  name.str("");
  name << "tmp/VTKL1Side" << index  <<".vtu";
  myfile.open (name.str());


  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL1s;
  std::vector<double> cellIndicesL1s;
  std::vector<std::vector<double> > positionsL1s;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==0 && cellData[cellIndex][L123Ind]==1){
      cellCounter++;
      connectivityL1s.resize(cellCounter);
      cellIndicesL1s.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL1s[cellCounter-1].resize(numCellVertices);
      positionsL1s.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL1s[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL1s[pointCounter-numCellVertices+intc].resize(3);
        positionsL1s[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL1s[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL1s[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL1s[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL1s[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL1s[d]][0]<<" "
           <<cellData[cellIndicesL1s[d]][1]<<" "
           <<cellData[cellIndicesL1s[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL1s[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL1s[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////

//////////////////////L2s/////////////////////////////////
  name.str("");
  name << "tmp/VTKL2Side" << index  <<".vtu";
  myfile.open (name.str());


  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL2s;
  std::vector<double> cellIndicesL2s;
  std::vector<std::vector<double> > positionsL2s;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==0 && cellData[cellIndex][L123Ind]==2){
      cellCounter++;
      connectivityL2s.resize(cellCounter);
      cellIndicesL2s.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL2s[cellCounter-1].resize(numCellVertices);
      positionsL2s.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL2s[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL2s[pointCounter-numCellVertices+intc].resize(3);
        positionsL2s[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL2s[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL2s[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL2s[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL2s[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL2s[d]][0]<<" "
           <<cellData[cellIndicesL2s[d]][1]<<" "
           <<cellData[cellIndicesL2s[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL2s[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL2s[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////

//////////////////////L3s/////////////////////////////////
  name.str("");
  name << "tmp/VTKL3Side" << index  <<".vtu";
  myfile.open (name.str());

  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL3s;
  std::vector<double> cellIndicesL3s;
  std::vector<std::vector<double> > positionsL3s;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==0 && cellData[cellIndex][L123Ind]==3){
      cellCounter++;
      connectivityL3s.resize(cellCounter);
      cellIndicesL3s.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL3s[cellCounter-1].resize(numCellVertices);
      positionsL3s.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL3s[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL3s[pointCounter-numCellVertices+intc].resize(3);
        positionsL3s[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL3s[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL3s[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL3s[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL3s[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL3s[d]][0]<<" "
           <<cellData[cellIndicesL3s[d]][1]<<" "
           <<cellData[cellIndicesL3s[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL3s[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL3s[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////

//////////////////////L1b/////////////////////////////////
  name.str("");
  name << "tmp/VTKL1Bottom" << index  <<".vtu";
  myfile.open (name.str());


  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL1b;
  std::vector<double> cellIndicesL1b;
  std::vector<std::vector<double> > positionsL1b;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==-1 && cellData[cellIndex][L123Ind]==1){
      cellCounter++;
      connectivityL1b.resize(cellCounter);
      cellIndicesL1b.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL1b[cellCounter-1].resize(numCellVertices);
      positionsL1b.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL1b[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL1b[pointCounter-numCellVertices+intc].resize(3);
        positionsL1b[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL1b[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL1b[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL1b[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL1b[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL1b[d]][0]<<" "
           <<cellData[cellIndicesL1b[d]][1]<<" "
           <<cellData[cellIndicesL1b[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL1b[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL1b[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////

//////////////////////L2b/////////////////////////////////
  name.str("");
  name << "tmp/VTKL2Bottom" << index  <<".vtu";
  myfile.open (name.str());

  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL2b;
  std::vector<double> cellIndicesL2b;
  std::vector<std::vector<double> > positionsL2b;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==-1 && cellData[cellIndex][L123Ind]==2){
      cellCounter++;
      connectivityL2b.resize(cellCounter);
      cellIndicesL2b.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL2b[cellCounter-1].resize(numCellVertices);
      positionsL2b.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL2b[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL2b[pointCounter-numCellVertices+intc].resize(3);
        positionsL2b[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL2b[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL2b[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL2b[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL2b[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL2b[d]][0]<<" "
           <<cellData[cellIndicesL2b[d]][1]<<" "
           <<cellData[cellIndicesL2b[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL2b[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL2b[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////

//////////////////////L3b/////////////////////////////////
  name.str("");
  name << "tmp/VTKL3Bottom" << index  <<".vtu";
  myfile.open (name.str());

  pointCounter=0; 
  cellCounter=0; 

  std::vector<std::vector<double> > connectivityL3b;
  std::vector<double> cellIndicesL3b;
  std::vector<std::vector<double> > positionsL3b;
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++) {
    
    if(cellData[cellIndex][tbsInd]==-1 && cellData[cellIndex][L123Ind]==3){
      cellCounter++;
      connectivityL3b.resize(cellCounter);
      cellIndicesL3b.push_back(cellIndex);
      size_t numCellVertices= T->cell(cellIndex).numVertex();
      pointCounter=pointCounter+numCellVertices;
      connectivityL3b[cellCounter-1].resize(numCellVertices);
      positionsL3b.resize(pointCounter);
      for(size_t intc=0; intc<numCellVertices; intc++){
        connectivityL3b[cellCounter-1][intc]=pointCounter-numCellVertices+intc;
        positionsL3b[pointCounter-numCellVertices+intc].resize(3);
        positionsL3b[pointCounter-numCellVertices+intc]
          =vertexData[T->cell(cellIndex).vertex(intc) -> index()];      
      }
      
    }
  }
  // VTK file header
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"<< std::endl
         << "<UnstructuredGrid>"<< std::endl
         << "<Piece  NumberOfPoints=\""<<  pointCounter
         << "\" NumberOfCells=\""<<        cellCounter<<"\">"<< std::endl
         << "<Points>"<< std::endl
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"  << std::endl;
  
 
  
  for (size_t vIndex=0 ; vIndex<pointCounter ; vIndex++){
        
    for( size_t d=0 ; d<dimension ; d++ )
      myfile << positionsL3b[vIndex][d] << " ";
    myfile << std::endl;
       
  }


  //connectivity
  myfile<< "</DataArray>"<<std::endl
        << "</Points>"<<std::endl
        << "<Cells>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl;

  for (size_t cIndex=0 ; cIndex<cellCounter ; cIndex++){
    size_t numCV= connectivityL3b[cIndex].size();
    
    for( size_t d=0 ; d<numCV ; d++ )
      myfile << connectivityL3b[cIndex][d] << " ";
    myfile << std::endl;
    
  }

  //off-sets 
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl;
  tmpC=0; 
  for( size_t d=0 ; d<cellCounter ; d++ ){
    tmpC +=connectivityL3b[d].size();
    myfile <<tmpC<<" ";
  }
  myfile << std::endl;
  
  //types
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<std::endl; 
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<7<<" ";
  myfile << std::endl;
  
 
  //cell vector 
  myfile<< "</DataArray>"<<std::endl
        << "</Cells>"<<std::endl
        << "<CellData Scalars=\"cell variable 3\" Vectors=\"cell vector\">"<<std::endl
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">"<<std::endl;
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL3b[d]][0]<<" "
           <<cellData[cellIndicesL3b[d]][1]<<" "
           <<cellData[cellIndicesL3b[d]][2]<<" ";
  myfile << std::endl;
  
  
  //cell vector length
  myfile<< "</DataArray>"<<std::endl
        << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\"> ";
  for( size_t d=0 ; d<cellCounter ; d++ )
    myfile <<cellData[cellIndicesL3b[d]][3]<<" ";
  myfile << std::endl;
 
  // rest of cell vector 
  
  for(size_t cellCompInd=4; cellCompInd<cellCompNum; cellCompInd++){

    myfile<< "</DataArray>"<<std::endl
          <<"<DataArray type=\"Float64\" Name=\" cell variable "<<cellCompInd
          <<" \" format=\"ascii\">"<<std::endl;
    for( size_t d=0 ; d<cellCounter ; d++ )
      myfile <<cellData[cellIndicesL3b[d]][cellCompInd]<<" ";
    myfile << std::endl;
  }
  
  myfile<< "</DataArray>"<<std::endl
        << "</CellData>"<<std::endl
        << "</Piece>"<<std::endl
        << "</UnstructuredGrid>"<<std::endl
        << "</VTKFile>"<<std::endl;

  myfile.close();
  /////////////////////////////////////////////////


  index++;
}


// initiation is done by an indipendent reaction now <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// void VertexFromTRBScenterTriangulationMT::
// initiate(Tissue &T,
// 	 DataMatrix &cellData,
// 	 DataMatrix &wallData,
// 	 DataMatrix &vertexData,
// 	 DataMatrix &cellDerivs,
// 	 DataMatrix &wallDerivs,
// 	 DataMatrix &vertexDerivs)
// {
//   size_t dimension=3; //Only implemented for 3D models
//   assert (dimension==vertexData[0].size());
//   size_t numVariable = T.cell(0).numVariable();
//   assert (numVariable==cellData[0].size());
//   // Create the new variables
//   if (variableIndex(1,0) != numVariable) {
//     std::cerr << "VertexFromTRBScenterTriangulation::initiate() "
// 	      << "Wrong index given as start index for additional variables."
// 	      << std::endl;
//     exit(-1);
//   }
//   size_t numCell = cellData.size();
//   assert (numCell==T.numCell());
//   std::vector<double> com(dimension);
  
//   for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
//     size_t numInternalWall = T.cell(cellIndex).numVertex();
//     cellData[cellIndex].resize(numVariable+dimension+numInternalWall);
//     cellDerivs[cellIndex].resize(numVariable+dimension+numInternalWall);
//     com = T.cell(cellIndex).positionFromVertex(vertexData);
//     // Set center position to com of the cell
//     for (size_t d=0; d<dimension; ++d)
//       cellData[cellIndex][numVariable+d] = com[d];    
//     // Set internal wall lengths to the distance btw com and the vertex
//     for (size_t wallindex=0; wallindex<numInternalWall; ++wallindex) {
//       Vertex *tmpVertex = T.cell(cellIndex).vertex(wallindex); 
//       size_t vertexIndex = tmpVertex->index();
//       double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
// 				   (com[0]-vertexData[vertexIndex][0])+
// 				   (com[1]-vertexData[vertexIndex][1])*
// 				   (com[1]-vertexData[vertexIndex][1])+
// 				   (com[2]-vertexData[vertexIndex][2])*
// 				   (com[2]-vertexData[vertexIndex][2]) );   
//       cellData[cellIndex][numVariable+dimension+wallindex] = distance;
//     }
//   }
// }

  
} // end namespace Hypocotyl3D
