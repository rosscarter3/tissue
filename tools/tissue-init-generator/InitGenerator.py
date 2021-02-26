# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import random
import math
import numpy as np

def createLineOfCells(numCells,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread):

    # Calculate number of walls and vertices and write to file
    numWalls = 3*numCells+1;
    numVertex = 2*numCells+2;
    chan = open(outputName,"w")
    chan.write(str(numCells)+" "+str(numWalls)+" "+str(numVertex)+"\n")

    #
    # Topology, wall connections
    #
    v1=0;
    v1Add=0;
    v2=1;
    v2Add=0;
    for iw in range(0, numWalls):
        c1=c2=-1;
        # c1
        c1 = int((iw-1)/3);
        if c1==-1:
           c1=0
        # c2 not -1
        if (not(iw % 3) and iw and iw<(numWalls-1)):
    	     c2 = iw/3;
        chan.write(str(iw)+" "+str(c1)+" "+str(c2)+" "+str(v1)+" "+str(v2)+"\n")

        #v1 
        if ( v1 % 2 ): 
            v1+=1
        elif ( v1Add!=0 ): 
            v1+=1
            v1Add=0
        else:
             v1Add+=1

        # v2
        if ( not(v2 % 2) or v2==1):
    	     v2+=1
        elif ( v2Add!=0 ):
    	     v2+=1
    	     v2Add=0
        else: 
            v2Add+=1

         
    #
    # Vertex positions
    #
    wallLength=1.0;
    oX = 0.0;
    oY = 0.0;
    chan.write(str(numVertex)+" "+str(3)+"\n") 

    
    for iv in range(0,numVertex): 
        x = oX + int(iv/2)*wallLength;
        y = oY;
        if( iv % 2 ):
    	     y = oY+wallLength;
        chan.write(str(x)+" "+str(y)+" "+str(0)+"\n")

   
    
    #
    # Wall variables
    #
    chan.write(str(numWalls)+" "+str(1)+" "+str(2*numWallVariables)+"\n")
    for iw in range(0, numWalls):
        chan.write(str(wallLength)+" ")
        for j in range(0, numWallVariables):
            wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
        chan.write("\n")
    
    #
    # Cell variables
    #
    chan.write(str(numCells)+" "+str(numCellVariables+4)+"\n")
    for i in range(0,numCells): 
        chan.write( "1 0 1 0.9 ")
        for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
        chan.write("\n")
    

 

def createCircleOfCells(numCells,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,wallsSpread):

    # Calculate number of walls and vertices and write to file
    numWalls = 3*numCells;
    numVertex = 2*numCells;
    chan = open(outputName,"w")
    chan.write(str(numCells)+" "+str(numWalls)+" "+str(numVertex)+"\n")

    #   
    # Topology, wall connections
    #
    v1=0;
    v1Add=0;
    v2=1;
    v2Add=0;
    for iw in range(0, numWalls):
        c1=c2=-1;
        # c1
        c1 = int((iw-1)/3);
        if c1==-1:
           c1=0

        # c2 not -1
        if (not (iw % 3) and iw and iw<(numWalls-1)):
	
    	     c2 = iw/3;
        elif ( not iw): 
	    c2 = numCells-1;
        chan.write(" "+str(iw)+" "+str(c1)+" "+str(c2)+" "+str(v1)+" "+str(v2)+"\n")
  
        #v1 
        if ( v1 % 2 ): 
            v1+=1
            v1 = v1%numVertex;
        elif ( v1Add!=0 ): 
            v1+=1
            v1 = v1%numVertex;
            v1Add=0
        else:
             v1Add+=1

        # v2
        if ( not(v2 % 2) or v2==1):
            v2+=1
            v2 = v2%numVertex;
        elif ( v2Add!=0 ):
            v2+=1
            v2 = v2%numVertex;
            v2Add=0
        else: 
            v2Add+=1
       
    #
    # Vertex positions
    #
    wallLength=1.0;
    oX = 0.0;
    oY = 0.0;
    PI = 3.14159265;
    r = numCells*wallLength/(2.0*PI);
    deltaPhi = 2.0*PI/numCells;
    phi=0.0;

    chan.write("\n"+str(numVertex)+" "+str(3)+"\n")
    
    for iv in range(0,numVertex): 
        x = oX + int(iv/2)*wallLength;
        y = oY;
        if( iv % 2 ):
            x = oX + (r+wallLength)*math.cos(phi);
            y = oY + (r+wallLength)*math.sin(phi);
            phi += deltaPhi;
        else:
	         x = oX + r*math.cos(phi);
	         y = oY + r*math.sin(phi);
        chan.write(" "+str(x)+" "+str(y)+" 0\n")
        
    
    #
    # Wall variables
    #
    chan.write(str(numWalls)+" "+str(1)+" "+str(2*numWallVariables)+"\n")
    for iw in range(0, numWalls):
        chan.write(str(wallLength)+" ")
        for j in range(0, numWallVariables):
            wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(wallsData1)+" "+str(wallsData2)+" ") 
        chan.write("\n")

    #
    # Cell variables
    #
    chan.write("\n"+str(numCells)+" "+str(numCellVariables+4)+"\n")
    for i in range(0,numCells): 
        chan.write( "1 0 1 0.9 ")
        for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
        chan.write("\n")
    
    
    
def createSquareGrid(nGrid,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,wallsSpread):  
    n=nGrid
    global nCells, nWALLS, nVERTEX
    c=0
    nCells= 0
    nWALLS=0
    nVERTEX=0
    chan = open(outputName,"w")
    nCells= n**2
    nWALLS=4+6*(n-1)+2*(n-1)**2
    nVERTEX=(n+1)**2
    ncount=n-1 
    vertexone=0
    vertextwo=1
    ncount3=3*n
    chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   #prints first line of init file
    while c<=ncount: #first row of vertical walls
       chan.write(" "+str(c)+" "+str(c)+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
       c=c+1
       vertexone=vertexone+1
       vertextwo=vertextwo+1
    
    
    cellone=-1
    celltwo=0
    vertexone=n+1
    vertextwo=celltwo
    wallcounter=1
    celloneA=0
    celltwoA=celltwo
    vertexoneA=n+1
    vertextwoA=n+2
    
    
    while wallcounter<=n-1:
     ncount=ncount+n+1
     chan.write(" "+str(c)+" "+str(celltwo)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
     c=c+1
     vertexone=vertexone+1
     vertextwo=vertextwo+1
     celltwo=celltwo+1
     cellone=cellone+1
     while c<=ncount-1: #horizontal walls
       chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
       c=c+1
       vertexone=vertexone+1
       vertextwo=vertextwo+1
       celltwo=celltwo+1
       cellone=cellone+1
     chan.write(" "+str(c)+" "+str(cellone)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
     vertexone=vertexone+1
     vertextwo=vertextwo+1
     ncount=ncount+n   
     c=c+1
     celltwoA=celltwo
    
     while c<=ncount: #internal vertical walls
      chan.write(" "+str(c)+" "+str(celloneA)+" "+str(celltwoA)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
      c=c+1
      vertextwoA+=1
      vertexoneA+=1
      celltwoA+=1
      celloneA+=1
     vertextwoA+=1
     vertexoneA+=1
     wallcounter+=1
    
    ncount=ncount+n+1
    chan.write(" "+str(c)+" "+str(celltwo)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
    c=c+1
    vertexone=vertexone+1
    vertextwo=vertextwo+1
    celltwo=celltwo+1
    cellone=cellone+1
    while c<=ncount-1: #horizontal walls
       chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
       c=c+1
       vertexone=vertexone+1
       vertextwo=vertextwo+1
       celltwo=celltwo+1
       cellone=cellone+1
    chan.write(" "+str(c)+" "+str(cellone)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
    vertexone=vertexone+1
    vertextwo=vertextwo+1
    ncount=ncount+n   
    c=c+1
    
    while c<=nWALLS-1:
      chan.write(" "+str(c)+" "+str(celloneA)+" "+str(-1)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
      c=c+1
      vertextwoA=vertextwoA+1
      vertexoneA=vertexoneA+1
      celltwoA=celltwoA+1
      celloneA=celloneA+1
    
    chan.write(" \n")
    chan.write(" "+str(nVERTEX)+" 3 \n")
    ycoord=0  
    xcoord=0
    while xcoord<=n:
      while ycoord<=n:
        chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
        ycoord+=1
      ycoord=0
      xcoord+=1
    
    chan.write("\n"+str(nWALLS)+" 1 "+str(numWallVariables*2)+"\n ")
    wallcounter=1
    wallLength=1
    while wallcounter<=nWALLS:
        chan.write(str(wallLength)+" ")
        for j in range(0, numWallVariables):
            wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
            
        chan.write("\n")

        wallcounter+=1
    chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
    wallcounter=1
    
    while wallcounter<=nCells:
    

        chan.write( "1 0 1 0.9 ")
        for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
            
        chan.write("\n")
        wallcounter+=1  
    chan.close()


def createOffsetSquareGrid(nGrid,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,wallsSpread):  
	n=nGrid
	global nCells, nWALLS, nVERTEX
	nCells= 0
	nWALLS=0
	nVERTEX=0


	chan = open(outputName,"w")

	#calculates number of cells, walls and vertices and prints first line of init file
	nCells= n**2
	nWALLS=4*n+(n-1)*(2*n+1)+n*(n+1)
	nVERTEX=2*(2*n+1)+(n-1)*(2*n+2)
	chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   


	#first row of vertical walls
	ncount=2*n-1 
	vertexone=0
	vertextwo=1
	cellTwo=0
	w=0

	while w<=ncount: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1

	#internal row(s) of vertical walls

	vertexone=2*n
	vertextwo=2*n+1
	for i in range(0,n-1,1):
		
		vertexone+=1
	   	vertextwo+=1
		if i==0:
			cellOne=i*n-0.5

		elif (i % 2 == 0): #even 
			cellOne=i*n-0.5
		else:
			cellOne=(i+1)*n-0.5
		if (i % 2 == 0): #even 
			counter=i+1
		else:
			counter=i
		cellTwo=counter*n
		chan.write(" "+str(w)+" "+str(int(cellTwo))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	




		w+=1
		cellTwo+=0.5
		cellOne+=0.5
	  	vertexone+=1
	   	vertextwo+=1
		while w<=(4+2*i)*n-1+i: 



		   chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
		   w+=1
		   cellTwo+=0.5
		   cellOne+=0.5
		   vertexone+=1
	   	   vertextwo+=1

		#internal row(s) of vertical walls
		chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")

		w+=1
		cellTwo+=0.5
		cellOne+=0.5
		vertexone+=1
	   	vertextwo+=1


	#final row of vertical walls
	vertexone+=1
	vertextwo+=1
	cellTwo=n*(n-1)
	while w<=nWALLS-(n*(n+1))-1: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1


	#horizontal walls
	vertexone=0
	vertextwo=2*n+2
	cellTwo=-1

	for i in range(0,n,1):
		
		k=1
		cellTwo=-1
		cellOne=i*n
		while k<=n+1:
			if k==n+1:
				cellOne=-1
			chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
			vertexone+=2
		   	vertextwo+=2
			cellTwo=cellOne
			cellOne+=1

			w+=1
			k+=1
		if (i ==n-2 and i % 2 == 1): #even 
				vertexone+=1
				vertextwo-=0
				
		elif (i ==n-2 and i % 2 == 0): #even 
				vertexone-=1
				vertextwo-=1
				

		elif (i !=n-2 and i % 2 == 0): #even 
			vertexone-=1
			vertextwo-=1
		else:
			vertexone+=1
			vertextwo+=1

	#vertex corrdinates
	chan.write(" \n")
	chan.write(" "+str(nVERTEX)+" 3 \n")
	ycoord=0  
	xcoord=0
	while xcoord<n:
	  while ycoord<=n:
	    chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
	    ycoord+=0.5
	  ycoord=-0.5
	  xcoord+=1
	while xcoord<=n:
	  if xcoord==n and n % 2 == 1:
		ycoord=0
		while ycoord<=n:
		    chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
		    ycoord+=0.5
	 	ycoord=-0.5
	  	xcoord+=1
	  elif xcoord==n and n % 2 == 0:
		ycoord=-0.5
		while ycoord<n:
	    		chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
	    		ycoord+=0.5
	 	ycoord=-0.5
	  	xcoord+=1




	chan.write("\n"+str(nWALLS)+" 1 "+str(2*numWallVariables)+" \n ")
	wallcounter=1
        wallLength=1
	while wallcounter<=nWALLS:
	    chan.write(str(wallLength)+" ")
            for j in range(0, numWallVariables):
              wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
           
            chan.write("\n")

            wallcounter+=1
        chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
	    
        wallcounter=1
        while wallcounter<=nCells:
	  chan.write( "1 0 1 0.9 ")
          for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
            
          chan.write("\n")
          wallcounter+=1  
        chan.close()







 
def createSquareGridNoise(nGrid,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,wallsSpread):  
    n=nGrid
    global nCells, nWALLS, nVERTEX
    c=0
    nCells= 0
    nWALLS=0
    nVERTEX=0
    chan = open(outputName,"w")
    nCells= n**2
    nWALLS=4+6*(n-1)+2*(n-1)**2
    nVERTEX=(n+1)**2
    ncount=n-1 
    vertexone=0
    vertextwo=1
    ncount3=3*n
    chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   #prints first line of init file
    while c<=ncount: #first row of vertical walls
       chan.write(" "+str(c)+" "+str(c)+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
       c=c+1
       vertexone=vertexone+1
       vertextwo=vertextwo+1
    
    
    cellone=-1
    celltwo=0
    vertexone=n+1
    vertextwo=celltwo
    wallcounter=1
    celloneA=0
    celltwoA=celltwo
    vertexoneA=n+1
    vertextwoA=n+2
    
    
    while wallcounter<=n-1:
     ncount=ncount+n+1
     chan.write(" "+str(c)+" "+str(celltwo)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
     c=c+1
     vertexone=vertexone+1
     vertextwo=vertextwo+1
     celltwo=celltwo+1
     cellone=cellone+1
     while c<=ncount-1: #horizontal walls
       chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
       c=c+1
       vertexone=vertexone+1
       vertextwo=vertextwo+1
       celltwo=celltwo+1
       cellone=cellone+1
     chan.write(" "+str(c)+" "+str(cellone)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
     vertexone=vertexone+1
     vertextwo=vertextwo+1
     ncount=ncount+n   
     c=c+1
     celltwoA=celltwo
    
     while c<=ncount: #internal vertical walls
      chan.write(" "+str(c)+" "+str(celloneA)+" "+str(celltwoA)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
      c=c+1
      vertextwoA+=1
      vertexoneA+=1
      celltwoA+=1
      celloneA+=1
     vertextwoA+=1
     vertexoneA+=1
     wallcounter+=1
    
    ncount=ncount+n+1
    chan.write(" "+str(c)+" "+str(celltwo)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
    c=c+1
    vertexone=vertexone+1
    vertextwo=vertextwo+1
    celltwo=celltwo+1
    cellone=cellone+1
    while c<=ncount-1: #horizontal walls
       chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
       c=c+1
       vertexone=vertexone+1
       vertextwo=vertextwo+1
       celltwo=celltwo+1
       cellone=cellone+1
    chan.write(" "+str(c)+" "+str(cellone)+" "+str(-1)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
    vertexone=vertexone+1
    vertextwo=vertextwo+1
    ncount=ncount+n   
    c=c+1
    
    while c<=nWALLS-1:
      chan.write(" "+str(c)+" "+str(celloneA)+" "+str(-1)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
      c=c+1
      vertextwoA=vertextwoA+1
      vertexoneA=vertexoneA+1
      celltwoA=celltwoA+1
      celloneA=celloneA+1
    
    chan.write(" \n")
    chan.write(" "+str(nVERTEX)+" 3 \n")
    ycoord=0  
    xcoord=0
    while xcoord<=n:
      while ycoord<=n:
        chan.write(""+str(xcoord+0.25*random.uniform(-1, 1))+" "+str(ycoord+0.25*random.uniform(-1, 1))+" 0\n")
        ycoord+=1
      ycoord=0
      xcoord+=1
    
    chan.write("\n"+str(nWALLS)+" 1 "+str(numWallVariables*2)+"\n ")
    wallcounter=1
    wallLength=1
    while wallcounter<=nWALLS:
        chan.write(str(wallLength)+" ")
        for j in range(0, numWallVariables):
            wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
            
        chan.write("\n")

        wallcounter+=1
    chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
    wallcounter=1
    
    while wallcounter<=nCells:
    

        chan.write( "1 0 1 0.9 ")
        for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
            
        chan.write("\n")
        wallcounter+=1  
    chan.close()




def createOffsetSquareGridNoise(nGrid,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,wallsSpread):  
	n=nGrid
	global nCells, nWALLS, nVERTEX
	nCells= 0
	nWALLS=0
	nVERTEX=0


	chan = open(outputName,"w")

	#calculates number of cells, walls and vertices and prints first line of init file
	nCells= n**2
	nWALLS=4*n+(n-1)*(2*n+1)+n*(n+1)
	nVERTEX=2*(2*n+1)+(n-1)*(2*n+2)
	chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   


	#first row of vertical walls
	ncount=2*n-1 
	vertexone=0
	vertextwo=1
	cellTwo=0
	w=0

	while w<=ncount: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1

	#internal row(s) of vertical walls

	vertexone=2*n
	vertextwo=2*n+1
	for i in range(0,n-1,1):
		
		vertexone+=1
	   	vertextwo+=1
		if i==0:
			cellOne=i*n-0.5

		elif (i % 2 == 0): #even 
			cellOne=i*n-0.5
		else:
			cellOne=(i+1)*n-0.5
		if (i % 2 == 0): #even 
			counter=i+1
		else:
			counter=i
		cellTwo=counter*n
		chan.write(" "+str(w)+" "+str(int(cellTwo))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	




		w+=1
		cellTwo+=0.5
		cellOne+=0.5
	  	vertexone+=1
	   	vertextwo+=1
		while w<=(4+2*i)*n-1+i: 



		   chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
		   w+=1
		   cellTwo+=0.5
		   cellOne+=0.5
		   vertexone+=1
	   	   vertextwo+=1

		#internal row(s) of vertical walls
		chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")

		w+=1
		cellTwo+=0.5
		cellOne+=0.5
		vertexone+=1
	   	vertextwo+=1


	#final row of vertical walls
	vertexone+=1
	vertextwo+=1
	cellTwo=n*(n-1)
	while w<=nWALLS-(n*(n+1))-1: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1


	#horizontal walls
	vertexone=0
	vertextwo=2*n+2
	cellTwo=-1

	for i in range(0,n,1):
		
		k=1
		cellTwo=-1
		cellOne=i*n
		while k<=n+1:
			if k==n+1:
				cellOne=-1
			chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
			vertexone+=2
		   	vertextwo+=2
			cellTwo=cellOne
			cellOne+=1

			w+=1
			k+=1
		if (i ==n-2 and i % 2 == 1): #even 
				vertexone+=1
				vertextwo-=0
				
		elif (i ==n-2 and i % 2 == 0): #even 
				vertexone-=1
				vertextwo-=1
			

		elif (i !=n-2 and i % 2 == 0): #even 
			vertexone-=1
			vertextwo-=1
		else:
			vertexone+=1
			vertextwo+=1

	#vertex corrdinates
	chan.write(" \n")
	chan.write(" "+str(nVERTEX)+" 3 \n")
	ycoord=0  
	xcoord=0
	while xcoord<n:
	  
	  while ycoord<=n:
	    chan.write(""+str(xcoord+0.2*random.uniform(-1, 1))+" "+str(ycoord+0.2*random.uniform(-1, 1))+" 0\n")
	    ycoord+=0.5
	  ycoord=-0.5
	  xcoord+=1
	while xcoord<=n:
	  if xcoord==n and n % 2 == 1:
		ycoord=0
		while ycoord<=n:
		    chan.write(""+str(xcoord+0.2*random.uniform(-1, 1))+" "+str(ycoord+0.2*random.uniform(-1, 1))+" 0\n")
		    ycoord+=0.5
	 	ycoord=-0.5
	  	xcoord+=1
	  elif xcoord==n and n % 2 == 0:
		ycoord=-0.5
		while ycoord<n:
	    		chan.write(""+str(xcoord+0.2*random.uniform(-1, 1))+" "+str(ycoord+0.2*random.uniform(-1, 1))+" 0\n")
	    		ycoord+=0.5
	 	ycoord=-0.5
	  	xcoord+=1




	chan.write("\n"+str(nWALLS)+" 1 "+str(2*numWallVariables)+" \n ")
	wallcounter=1
        wallLength=1
	while wallcounter<=nWALLS:
	    chan.write(str(wallLength)+" ")
            for j in range(0, numWallVariables):
              wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
           
            chan.write("\n")

            wallcounter+=1
        chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
	    
        wallcounter=1
        while wallcounter<=nCells:
	  chan.write( "1 0 1 0.9 ")
          for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
            
          chan.write("\n")
          wallcounter+=1  
        chan.close()



def createHexagonalGrid(nGrid,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,wallsSpread):  
	n=nGrid
	global nCells, nWALLS, nVERTEX
	nCells= 0
	nWALLS=0
	nVERTEX=0


	chan = open(outputName,"w")

	#calculates number of cells, walls and vertices and prints first line of init file
	nCells= n**2
	nWALLS=4*n+(n-1)*(2*n+1)+n*(n+1)
	nVERTEX=2*(2*n+1)+(n-1)*(2*n+2)
	chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   


	#first row of vertical walls
	ncount=2*n-1 
	vertexone=0
	vertextwo=1
	cellTwo=0
	w=0

	while w<=ncount: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1

	#internal row(s) of vertical walls

	vertexone=2*n
	vertextwo=2*n+1
	for i in range(0,n-1,1):

		vertexone+=1
	   	vertextwo+=1
		if i==0:
			cellOne=i*n-0.5

		elif (i % 2 == 0): #even 
			cellOne=i*n-0.5
		else:
			cellOne=(i+1)*n-0.5
		if (i % 2 == 0): #even 
			counter=i+1
		else:
			counter=i
		cellTwo=counter*n
		chan.write(" "+str(w)+" "+str(int(cellTwo))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	




		w+=1
		cellTwo+=0.5
		cellOne+=0.5
	  	vertexone+=1
	   	vertextwo+=1
		while w<=(4+2*i)*n-1+i: 



		   chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
		   w+=1
		   cellTwo+=0.5
		   cellOne+=0.5
		   vertexone+=1
	   	   vertextwo+=1

		#internal row(s) of vertical walls
		chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")

		w+=1
		cellTwo+=0.5
		cellOne+=0.5
		vertexone+=1
	   	vertextwo+=1


	#final row of vertical walls
	vertexone+=1
	vertextwo+=1
	cellTwo=n*(n-1)
	while w<=nWALLS-(n*(n+1))-1: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1


	#horizontal walls
	vertexone=0
	vertextwo=2*n+2
	cellTwo=-1

	for i in range(0,n,1):

		k=1
		cellTwo=-1
		cellOne=i*n
		while k<=n+1:
			if k==n+1:
				cellOne=-1
			chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
			vertexone+=2
		   	vertextwo+=2
			cellTwo=cellOne
			cellOne+=1

			w+=1
			k+=1
		if (i ==n-2 and i % 2 == 1): #even 
				vertexone+=1
				vertextwo-=0

		elif (i ==n-2 and i % 2 == 0): #even 
				vertexone-=1
				vertextwo-=1
				
		elif (i !=n-2 and i % 2 == 0): #even 
			vertexone-=1
			vertextwo-=1
		else:
			vertexone+=1
			vertextwo+=1
	#vertex corrdinates
	chan.write(" \n")
	chan.write(" "+str(nVERTEX)+" 3 \n")
	ycoord=0  
	xcoord=0
	counter2=0
	while counter2<n:
	  counter=0

	  while ycoord<=n:
	    chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
	    ycoord+=0.5
            if counter2==0:
		    if counter%2==0:
			xcoord-=0.25
		    else:
			xcoord+=0.25	
		    counter+=1		
            elif counter2%2==0:
		    if counter%2==0:
			xcoord+=0.25
		    else:
			xcoord-=0.25	
		    counter+=1		
            else:
		    if counter%2==0:
			xcoord-=0.25
		    else:
			xcoord+=0.25	
		    counter+=1	
	  ycoord=-0.5


	  #start next column
          if counter2==0:
		xcoord+=1.25

          elif counter2%2==0:
		xcoord+=1.25

	  else:
		  xcoord+=0.75
	  counter2+=1
	while counter2<=(n):
	  if counter2==n and n % 2 == 1:
		ycoord=0
	        counter=0
                xcoord-=0.25
		while ycoord<=n:
		
		    #chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
		    counter=0

	            while ycoord<=n:
		        
	    		chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
	   		ycoord+=0.5

                        if counter2%2==0:
		    		if counter%2==0:
					xcoord-=0.25
		    		else:
					xcoord+=0.25	
		    		counter+=1		
           	        else:
		    	   if counter%2==0:
				xcoord+=0.25
		   	   else:
				xcoord-=0.25	
		    	   counter+=1	

	  	xcoord+=1
                counter2+=1


	  elif counter2==n and n % 2 == 0:
		ycoord=-0.5
		while ycoord<(n-1):
	    	    counter=0
		    ycoord=-0.5
	            while ycoord<=(n-0.5):
		        
	    		chan.write(""+str(xcoord)+" "+str(ycoord)+" 0\n")
	   		ycoord+=0.5

                        if counter2%2==0:
		    		if counter%2==0:
					xcoord+=0.25
		    		else:
					xcoord-=0.25	
		    		counter+=1		
           	        else:
		    	   if counter%2==0:
				xcoord-=0.25
		   	   else:
				xcoord+=0.25	
		    	   counter+=1	

	 	ycoord=-0.5
	  	xcoord+=0.75
                counter2+=1




	chan.write("\n"+str(nWALLS)+" 1 "+str(2*numWallVariables)+" \n ")
	wallcounter=1
        wallLength=1
	while wallcounter<=nWALLS:
	    chan.write(str(wallLength)+" ")
            for j in range(0, numWallVariables):
              wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
           
            chan.write("\n")

            wallcounter+=1
        chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
	    
        wallcounter=1
        while wallcounter<=nCells:
	  chan.write( "1 0 1 0.9 ")
          for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
            
          chan.write("\n")
          wallcounter+=1  
        chan.close()



def createCylinderOfCells(numCells,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread):

        n=numCells
	numCells=n*n
        # Calculate number of walls and vertices and write to file

        numWalls = 3*n+(2*n)*(n-1);
        numVertex = (n+1)*n;
        chan = open(outputName,"w")
        chan.write(str(numCells)+" "+str(numWalls)+" "+str(numVertex)+"\n")
    
	############################################################################
	#
	# Wall/Cell/Vertex connections
	#
	############################################################################


	NoWalls=0

	# Walls along the top
	for j in range(0,n-1):
	 chan.write(str(NoWalls)+" "+str(j)+" -1 "+str(j)+" "+str(j+1)+" \n")
	 NoWalls+=1

	chan.write(str(NoWalls)+" "+str(n-1)+" -1 "+str(n-1)+" 0 \n")
	NoWalls+=1
	
        endValA=n-1
	endValB=endValA+n
	endValC=n

	# second horizontal set of walls
	rangeStart=n
	rangeEnd=2*n-1
	for j in range(rangeStart,rangeEnd):
	 chan.write(str(NoWalls)+" "+str(j-n)+" "+str(j)+" "+str(j)+" "+str(j+1)+" \n")
	 NoWalls+=1
	chan.write(str(NoWalls)+" "+str(endValA)+" "+str(endValB)+" "+str(endValB)+" "+str(endValC)+"\n")
	NoWalls+=1
	# shift of the previous set of wall data
	for i in range(1,n-1):
		rangeStart+=n
		rangeEnd+=n
		endValA+=n
		endValB+=n
		endValC+=n
		for j in range(rangeStart,rangeEnd):
		 chan.write(str(NoWalls)+" "+str(j-n)+" "+str(j)+" "+str(j)+" "+str(j+1)+" \n")
		 NoWalls+=1
		chan.write(str(NoWalls)+" "+str(endValA)+" "+str(endValB)+" "+str(endValB)+" "+str(endValC)+"\n")
		NoWalls+=1
	# walls down the RHS
	for i in range(n*n,n*n+n-1):
	 chan.write(str(NoWalls)+" "+str(i-n)+" -1 "+str(i)+" "+str(i+1)+" \n")
	 NoWalls+=1
	chan.write(str(NoWalls)+" "+str(n*n-1)+" -1 "+str(n*n+n-1)+" "+str(n*n)+" \n")
	NoWalls+=1
	# walls along the joining edge
	for i in range(1,n+1):
	 chan.write(str(NoWalls)+" "+str(n*(i-1))+"  "+str(n*(i-1)+n-1)+" "+str(n*(i-1))+" "+str(n*(i))+" \n")
	 NoWalls+=1
	# second horizontal set of walls
	for i in range(1,n+1):
	 chan.write(str(NoWalls)+" "+str(n*(i-1))+"  "+str(1+(n*(i-1)))+" "+str(1+(n*(i-1)))+" "+str(1+(n*i))+" \n")
	 NoWalls+=1
	# shift of the previous set of wall data
	for k in range(1,n-1):
		endValA=1+(k-1)
		endValB=2+(k-1)
		endValC=n+2+(k-1)
		
		for i in range(1,n+1):
			rangeStart+=n
			rangeEnd+=n

			chan.write(str(NoWalls)+" "+str(endValA)+" "+str(endValB)+" "+str(endValB)+" "+str(endValC)+"\n")
			NoWalls+=1
			endValA+=n
			endValB+=n
			endValC+=n


	
	############################################################################
	#
	# Vertex Coordinates
	#
	############################################################################
        chan.write(str(numVertex)+' 3\n')

	zz=50
	CoorStart=[[0],[5],[zz]]; # Initial Coordinate Point: Top of top circle of points
	theta=(math.pi/180)*(360-(360/n));                                      
	RotMat=np.matrix([[math.cos(theta), -math.sin(theta), 0],[ math.sin(theta), math.cos(theta), 0],[0, 0, 1]])
	Rotated=RotMat*CoorStart
        zCoord=math.sqrt((CoorStart[0]-Rotated[0])**2+(CoorStart[1]-Rotated[1])**2)

	for j in range(1,n+2):
	  x=str(CoorStart[0]).replace('[','').replace(']','')
	  y=str(CoorStart[1]).replace('[','').replace(']','')
	  z=str(CoorStart[2]).replace('[','').replace(']','')
	  chan.write(""+str(x)+" "+str(y)+" "+str(z)+"\n")

	# Gather values for the top ring of coordinates


	  for i in range(1, n):
		                                               
	    theta=(math.pi/180)*(360-(360/n)*i);                                      
	    RotMat=np.matrix([[math.cos(theta), -math.sin(theta), 0],[ math.sin(theta), math.cos(theta), 0],[0, 0, 1]]) 
	    Rotated=RotMat*CoorStart
	    Coor=[CoorStart] 
	    x= str(Rotated[0]).replace('[','').replace(']','')
	    y= str(Rotated[1]).replace('[','').replace(']','')
	    z= str(Rotated[2]).replace('[','').replace(']','')
	    chan.write(""+str(x)+" "+str(y)+" "+str(z)+"\n")

	  zz-=zCoord
	  CoorStart=[[0],[5],[zz]];
	wallLength=1




	############################################################################
	#
	# Wall variables
	#
	############################################################################
	chan.write(str(numWalls)+" "+str(1)+" "+str(2*numWallVariables)+"\n")
	for iw in range(0, numWalls):
	    chan.write(str(wallLength)+" ")
	    for j in range(0, numWallVariables):
		wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
		wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
		chan.write(str(wallsData1)+" "+str(wallsData2)+" ") 
	    chan.write("\n")

	############################################################################
	#
	# Cell variables
	#
	############################################################################
	chan.write("\n"+str(numCells)+" "+str(numCellVariables+4)+"\n")
	for i in range(0,numCells): 
	    chan.write( "1 0 1 0.9 ")
	    for j in range(0, numCellVariables):
		cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
		chan.write(str(cellsData)+" ")
	    chan.write("\n")

def createCylinderOfOfsetCells(nGrid,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread):

	n=nGrid
	global nCells, nWALLS, nVERTEX
	nCells= 0
	nWALLS=0
	nVERTEX=0


	chan = open(outputName,"w")

	#calculates number of cells, walls and vertices and prints first line of init file
	nCells= n**2
	nWALLS=n*(2*n+1)+(n+1)*n
	nVERTEX=2*(2*n+1)+(n-1)*(2*n+2)-2*n
	chan.write(""+str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")   

	############################################################################
	#
	# Wall/Cell/Vertex connections
	#
	############################################################################

	#first row of vertical walls
	ncount=2*n 
	vertexone=0
	vertextwo=1
	cellTwo=0
	w=0	
        cellOne=nCells-n
        chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellOne))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	cellOne+=0.5
	w+=1
	vertexone+=1
	vertextwo+=1

	while w<=ncount-1: 
	   chan.write(" "+str(w)+" "+str(int(cellOne))+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellOne+=0.5
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1

        chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	cellTwo+=0.5
	w+=1
	vertexone+=1
	vertextwo+=1
	
	#internal row(s) of vertical walls

	vertexone=2*n+1
	vertextwo=2*n+2
	for i in range(0,n-1,1):
	        wallTracker=0
	
		vertexone+=1
	   	vertextwo+=1
		if i==0:
			cellOne=i*n-0.5

		elif (i % 2 == 0): #even 
			cellOne=i*n-0.5
		else:
			cellOne=(i+1)*n-0.5
		if (i % 2 == 0): #even 
			counter=i+1
		else:
			counter=i
		cellTwo=counter*n
		chan.write(" "+str(w)+" "+str(int(cellTwo))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	

                wallTracker+=1


		w+=1
		cellTwo+=0.5
		cellOne+=0.5
	  	vertexone+=1
	   	vertextwo+=1
		while wallTracker<=2*n-1: 



		   chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
		   w+=1
		   cellTwo+=0.5
		   cellOne+=0.5
		   vertexone+=1
	   	   vertextwo+=1
                   wallTracker+=1

		#internal row(s) of vertical walls
		chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(-1)+" "+str(vertexone)+" "+str(vertextwo)+"\n")

		w+=1
		cellTwo+=0.5
		cellOne+=0.5
		vertexone+=1
	   	vertextwo+=1

	#final row of vertical walls
	vertexone+=1
	vertextwo+=1
	cellTwo=n*(n-1)
	while w<=nWALLS-(n*(n+1))-1: 
	   chan.write(" "+str(w)+" "+str(-1)+"  "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
	   cellTwo+=0.5
	   w+=1
	   vertexone+=1
	   vertextwo+=1


	#horizontal walls
	vertexone=1
	vertextwo=2*n+3
	cellTwo=-1

	for i in range(0,n,1):
		if i==n-1:
			vertextwo=0
		k=1
		cellTwo=-1
		cellOne=i*n
		while k<=n+1:
			if k==n+1:
				cellOne=-1
			chan.write(" "+str(w)+" "+str(int(cellOne))+" "+str(int(cellTwo))+" "+str(vertexone)+" "+str(vertextwo)+"\n")
			vertexone+=2
		   	vertextwo+=2
			cellTwo=cellOne
			cellOne+=1

			w+=1
			k+=1
		if (i ==n-2 and i % 2 == 1): #even 
				vertexone+=1
				vertextwo-=0
				
		elif (i ==n-2 and i % 2 == 0): #even 
				vertexone-=1
				vertextwo-=1
			

		elif (i !=n-2 and i % 2 == 0): #even 
			vertexone-=1
			vertextwo-=1
		else:
			vertexone+=1
			vertextwo+=1

	############################################################################
	#
	# Vertex Coordinates
	#
	############################################################################
        chan.write(str(nVERTEX)+' 3\n')

	zz=50
	CoorStart=[[0],[0.75],[zz]]; # Initial Coordinate Point: Top of top circle of points
	theta=(math.pi/180)*(360-(360/n));                                      
	RotMat=np.matrix([[math.cos(theta), -math.sin(theta), 0],[ math.sin(theta), math.cos(theta), 0],[0, 0, 1]])
	Rotated=RotMat*CoorStart
        zCoord=math.sqrt((CoorStart[0]-Rotated[0])**2+(CoorStart[1]-Rotated[1])**2)

	for j in range(1,2):
	  x=str(CoorStart[0]).replace('[','').replace(']','')
	  y=str(CoorStart[1]).replace('[','').replace(']','')
	  z=str(CoorStart[2]).replace('[','').replace(']','')
	  chan.write(""+str(x)+" "+str(z)+" "+str(y)+"\n")
     
	# Gather values for the top ring of coordinates
          z=float(z)+(zCoord/2)
          for i in range (0,2*n+1): 
	    chan.write(""+str(x)+" "+str(z)+" "+str(y)+"\n")
            z+=(zCoord/2)
	  for i in range(1, n):
		                                               
	    theta=(math.pi/180)*(360-(360/n)*i);                                      
	    RotMat=np.matrix([[math.cos(theta), -math.sin(theta), 0],[ math.sin(theta), math.cos(theta), 0],[0, 0, 1]]) 
	    Rotated=RotMat*CoorStart
	    Coor=[CoorStart] 
	    x= str(Rotated[0]).replace('[','').replace(']','')
	    y= str(Rotated[1]).replace('[','').replace(']','')
	    z= str(Rotated[2]).replace('[','').replace(']','')
	    chan.write(""+str(x)+" "+str(z)+" "+str(y)+"\n")
            z=float(z)+(zCoord/2)
            for i in range (0,2*n+1): 
	      chan.write(""+str(x)+" "+str(z)+" "+str(y)+"\n")
              z+=zCoord/2

	  zz-=zCoord
	  CoorStart=[[0],[0.75],[zz]];
	wallLength=1






	############################################################################
	#
	# Wall variables
	#
	############################################################################


	chan.write("\n"+str(nWALLS)+" 1 "+str(2*numWallVariables)+" \n ")
	wallcounter=1
        wallLength=1
	while wallcounter<=nWALLS:
	    chan.write(str(wallLength)+" ")
            for j in range(0, numWallVariables):
              wallsData1 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              wallsData2 = wallsMean[j] + wallsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
              chan.write(str(wallsData1)+" "+str(wallsData2)+" ")
           
            chan.write("\n")

            wallcounter+=1

	############################################################################
	#
	# Cell variables
	#
	############################################################################
        chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
	    
        wallcounter=1
        while wallcounter<=nCells:
	  chan.write( "1 0 1 0.9 ")
          for j in range(0, numCellVariables):
            cellsData = cellsMean[j] + cellsSpread[j]*(2.0*random.uniform(0, 1)-1.0);
            chan.write(str(cellsData)+" ")
            
          chan.write("\n")
          wallcounter+=1  
        chan.close()




 
 
def createInit(numCells,gridType,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread):
    if gridType=='circle':
         createCircleOfCells(numCells,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread)   
    elif gridType=='line':
         createLineOfCells(numCells,outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread)
    elif gridType=='squareGrid':
        createSquareGrid(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread)    
    elif gridType=='squareGridNoise':
        createSquareGridNoise(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread)     
    elif gridType=='offsetSquareGrid':
        createOffsetSquareGrid(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread) 
    elif gridType=='offsetSquareGridNoise':
        createOffsetSquareGridNoise(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread) 
    elif gridType=='hexagonalGrid':
        createHexagonalGrid(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread) 
    elif gridType=='cylinderSquare':
        createCylinderOfCells(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread) 
    elif gridType=='cylinderOffsetSquare':
        createCylinderOfOfsetCells(numCells, outputName,numCellVariables,cellsMean,
    cellsSpread,numWallVariables,wallsMean,
    wallsSpread) 



    else:
        print('WARNING: unknown gridType')
         



