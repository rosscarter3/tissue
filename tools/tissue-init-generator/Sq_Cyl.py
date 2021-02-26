import math
import numpy as np
import random

chan = open("test.init","w")
chan.write("400 820 420 \n")
chan.write(" 0 0 -1 0 1 \n")
NoWalls=1
for j in range(1,19):
 chan.write(str(NoWalls)+" "+str(j)+" -1 "+str(j)+" "+str(j+1)+" \n")
 NoWalls+=1
chan.write(str(NoWalls)+" 19 -1 19 0 \n")
NoWalls+=1
endValA=19
endValB=39
endValC=20
rangeStart=20
rangeEnd=39

for j in range(rangeStart,rangeEnd):
 chan.write(str(NoWalls)+" "+str(j-20)+" "+str(j)+" "+str(j)+" "+str(j+1)+" \n")
 NoWalls+=1
chan.write(str(NoWalls)+" "+str(endValA)+" "+str(endValB)+" "+str(endValB)+" "+str(endValC)+"\n")
NoWalls+=1
for i in range(1,19):
	rangeStart+=20
	rangeEnd+=20
	endValA+=20
	endValB+=20
	endValC+=20
	for j in range(rangeStart,rangeEnd):
	 chan.write(str(NoWalls)+" "+str(j-20)+" "+str(j)+" "+str(j)+" "+str(j+1)+" \n")
	 NoWalls+=1
	chan.write(str(NoWalls)+" "+str(endValA)+" "+str(endValB)+" "+str(endValB)+" "+str(endValC)+"\n")
	NoWalls+=1

for i in range(400,419):
 chan.write(str(NoWalls)+" "+str(i-20)+" -1 "+str(i)+" "+str(i+1)+" \n")
 NoWalls+=1
chan.write(str(NoWalls)+" 399 -1 419 400 \n")
NoWalls+=1
for i in range(1,21):
 chan.write(str(NoWalls)+" "+str(20*(i-1))+"  "+str(20*(i-1)+19)+" "+str(20*(i-1))+" "+str(20*(i))+" \n")
 NoWalls+=1
for i in range(1,21):
 chan.write(str(NoWalls)+" "+str(20*(i-1))+"  "+str(1+20*(i-1))+" "+str(1+20*(i-1))+" "+str(1+20*(i))+" \n")
 NoWalls+=1
for k in range(1,19):
	endValA=1+(k-1)
	endValB=2+(k-1)
	endValC=22+(k-1)
	print NoWalls
	for i in range(1,21):
		rangeStart+=20
		rangeEnd+=20

		chan.write(str(NoWalls)+" "+str(endValA)+" "+str(endValB)+" "+str(endValB)+" "+str(endValC)+"\n")
		NoWalls+=1
		endValA+=20
		endValB+=20
		endValC+=20


chan.write('420 3\n')

zz=10
CoorStart=[[0],[5],[zz]]; # Initial Coordinate Point: Top of top circle of points
for j in range(1,22):
  x=str(CoorStart[0]).replace('[','').replace(']','')
  y=str(CoorStart[1]).replace('[','').replace(']','')
  z=str(CoorStart[2]).replace('[','').replace(']','')
  chan.write(""+str(x)+" "+str(y)+" "+str(z)+"\n")
  for i in range(1, 20):
                                                       #
    theta=(math.pi/180)*(360-18*i);                                       #
    RotMat=np.matrix([[math.cos(theta), -math.sin(theta), 0],[ math.sin(theta), math.cos(theta), 0],[0, 0, 1]]) # Gather values for
    Rotated=RotMat*CoorStart
    Coor=[CoorStart] 
    x= str(Rotated[0]).replace('[','').replace(']','')
    y= str(Rotated[1]).replace('[','').replace(']','')
    z= str(Rotated[2]).replace('[','').replace(']','')
    chan.write(""+str(x)+" "+str(y)+" "+str(z)+"\n")

  zz-=1
  CoorStart=[[0],[5],[zz]];
wallLength=1
numWalls=820
numCells=400
numCellVariables=5
cellsMean=[0.5,0.5,0.5,0,0]
cellsSpread=[0.02,0.02,0.02,0,0]
numWallVariables=2
wallsMean=[0, 0]
wallsSpread=[0,0]

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
    
