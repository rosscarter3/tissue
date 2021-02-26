# -*- coding: utf-8 -*-
"""
Simple 2.5D cell wall and nucleus. Two shapes of the cell Rectangle and Ellipse
while the nuclei always is shaped as an ellipse.

In the main function (and the generating functions for vertex positions)
size parameters (Lx,Ly for cell and Dx,Dy for nucleus) and number vertices
(equal numbers for wall and nuclear membrane) can be set. If rectangular
shape is used this should be a factor of 4.

If springs are to be connected from wall to nucleus vertices, they are aligned in
space so conections between {0,numV-1} to {numV,2*numV-1} gives an even spread.

Author: henrik.jonsson@slcu.cam.ac.uk
"""
import random
import math
import numpy as np

def createSquareCellNucleus(numVertexPerLayer,Lx,Ly,Dx,Dy,Dimension=3):
    # ---
    # Initiate some standard numbers (temporary until used as arguments for function)
    # ---
    Rx = 0.5*Dx
    Ry = 0.5*Dy
    if (numVertexPerLayer%4 != 0) :
        print("createSquareCellNucleus: Error, number of vertices not factor of 4.")
        exit()
    numVertexSide = numVertexPerLayer/4 
    dx = Lx/(numVertexSide)
    dy = Ly/(numVertexSide)

    numVertex = 2*numVertexPerLayer
    numWall = numVertex
    numCell = 2 #cytosol and nucleus
    
    # Create structure, layer by layer
    #maxZI = numZ-1
    #for zI in range(0,maxZI):
    # z = zI*dZ
    z = 0.0
    
    #
    # Create vertex positions for cell (rectangle)
    #
    # Start from bottom left (at 0,0)
    vertexIndex = 0
    x = 0.0
    y = 0.0
    print(str(x)+" "+str(y)+" "+str(z))
    #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
    vertexIndex = vertexIndex+1
    if (numVertexSide>1) :
        for dVI in range(2,numVertexSide+1):
            # x = x
            y = y + dy
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            vertexIndex = vertexIndex+1
    # From top left
    x = 0.0
    y = Ly
    print(str(x)+" "+str(y)+" "+str(z))
    #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
    vertexIndex = vertexIndex+1
    if (numVertexSide>1) :
        for dVI in range(2,numVertexSide+1):
            x = x + dx
            # y = y
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            vertexIndex = vertexIndex+1
    # From top right
    x = Lx
    y = Ly
    print(str(x)+" "+str(y)+" "+str(z))
    #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
    vertexIndex = vertexIndex+1
    if (numVertexSide>1) :
        for dVI in range(2,numVertexSide+1):
            # x = x
            y = y - dy
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            vertexIndex = vertexIndex+1
    # From bottom right
    x = Lx
    y = 0.0
    print(str(x)+" "+str(y)+" "+str(z))
    #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
    vertexIndex = vertexIndex+1
    if (numVertexSide>1) :
        for dVI in range(2,numVertexSide+1):
            x = x - dx
            # y = y
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            vertexIndex = vertexIndex+1
    #print("\n")
    #        
    # Create vertex positions for nucleus (ellipse)
    #
    oX = 0.5*Lx
    oY = 0.5*Ly
    PI = 3.14159265
    #r = numCells*wallLength/(2.0*PI)
    deltaPhi = -2.0*PI/numVertexPerLayer
    phi=5.0*PI/4.0; # bottom left corner
    for vertexIndex in range(numVertexPerLayer,numVertex): 
            x = oX + Rx*math.cos(phi)
            y = oY + Ry*math.sin(phi)
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            phi = phi + deltaPhi
    print("\n")

def createEllipseCellNucleus(numVertexPerLayer,Lx,Ly,Dx,Dy,Dimension=3):
    # ---
    # Initiate some standard numbers (temporary until used as arguments for function)
    # ---
    Rx = 0.5*Dx # x-radius for nucleus
    Ry = 0.5*Dy # y-radius for nucleus
    #numZ = 1 # sets the z-extension of the object
    #dZ = 1.0 # oly relevant if numZ>1
    
    #numVertexPerLayer = 30 # cell vertices and same number of nucleus ones
    numVertex = 2*numVertexPerLayer
    numWall = numVertex
    numCell = 2 #cytosol and nucleus
    
    # Create structure, layer by layer
    #maxZI = numZ-1
    #for zI in range(0,maxZI):
    # z = zI*dZ
    z = 0.0
    
    # ---        
    # Create vertex positions for cell (ellipse)
    # ---
    oX = 0.0
    oY = 0.0
    PI = 3.14159265
    deltaPhi = 2.0*PI/numVertexPerLayer
    phi=0.0
    for vertexIndex in range(0,numVertexPerLayer): 
            x = oX + Lx*math.cos(phi)
            y = oY + Ly*math.sin(phi)
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            phi = phi + deltaPhi
            
    #        
    # Create vertex positions for nucleus (ellipse)
    #
    oX = 0.0
    oY = 0.0
    PI = 3.14159265
    #r = numCells*wallLength/(2.0*PI)
    deltaPhi = 2.0*PI/numVertexPerLayer
    phi=0.0
    for vertexIndex in range(numVertexPerLayer,numVertex): 
            x = oX + Rx*math.cos(phi)
            y = oY + Ry*math.sin(phi)
            print(str(x)+" "+str(y)+" "+str(z))
            #print(str(vertexIndex)+" "+str(x)+" "+str(y)+" "+str(z))
            phi = phi + deltaPhi
    print("\n")

# Print the edge connection data (first part of init file)
# 
def printEdgeConnections(numWall):
    layerNum = numWall/2
    # cytosol edges (cytosol,background) = (0,-1)
    for wallIndex in range(0,layerNum):
        nextVertexIndex=(wallIndex+1)%layerNum
        print(str(wallIndex)+" "+str(0)+" "+str(-1)+" "+str(wallIndex)+" "+str(nextVertexIndex))
    # nucleus edges (cytosol,background) = (0,1)
    for wallIndex in range(layerNum,numWall):
        nextVertexIndex=(wallIndex+1)
        if (nextVertexIndex>=numWall):
            nextVertexIndex = layerNum
        # two last ones swapped to get reversed order of listed vertex indices
        # to sort in opposite direction
        # compared to the cell wall edges (still under test)
        #print(str(wallIndex)+" "+str(1)+" "+str(0)+" "+str(nextVertexIndex)+" "+str(wallIndex))
        # version connecting both cytosol and nucleus edges to background
        # (instead of nuclear membrane to nucleus and cytosol)
        print(str(wallIndex)+" "+str(1)+" "+str(-1)+" "+str(nextVertexIndex)+" "+str(wallIndex))
    print("\n")
    
# Print edge length and variables, flag for wall/nuclear
# right now with all lengths=1, so need to run with reaction initiateWallLength
#
def printEdgeData(numWall):
    layerNum = numWall/2
    numLength = 1
    numVar = 4
    print(str(numWall)+" "+str(numLength)+" "+str(numVar)+" # numWall numL numVar")
    # Cell Wall edges
    for wallIndex in range(0,layerNum):
        print("1.0 1 0 0 0")
    # Nuclear edges
    for wallIndex in range(layerNum,numWall):
        print("1.0 0 1 0 0")
    print("\n")

# Print cell variables (cytosol and nucleus)
# First two flags/vars seperating cytosol and nucleus
def printCellData(numCell):
    numVar = 4
    print(str(numCell)+" "+str(numVar)+" # numCell numVar")
    # Cytosol
    print("1 0 0 0")
    # Nucleus
    print("0 1 0 0")
    print("\n")
    
def main():

    Dimension = 3 # If 3, all z positions are 0.0 (currently only working dimension) 
    Lx = 2.0 # x-size for cell
    Ly = 2.0 # y-size for cell
    Dx = 1.0 # x-diagonal for nucleus
    Dy = 1.0 # y-diagonal for nucleus
    numVertexPerLayer = 8 # number of vertices for cell walls and for nucleus
    numVertex = numWall = 2*numVertexPerLayer
    numCell = 2 # Always 2, index 0 -> cytosol, index 1 -> nucleus
    
    # Print some init comments
    print("#")
    print("# Generated by cellPlusNucleus.py script")
    print("# See gitlab.com/slcu/teamhj/laura/tissue-init-generator")
    print("#")
    
    # Print size values for init file
    print(str(numCell)+" "+str(numWall)+" "+str(numVertex)+" # numCell numWall numVertex")

    # Print connections between edges and cells and vertices
    # Assumes two layers of (wall/nucleus) connected edges 
    printEdgeConnections(numWall)

    # Create and print vertex positions
    print(str(numVertex)+" "+str(Dimension)+" # numVertex numDimension")
    #createSquareCellNucleus(numVertexPerLayer,Lx,Ly,Dx,Dy)
    createEllipseCellNucleus(numVertexPerLayer,Lx,Ly,Dx,Dy)

    # Print Edge data
    printEdgeData(numWall)

    # Print Cell data
    printCellData(numCell)
    
if __name__== "__main__":
    main()
