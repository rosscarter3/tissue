# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 15:55:53 2017

@author: laura
"""
from InitGenerator import *

numCellVariables=4
cellsMean=[1,2,3,4]
cellsSpread=[0.02,0.02,0.02,0.02]
numWallVariables=2
wallsMean=[0, 0]
wallsSpread=[0,0]




gridType='hexagonalGrid'
numCells=7
outputName='hex.init'
createInit(numCells,gridType,outputName,
           numCellVariables,cellsMean,cellsSpread,numWallVariables,
           wallsMean,wallsSpread)


