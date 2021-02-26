# Tissue Init Generator

Code for generating lines or grids of cells for use in tissue simulations

## Types of init file

* circle 
* line
* square grid
* offset square grid
* square grid with noise in vertex positions
* offset square grid with noise in vertex positions


## How to run
`createInit(numCells,gridType,outputName,numCellVariables,cellsMean,cellsSpread,
numWallVariables,wallsMean,wallsSpread)`

Input Variables: 

* `numCells` - numCells by numCells grid, or number of cells in line or circle. 
* `gridType` - can be one of the following: 'circle', 'line', 'squareGrid', 
'squareGridNoise', 'offsetSquareGrid', 'offsetSquareGridNoise','hexagonalGrid'
* `outputName` - the name of the file for the output eg 'circle5.init'
* `numCellVariables` - Number of variables in each cell
* `cellsMean` - Mean value for each cell variable eg `cellsMean=[1,1,1,1]`. Note need to have len(cellsMean)=numCellVariables 
* `cellsSpread` - for adding noise to cell variables 
* `numWallVariables`- Number of variables in each wall
* `wallsMean`- Mean value for each wall variable eg `wallsMean=[1,1,1,1]`. Note need to have len(wallsMean)=numWallVariables 
* `wallsSpread` - for adding noise to cell variables 

## Hypocotyl Geometry
matalb script to generate hypocotyl geometry files, by Behruz, July 2017


## Authors

* **Laura Brown**  [GitLab](https://gitlab.com/slcu/teamHJ/laura)


