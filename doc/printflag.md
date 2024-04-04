# Description of the outputs specified by the printFlag option provided in the solver file

## Introduction

This page describes the output specified by the printFlag parameter in the solver file. The prinFlag is an index specifying the output for a simulation and is provided in addition to the number of time points to print the variables.

## Special printFlag indices

There are 9 'standard' printFlag  values implemented. These are:

<p><b>(0)</b> Standard output for openGL developed plotting of cell and wall variables </p>
<p><b>(1)</b> Standard vtu output assuming single wall compartment for wall variables </p>
<p><b>(2)</b> Standard vtu output assuming two wall compartment for wall variables (except for initial length) </p>
<p><b>(3)</b> Standard output for openGL developed plotting of cell variables </p>
<p><b>(4)</b> Standard output for openGL developed plotting of wall variables </p>
<p><b>(5)</b> Output that can be used for plotting in gnuplot. </p>
<p><b>(6)</b> PLY format </p>
<p><b>(7)</b> PLY format assuming centraltriangulation </p>
<p><b>(10)</b> As (2), and in addition a file is generated (tissue.idata) storing the states in init format.
The time points are divided by a line '#tCount = value' to find individual time points in file.
In addition there are several methods for plotting also membrane data (e.g. PIN1). </p>

## AdHoc printFlag indices

After the standard printFlag indices described above, a number of more adhoc printing options are provided. In this section some useful ones are mentioned: 

<p><b>(26)</b> Standard vtu output assuming single wall compartment for wall variables and in addition printing cell variable values for statistical analysis (or gnuplot plotting) in the file tissue.gdata .</p>

For knowing some of this, the source code BaseSolver::Print() might need to be read.


## Implementing new printFlag options

To implement a new way of printing specific data during a simulation, this can be done by editing the BaseSolver::print(). Add your new prinFlag index in the list and potentially describe it in the list above.



