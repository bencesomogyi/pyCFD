pyCFD
=====

This code and its documentation is a result of solving the programming example for the lecture Numerical Methods in Fluid Mechanics and Heat Transfer (LV-Nr.: 321.023 VO, WS 2012/13) @ Graz university of Technology.

About the code
==============

Prior to solving the example the author has set a list of criterion about what
and how the solution should serve learning and self developement. These
criterion resulted in the following list of features:

* a general purpose 3D code/library was written

* colocated grid arrangement was used

* code was written in python in an object oriented fashion

* bottlenecks of the code were re-implemented in C via Cython

* the library supports writing the mathematical operators in the governing
  equations in a style ispired by OpenFOAM

* the output is provided in the file formats of the vtk library allowing the use
  of Paraview for post-processing
