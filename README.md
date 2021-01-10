# yosoku

Creates a Miller set given spacegroup, unit cell and resolution limit. 
A sequence or number of scatterers is taken as input - if a sequence is given then the suggested number of molecules in the asu is calculated using Matthews/Rupp.
Using the above information, the number of reflections per scatterer is calculated and a prediction of S phasing likelyhood is given. A graph is produced to indicate how phasing can be achieved by increasing the resolution.

Dependencies: cctbx, matplotlib, numpy, scipy, warnings

Inputs: spacegroup - number or text, unit cell - a b c al be ga (commas not important), asu per mol - number (suggestion given if sequence supplied), resolution - interger or float, sequence or scatterers - sequence as string or number.
