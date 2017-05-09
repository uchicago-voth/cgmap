## CGMap: a simple tool to apply linear maps to molecular trajectories

![Build Status](https://travis-ci.org/uchicago-voth/cgmap.svg?branch=master)

CGMap provides a tool to map fine-grained (FG) molecular dynamics trajectories
to coarse-grained (CG) trajectories. The map between FG and CG trajectories is
constrained to be linear and time invariant, with weights defined as either
corresponding to center of points, center of mass, or center of charge, although
the user is able to modify this in a manual way (see test cases). The code is
able to automatically detect and modify the contributing weight of atoms shared
between CG sites in a single molecule if the appropriate options are enabled
(see test cases).

This code is in beta. Please report bugs on the GitHub bug tracker. Regression
tests are provided. Note that numerical differences in underlying libraries
(through numpy) often create slightly different outputs.

Periodic boundary issues are handled: atoms are moved to the smallest local
image before mapping, and moved back to the image where the most contributing
atoms were present. Alternatively, options are available to atomatically wrap 
the CG representation after it is produced in the local minimal image.

See the test directory for many cases of example usage. No command line
interface is present, nor planned.
