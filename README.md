## CGMap: a simple tool to apply linear maps to molecular trajectories

![Build Status](https://travis-ci.org/uchicago-voth/cgmap.svg?branch=master)

CGMap provides a tool to map fine-grained (FG) molecular dynamics trajectories
to coarse-grained (CG) trajectories. The map between FG and CG trajectories are
constrained to be linear, with weights defined as either corresponding to
center of points, center of mass, or center of charge.

This code is in beta. Please report bugs on the GitHub bug tracker. Regression
tests are provided. Note that numerical differences in underlying libraries
(through numpy) often create slightly different outputs.

Periodic boundary issues are not currently validated. Generally, wrapping
keeping molecules whole, mapping, and then remapping is the best course of
action, although some effort in the code is made to automatically resolve this.

See the test directory for many cases of example usage. No command line
interface is present, nor planned.
