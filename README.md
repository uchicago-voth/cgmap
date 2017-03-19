## CGMap: a simple tool to apply linear maps to molecular trajectories

[![Build Status](https://travis-ci.org/alekepd/cgmap.svg?branch=master)](https://travis-ci.org/alekepd/cgmap)

CGMap provides a tool to map fine-grained (FG) molecular dynamics trajectories
to coarse-grained (CG) trajectories. The map between FG and CG trajectories are
constrainted to be linear, with weights defined as either corresponding to
center of points, center of mass, or center of charge.

This code is in beta. Plese report bugs on the github bug tracker.

Periodic boundary issues are not currently addressed. Generally, wrapping
keeping molecules whole, mapping, and then remapping is the best course of
action.

See the test directory for example usage. No command line interface is present,
nor planned.
