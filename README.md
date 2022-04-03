# OpenDLO

Open Discontinuity Layout Optimisation - C++ static library.

[Download the test application](https://github.com/reniercloete/OpenDLO/releases/download/1.0.0/OpenDLO_1.0.0.zip) to test it.

## Introduction
OpenDLO is C++ implementation of discontinuity layout optimisation (DLO) and can be used to estimate the ultimate load carrying capacity of plates subjected to arbitrary out-of-plane loads.  The library currently only supports convex outlines without openings, fixed, free and simply supported edges as well as uniformly distributed loads.  Planned work includes arbitrary outlines with openings, internal supports, point loads and line loads.

OpenDLO could form the basis of software for designing structural steel connections, concrete slabs and masonry wall panels.

## Requirements
OpenDLO can make use of either the [Coin-OR Linear Programming (CLP) library](https://www.coin-or.org/Tarballs/Clp/Clp-1.17.6.zip) or [Mosek](https://www.mosek.com/).  Mosek is not open source, but trial and academic licences are available.  OpenDLO uses a version of CLP that is accelerated with Intel MKL.  The test application makes use of [GLFW](https://www.glfw.org/).

## Test application

The test application defines a 1x1 isotropic square plate with fixed edges and a unit UDL.  With a 0.25 interval discritization (4 nodes per side) the OpenDLO calculates a load factor of 44.24 (analytically exact is 42.851) and a yield line layout as shown below.

![alt text](https://user-images.githubusercontent.com/95902170/161394313-699904d6-f258-4e51-a9a3-313e29a5f9c0.jpeg)

OpenDL was implemented using: Gilbert M., He L., Smith C.C., Le C.V. Automatic yield-line analysis of slabs using discontinuity layout optimization Proc. R. Soc. A, 470 (2168) (2014), p. 20140071
