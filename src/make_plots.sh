#!/bin/bash

# Create plots that show the relative deviation between reference data
# and the outputs from the various demo programs.

# Data files with test knot positions are hardcoded in plot scripts,
# so we need to execute the plot scripts from the correct directory.
pushd ./test/python

#### Java ####

# Straight Wire Segment

./plotStraightWireSegmentVsRef.py \
    ../resources/StraightWireSegment_A_z_ref.dat \
    ../../../data/StraightWireSegment_A_z_Java.dat \
    ../../../article/img/StraightWireSegment_A_z_Java.pdf \
    "\$A_z\$ of Straight Wire Segment (Java)"

./plotStraightWireSegmentVsRef.py \
    ../resources/StraightWireSegment_B_phi_ref.dat \
    ../../../data/StraightWireSegment_B_phi_Java.dat \
    ../../../article/img/StraightWireSegment_B_phi_Java.pdf \
    "\$B_\\varphi\$ of Straight Wire Segment (Java)"

# Circular Wire Loop

./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_A_phi_ref.dat \
    ../../../data/CircularWireLoop_A_phi_Java.dat \
    ../../../article/img/CircularWireLoop_A_phi_Java.pdf \
    "\$A_\\varphi\$ of Circular Wire Loop (Java)"

./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_B_rho_ref.dat \
    ../../../data/CircularWireLoop_B_rho_Java.dat \
    ../../../article/img/CircularWireLoop_B_rho_Java.pdf \
    "\$B_\\rho\$ of Circular Wire Loop (Java)"

./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_B_z_ref.dat \
    ../../../data/CircularWireLoop_B_z_Java.dat \
    ../../../article/img/CircularWireLoop_B_z_Java.pdf \
    "\$B_z\$ of Circular Wire Loop (Java)"

popd # back from ./test/python
