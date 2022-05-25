#!/bin/bash

# Create plots that show the relative deviation between reference data
# and the outputs from the Java implementation of ABSCAB.
# All internal methods are evaluated at all test points.
# This is done to illustrate the reasoning that went into
# switching between different implementations.

# Data files with test knot positions are hardcoded in plot scripts,
# so we need to execute the plot scripts from the correct directory.
pushd ../src/test/python || exit -1

#### Java ####

# Straight Wire Segment

for i in A_z_along_rhoP_0 \
         A_z_along_zP_0_or_1 \
         A_z_6a \
         A_z_6b \
         A_z_6c \
         A_z_1
do
  ./plotStraightWireSegmentVsRef.py \
    ../resources/StraightWireSegment_A_z_ref.dat \
    ../../../data/StraightWireSegment_${i}_Java.dat \
    ../../../more_plots/StraightWireSegment_${i}_Java.pdf #\
    #"${i} of Straight Wire Segment (Java)"
done

for i in B_phi_3 \
         B_phi_4 \
         B_phi_5
do
  ./plotStraightWireSegmentVsRef.py \
    ../resources/StraightWireSegment_B_phi_ref.dat \
    ../../../data/StraightWireSegment_${i}_Java.dat \
    ../../../more_plots/StraightWireSegment_${i}_Java.pdf #\
    #"${i} of Straight Wire Segment (Java)"
done

# Circular Wire Loop

for i in A_phi_1 \
         A_phi_6 \
         A_phi_5
do
  ./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_A_phi_ref.dat \
    ../../../data/CircularWireLoop_${i}_Java.dat \
    ../../../more_plots/CircularWireLoop_${i}_Java.pdf #\
    #"${i} of Circular Wire Loop (Java)"
done

for i in B_rho_3 \
         B_rho_1 \
         B_rho_4
do
  ./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_B_rho_ref.dat \
    ../../../data/CircularWireLoop_${i}_Java.dat \
    ../../../more_plots/CircularWireLoop_${i}_Java.pdf #\
    #"${i} of Circular Wire Loop (Java)"
done

for i in B_z_1 \
         B_z_2 \
         B_z_4 \
         B_z_5 \
         B_z_6
do
  ./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_B_z_ref.dat \
    ../../../data/CircularWireLoop_${i}_Java.dat \
    ../../../more_plots/CircularWireLoop_${i}_Java.pdf #\
    #"${i} of Circular Wire Loop (Java)"
done

popd # back from ../src/test/python
