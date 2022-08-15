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

for i in A_z_ax \
         A_z_rad \
         A_z_f \
         A_z_n
do
  ./plotStraightWireSegmentVsRef.py \
    ../resources/StraightWireSegment_A_z_ref.dat \
    ../../../data/StraightWireSegment_${i}_Java.dat \
    ../../../more_plots/StraightWireSegment_${i}_Java.pdf #\
#    "${i} of Straight Wire Segment (Java)"
done

for i in B_phi_rad \
         B_phi_f \
         B_phi_n
do
  ./plotStraightWireSegmentVsRef.py \
    ../resources/StraightWireSegment_B_phi_ref.dat \
    ../../../data/StraightWireSegment_${i}_Java.dat \
    ../../../more_plots/StraightWireSegment_${i}_Java.pdf #\
#    "${i} of Straight Wire Segment (Java)"
done

# Circular Wire Loop

for i in A_phi_f \
         A_phi_n \
         A_phi_v
do
  ./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_A_phi_ref.dat \
    ../../../data/CircularWireLoop_${i}_Java.dat \
    ../../../more_plots/CircularWireLoop_${i}_Java.pdf #\
#    "${i} of Circular Wire Loop (Java)"
done

for i in B_rho_f \
         B_rho_n \
         B_rho_v
do
  ./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_B_rho_ref.dat \
    ../../../data/CircularWireLoop_${i}_Java.dat \
    ../../../more_plots/CircularWireLoop_${i}_Java.pdf #\
#    "${i} of Circular Wire Loop (Java)"
done

for i in B_z_f1 \
         B_z_f2 \
         B_z_n \
         B_z_v
do
  ./plotCircularWireLoopVsRef.py \
    ../resources/CircularWireLoop_B_z_ref.dat \
    ../../../data/CircularWireLoop_${i}_Java.dat \
    ../../../more_plots/CircularWireLoop_${i}_Java.pdf #\
#    "${i} of Circular Wire Loop (Java)"
done

popd # back from ../src/test/python
