#!/bin/bash

python plot_CircularWireLoop_results.py \
  ../../test/resources/CircularWireLoop_A_phi_ref.dat \
  ../../../data/CircularWireLoop_A_phi_TEAL.dat \
  ../../../article/img/CircularWireLoop_A_phi_TEAL.pdf \
  "CWL: $\\tilde{A}_\\varphi$"

python plot_CircularWireLoop_results.py \
  ../../test/resources/CircularWireLoop_B_rho_ref.dat \
  ../../../data/CircularWireLoop_B_rho_TEAL.dat \
  ../../../article/img/CircularWireLoop_B_rho_TEAL.pdf \
  "CWL: $\\tilde{B}_\\rho$"

python plot_CircularWireLoop_results.py \
  ../../test/resources/CircularWireLoop_B_z_ref.dat \
  ../../../data/CircularWireLoop_B_z_TEAL.dat \
  ../../../article/img/CircularWireLoop_B_z_TEAL.pdf \
  "CWL: $\\tilde{B}_z$"

