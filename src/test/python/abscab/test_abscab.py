import os
import sys
import numpy as np

abscab_path = os.path.abspath("../../../main/python")
if not abscab_path in sys.path:
    sys.path.append(abscab_path)

from abscab import straightWireSegment_A_z,   \
                   straightWireSegment_B_phi, \
                   circularWireLoop_A_phi,    \
                   circularWireLoop_B_rho,    \
                   circularWireLoop_B_z


