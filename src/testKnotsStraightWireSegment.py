#!/usr/bin/env python

import numpy as np

# machine precision: ca. 2.22e-16 for 64-bit double precision
eps = np.finfo(np.float64).eps

print("machine precision: ", eps)

###################################################
# assemble list of test knots in radial direction #
###################################################
testKnotsRp = []

# 0
testKnotsRp.append(0.0)

# 1e-30 ... 1e30
for exponent in range(-30, 30+1):
  testKnotsRp.append(np.power(10.0, exponent))

print("number of test knots in  radial  direction: ",
      len(testKnotsRp))

#####################################################
# assemble list of test knots in vertical direction #
#####################################################
testKnotsZp = []

# -1e30 ... -1e-30
for exponent in range(30, -30-1, -1):
  testKnotsZp.append(-np.power(10.0, exponent))

# 0
testKnotsZp.append(0.0)

# 1e-30 ... 1e-1
for exponent in range(-30, 0):
  testKnotsZp.append(np.power(10.0, exponent))

# 1/2
testKnotsZp.append(0.5)

# 1 - 1e-1 ... 1 - 1e-15
for exponent in range(-1, -15-1, -1):
  testKnotsZp.append(1 - np.power(10.0, exponent))

# 1 -  eps/2 (next lower double precision number from 1)
testKnotsZp.append(1.0 - eps/2.0)

# 1
testKnotsZp.append(1.0)

# 1 + eps (next higher double precision number from 1)
testKnotsZp.append(1.0 + eps)

# 1 + 1e-15 ... 1 + 1e-1
for exponent in range(-15, -1+1):
  testKnotsZp.append(1 + np.power(10.0, exponent))

# 2
testKnotsZp.append(2.0)

# 1e1 ... 1e30
for exponent in range(1, 30+1):
  testKnotsZp.append(np.power(10.0, exponent))

print("number of test knots in vertical direction: ",
      len(testKnotsZp))

######################################################################
# assemble test points from knots; exclude locations on wire segment #
######################################################################








