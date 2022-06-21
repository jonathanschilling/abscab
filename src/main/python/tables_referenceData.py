#!/usr/bin/env python3

# Extract a suitable subset of test points and the corresponding reference data
# from the tables provided along with the ABSCAB source code.
# These tables are included at the end of the ABSCAB article.

import numpy as np

dataFolder = "../../test/resources/"

####### Straight Wire Segment #######

# rho', z' knots to include in reference tables for straight wire segment
swsRefSetRp    = [0.0,   1.0e-30,      1.0e-15,      1.0,   1.0e15,      1.0e30     ]
swsRefSetRpTeX = ["$0$", "$10^{-30}$", "$10^{-15}$", "$1$", "$10^{15}$", "$10^{30}$"]

swsRefSetZp    = [-1.0e30, -1.0e15, -1.0, -1.0e-15, -1.0e-30,
                  0.0, 1.0e-30, 1.0e-15, 0.5, 1.0, 2.0, 1.0e15, 1.0e30]
swsRefSetZpTeX = ["$-10^{30}$", "$-10^{15}$", "$-1$", "$-10^{-15}$", "$-10^{-30}$",
                  "$0$", "$10^{-30}$", "$10^{-15}$", "$1/2$", "$1$", "$2$", "$10^{15}$", "$10^{30}$"]

# s, E, M for rho' and z'
swsTestPoints = np.loadtxt(dataFolder+"testPointsStraightWireSegment.dat", dtype=int)

# %.17e for rho' and z'
swsTestPointsRP = np.loadtxt(dataFolder+"testPointsRpStraightWireSegment.dat")
swsTestPointsZP = np.loadtxt(dataFolder+"testPointsZpStraightWireSegment.dat")

numCases = len(swsTestPointsRP)

texLines = []

swsNumRefCases = 0
for i in range(numCases):

    rP = swsTestPointsRP[i]
    zP = swsTestPointsZP[i]

    if rP in swsRefSetRp and zP in swsRefSetZp and (rP > 0.0 or zP < 0.0 or zP > 1.0):

        sR = swsTestPoints[i, 0]
        eR = swsTestPoints[i, 1]
        mR = swsTestPoints[i, 2]

        sZ = swsTestPoints[i, 3]
        eZ = swsTestPoints[i, 4]
        mZ = swsTestPoints[i, 5]

##        print("%2d %+.3e %+.3e %d %4d %16d %d %4d %16d" %
##              (swsNumRefCases, rP, zP, sR, eR, mR, sZ, eZ, mZ))

        jR = swsRefSetRp.index(rP)
        jZ = swsRefSetZp.index(zP)

        texLine = "    %2d & %10s & %10s & %d & %4d & %16d & %d & %4d & %16d \\\\" % \
                  (swsNumRefCases, swsRefSetRpTeX[jR], swsRefSetZpTeX[jZ],
                   sR, eR, mR, sZ, eZ, mZ)
        print(texLine)
                  
        texLines.append(texLine)

        swsNumRefCases += 1
    
print("number of reference cases for SWS: %d"%(swsNumRefCases,))

with open("", "w") as texFile:
    for texLine in texLines:
        texFile.write(texLine + "\n")

####### Circular Wire Loop #######

# rho', z' knots to include in reference tables for circular wire loop
cwlRefSetRp = [0.0, 1.0e-30, 1.0e-15, 0.5, 1.0, 2.0, 1.0e15, 1.0e30]
cwlRefSetZp = [0.0, 1.0e-30, 1.0e-15, 1.0, 1.0e15, 1.0e30]
