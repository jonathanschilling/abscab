import os
import sys
import numpy as np

abscab_path = os.path.abspath("../../../main/python")
if not abscab_path in sys.path:
    sys.path.append(abscab_path)

from abscab import MU_0,                         \
                   straightWireSegment_A_z,      \
                   straightWireSegment_B_phi,    \
                   circularWireLoop_A_phi,       \
                   circularWireLoop_B_rho,       \
                   circularWireLoop_B_z,         \
                   magneticFieldPolygonFilament, \
                   magneticFieldCircularFilament, \
                   magneticFieldPolygonFilamentVertexSupplier

from abscab_util import errorMetric

omega = 0.0

radius = 0.0
def vertexSupplierStd(idxVertex):
    phi = omega * idxVertex
    point = np.zeros(3)
    point[0] = radius * np.cos(phi)
    point[1] = radius * np.sin(phi)
    # point[2] = 0.0
    return point

rCorr = 0.0
def vertexSupplierMcG(idxVertex):
    phi = omega * idxVertex
    point = np.zeros(3)
    point[0] = rCorr * np.cos(phi)
    point[1] = rCorr * np.sin(phi)
    # point[2] = 0.0
    return point

def demoMcGreivy():
    global omega, radius, rCorr

    radius  = 1.23 # m
    current = 17.0 # A

    center = np.array([ 0.0, 0.0, 0.0 ])
    normal = np.array([ 0.0, 0.0, 1.0 ])

    evalPos = np.array([
            [10.0, 5.0, 0.0],
    ])

    magneticField = magneticFieldCircularFilament(center, normal, radius, current, evalPos)
    bZRef = magneticField[0, 2]
    print("ref B_z = %.3e"%(bZRef,))

    # mimic circular wire loop as:
    # a) Polygon with points on the circule to be mimiced
    # b) Polygon with points slightly offset radially outward (McGreivy correction)
    # --> a) should have 2nd-order convergence;
    #     b) should have 4th-order convergence wrt. number of Polygon points

    allNumPhi = [
            10, 30, 100, 300, 1000, 3000,
            10000, 30000, 100000, 300000,
            1000000, 3000000,
            10000000, 30000000,
            100000000, 300000000, 1000000000
    ]

    # allNumPhi = [10, 30, 100, 300, 1000, 3000]
    
    numCases = len(allNumPhi)

    allBzStdErr = np.zeros(numCases)
    allBzMcGErr = np.zeros(numCases)

    resultTable = np.zeros((numCases, 3))

    useCompensatedSummation = True

    for i in range(numCases):

        numPhi = allNumPhi[i]
        print("case %2d/%2d: numPhi = %d"%(i+1, numCases, numPhi))

        omega = 2.0 * np.pi / (numPhi-1)

        magneticField = magneticFieldPolygonFilamentVertexSupplier(numPhi, vertexSupplierStd, current, evalPos, useCompensatedSummation)
        bZStd = magneticField[0, 2]

        # double[][] verticesStd = polygonCircleAround0(radius, numPhi);
        # double bZStd = ABSCAB.magneticFieldPolygonFilament(verticesStd, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

        allBzStdErr[i] = errorMetric(bZRef, bZStd)
        print("ABSCAB B_z = %.3e (err %g)"%(bZStd, allBzStdErr[i]))

        # McGreivy radius correction
        dPhi = 2.0 * np.pi / (numPhi - 1) # spacing between points

        # TODO: understand derivation of alpha for special case of closed circle
        # |dr/ds| = 2*pi
        # --> alpha = 1/R * (dr)^2 / 12
        # == 4 pi^2 / (12 R)
        rCorr = radius * (1.0 + dPhi*dPhi/ 12)

        magneticField = magneticFieldPolygonFilamentVertexSupplier(numPhi, vertexSupplierMcG, current, evalPos, useCompensatedSummation)
        bZMcG = magneticField[0, 2]

        # double[][] verticesMcG = polygonCircleAround0(rCorr, numPhi);
        # double bZMcG = ABSCAB.magneticFieldPolygonFilament(verticesMcG, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

        allBzMcGErr[i] = errorMetric(bZRef, bZMcG)
        print("McGrvy B_z = %.3e (err %g)"%(bZMcG, allBzMcGErr[i]))

        resultTable[i, 0] = numPhi
        resultTable[i, 1] = allBzStdErr[i]
        resultTable[i, 2] = allBzMcGErr[i]

    if useCompensatedSummation:
        np.savetxt("convergenceMcGreivy_CompensatedSummation.dat", resultTable)
    else:
        np.savetxt("convergenceMcGreivy_StandardSummation.dat", resultTable)




if __name__ == "__main__":
    demoMcGreivy()
    
    