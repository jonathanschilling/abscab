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
                   circularWireLoop_B_z,      \
                   MU_0, \
                   magneticFieldPolygonFilament
                   
from abscab_util import assertRelAbsEquals

def testStraightWireSegment():
    
    test_points_rp = np.loadtxt("../../resources/testPointsRpStraightWireSegment.dat")
    test_points_zp = np.loadtxt("../../resources/testPointsZpStraightWireSegment.dat")
    
    numCases = len(test_points_rp)
    if len(test_points_zp) != numCases:
        raise RuntimeError("number of test point coordinates needs to agree (is numR=%d, numZ=%d)"%(numCases, len(test_points_zp)))
        
    aZRef   = np.loadtxt("../../resources/StraightWireSegment_A_z_ref.dat")
    bPhiRef = np.loadtxt("../../resources/StraightWireSegment_B_phi_ref.dat")
    
    toleranceAZ   = 1.0e-15
    toleranceBPhi = 1.0e-15
    
    status = 0
    for i in range(numCases):

        rp = test_points_rp[i]
        zp = test_points_zp[i]

        # compute values using implementation to test
        aZ   = straightWireSegment_A_z(rp, zp)
        bPhi = straightWireSegment_B_phi(rp, zp)

        aZStatus = assertRelAbsEquals(aZRef[i], aZ, toleranceAZ)
        if aZStatus != 0:
            print("error: mismatch at Straight Wire Segment A_z test case %d"%(i,));
            print("     rho' = %.17e"%(rp,));
            print("       z' = %.17e"%(zp,));
            print("  ref A_z = %+.17e"%(aZRef[i],));
            print("  act A_z = %+.17e"%(aZ,));
        status |= aZStatus

        bPhiStatus = assertRelAbsEquals(bPhiRef[i], bPhi, toleranceBPhi);
        if bPhiStatus != 0:
            print("error: mismatch at Straight Wire Segment B_phi test case %d"%(i,));
            print("       rho' = %.17e"%(rp,));
            print("         z' = %.17e"%(zp,));
            print("  ref B_phi = %+.17e"%(bPhiRef[i],));
            print("  act B_phi = %+.17e"%(bPhi,));
        status |= bPhiStatus

        if status != 0:
            break

    return status
        
def testCircularWireLoop():
    
    test_points_rp = np.loadtxt("../../resources/testPointsRpCircularWireLoop.dat")
    test_points_zp = np.loadtxt("../../resources/testPointsZpCircularWireLoop.dat")
    
    numCases = len(test_points_rp)
    if len(test_points_zp) != numCases:
        raise RuntimeError("number of test point coordinates needs to agree (is numR=%d, numZ=%d)"%(numCases, len(test_points_zp)))
        
    aPhiRef = np.loadtxt("../../resources/CircularWireLoop_A_phi_ref.dat")
    bRhoRef = np.loadtxt("../../resources/CircularWireLoop_B_rho_ref.dat")
    bZRef   = np.loadtxt("../../resources/CircularWireLoop_B_z_ref.dat")
    
    toleranceAPhi = 1.0e-15
    toleranceBRho = 1.0e-13
    toleranceBZ   = 1.0e-14
    
    status = 0
    for i in range(numCases):

        rp = test_points_rp[i]
        zp = test_points_zp[i]

        # compute values using implementation to test
        aPhi = circularWireLoop_A_phi(rp, zp)
        bRho = circularWireLoop_B_rho(rp, zp)
        bZ   = circularWireLoop_B_z(rp, zp)

        aPhiStatus = assertRelAbsEquals(aPhiRef[i], aPhi, toleranceAPhi)
        if aPhiStatus != 0:
            print("error: mismatch at Circular Wire Loop A_phi test case %d"%(i,));
            print("       rho' = %.17e"%(rp,));
            print("         z' = %.17e"%(zp,));
            print("  ref A_phi = %+.17e"%(aPhiRef[i],));
            print("  act A_phi = %+.17e"%(aPhi,));
        status |= aPhiStatus

        bRhoStatus = assertRelAbsEquals(bRhoRef[i], bRho, toleranceBRho);
        if bRhoStatus != 0:
            print("error: mismatch at Circular Wire Loop B_rho test case %d"%(i,));
            print("       rho' = %.17e"%(rp,));
            print("         z' = %.17e"%(zp,));
            print("  ref B_rho = %+.17e"%(bRhoRef[i],));
            print("  act B_rho = %+.17e"%(bRho,));
        status |= bRhoStatus
        
        bZStatus = assertRelAbsEquals(bZRef[i], bZ, toleranceBZ);
        if bZStatus != 0:
            print("error: mismatch at Circular Wire Loop B_z test case %d"%(i,));
            print("     rho' = %.17e"%(rp,));
            print("       z' = %.17e"%(zp,));
            print("  ref B_z = %+.17e"%(bZRef[i],));
            print("  act B_z = %+.17e"%(bZ,));
        status |= bZStatus

        if status != 0:
            break
        
    return status

def testMagneticFieldInfiniteLineFilament():
    tolerance = 1.0e-15

    # Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
    # B(r) = mu_0 * I / (2 pi r)
    # Test this here with:
    # I = 123.0 A
    # r = 0.132 m
    # => B = 0.186 mT
    current = 123.0
    r = 0.132
    bPhiRef = MU_0 * current / (2.0 * np.pi * r)
    # print("ref bPhi = %.5e"%(bPhiRef,))

    vertices = np.array([
            [0.0, 0.0, -1.0e6],
            [0.0, 0.0,  1.0e6]
    ])

    evalPos = np.array([
            [r, 0.0, 0.0]
    ])

    # y component is B_phi
    magneticField = magneticFieldPolygonFilament(vertices, current, evalPos)
    bPhi = magneticField[0, 1]
    # print("act bPhi = %.5e"%(bPhi,))

    relAbsErr = np.abs(bPhi - bPhiRef) / (1.0 + np.abs(bPhiRef))
    # print("raErr = %.5e"%(relAbsErr,))

    return assertRelAbsEquals(bPhiRef, bPhi, tolerance)

def testBPhiInfiniteLineFilament():
    tolerance = 1.0e-15

    # Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
    # B(r) = mu_0 * I / (2 pi r)
    # Test this here with:
    # I = 123.0 A
    # r = 0.132 m
    # => B = 0.186 mT
    current = 123.0
    r = 0.132
    bPhiRef = MU_0 * current / (2.0 * np.pi * r)
    # print("ref bPhi = %.5e"%(bPhiRef,))

    # half the length of the wire segment
    halfL = 1e6
    L = 2*halfL
    rhoP = r / L
    zP = halfL / L
    bPhi = MU_0 * current / (4.0 * np.pi * L) * straightWireSegment_B_phi(rhoP, zP)
    # print("act bPhi = %.5e"%(bPhi,))

    relAbsErr = np.abs(bPhi - bPhiRef) / (1.0 + np.abs(bPhiRef))
    # print("raErr = %.5e"%(relAbsErr,))

    return assertRelAbsEquals(bPhiRef, bPhi, tolerance)

def testMagneticFieldInsideLongCoil():
    tolerance = 1.0e-4

    # Demtroeder 2, Sec. 3.2.3 ("Magnetic field of a long coil")
    # B_z = mu_0 * n * I
    # where n is the winding density: n = N / L
    # of a coil of N windings over a length L
    # Example (which is tested here):
    # n = 1e3 m^{-1}
    # I = 10 A
    # => B = 0.0126T
    bZRef = 0.0126

    N       = 50000 # windings
    L       = 50.0  # total length of coil in m
    current = 10.0  # A
    radius  = 1.0   # m
    
    n = N/L   # winding density: windings per m

    bZ = 0.0;
    for i in range(N):

        # axial position of coil
        z0 = -L/2.0 + (i + 0.5) / n

        # compute magnetic field
        prefac = MU_0 * current / (np.pi * radius)
        bZContrib = prefac * circularWireLoop_B_z(0.0, z0)

        # print("coil %d at z0 = % .3e => contrib = %.3e"%(i, z0, bZContrib,))

        bZ += bZContrib
    # print("B_z = %.5e"%(bZ,))

    relAbsErr = np.abs(bZ - bZRef) / (1.0 + np.abs(bZRef))
    # print("raErr = %.5e"%(relAbsErr,))

    return assertRelAbsEquals(bZRef, bZ, tolerance)
        
if __name__ == "__main__":
    status = 0
    status |= testStraightWireSegment()
    status |= testCircularWireLoop()
    status |= testMagneticFieldInfiniteLineFilament()
    status |= testBPhiInfiniteLineFilament()
    status |= testMagneticFieldInsideLongCoil()
    sys.exit(status)
    