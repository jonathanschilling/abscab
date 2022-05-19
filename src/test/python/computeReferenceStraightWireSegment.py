#!/usr/bin/env python3

import os
import sys
import mpmath as mp

# default number of digits of precision for calculations
mp.mp.dps = 300

# frequently-used constants with arbitrary precision
zero      = mp.mpf("0")
one       = mp.mpf("1")
two       = mp.mpf("2")
three     = mp.mpf("3")
four      = mp.mpf("4")
half      = one/two
threeHalf = three/two
pi        = mp.pi

# constants used in the definition of IEEE754 double precision variables
exponentBias  = mp.mpf("1023")
mantissaScale = mp.power(two, mp.mpf("52"))

def Ri(rp, zp):
    return mp.sqrt(mp.power(rp, two) + mp.power(zp, two))

def Rf(rp, zp):
    return mp.sqrt(mp.power(rp, two) + mp.power(one - zp, two))

def swlEps(rp, zp):
    return one / (Ri(rp, zp) + Rf(rp, zp))

def A_z(rp, zp):
    return mp.atanh(swlEps(rp, zp))

def B_phi(rp, zp):
    ri = Ri(rp, zp)
    rf = Rf(rp, zp)
    ri_p_rf = ri + rf
    return ri_p_rf / (rf * (mp.power(ri_p_rf, two) - one))

def ieee754_to_arb(s, E, M):

    # capture exact zero
    if E == 0 and M == 0:
        return zero

    # parse sign bit
    # 0 means positive, 1 means negative, as in (-1)^s
    sign = None
    if s == 0:
        sign = one
    else:
        sign = -one

    # compute exact representation of double precision variable
    return sign * mp.power(two, E - exponentBias) * (one + M/mantissaScale)
    
if __name__ == "__main__":

    testPointsFile = "../resources/testPointsStraightWireSegment.dat"
    outFilenameAZ   = "../resources/StraightWireSegment_A_z_ref.dat"
    outFilenameBPhi = "../resources/StraightWireSegment_B_phi_ref.dat"

    with open(testPointsFile, "r") as f:
        lines = f.readlines()

        with open(outFilenameAZ, "w") as outFile:
            outFile.write("")
        with open(outFilenameBPhi, "w") as outFile:
            outFile.write("")

        numLines = len(lines)
        
        for i,line in enumerate(lines):

            # skip comment lines
            if line[0] == "#":
                continue

            if (i+1)%100 == 0:
                print("line %d / %d"%(i+1, numLines))
            
            parts = line.strip().split()

            # construct rp from sign, exponent and mantissa
            signRp     = int(parts[0].strip())
            exponentRp = int(parts[1].strip())
            mantissaRp = int(parts[2].strip())
            rp = ieee754_to_arb(signRp, exponentRp, mantissaRp)

            # construct zp from sign, exponent and mantissa
            signZp     = int(parts[3].strip())
            exponentZp = int(parts[4].strip())
            mantissaZp = int(parts[5].strip())
            zp = ieee754_to_arb(signZp, exponentZp, mantissaZp)        

            # compute magnetostatic quantities: A_z, B_phi
            aZ   = A_z(rp, zp)
            bPhi = B_phi(rp, zp)

            with mp.workdps(20):
                with open(outFilenameAZ, "a") as outFile:
                    outFile.write(str(aZ) + "\n")
                with open(outFilenameBPhi, "a") as outFile:
                    outFile.write(str(bPhi) + "\n")

    sys.exit(0)
