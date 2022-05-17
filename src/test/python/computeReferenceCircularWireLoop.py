#!/usr/bin/env python

import os
import sys
import mpmath as mp

# default number of digits of precision for calculations
mp.mp.dps = 100

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

def kSq(rp, zp):
    return (four * rp) / \
           (mp.power(zp, two) + mp.power(one + rp, two))

def kCSq(rp, zp):
    return (mp.power(zp, two) + mp.power(one - rp, two)) / \
           (mp.power(zp, two) + mp.power(one + rp, two))

def aPhiIntegrand(kcsq, phi):
    return (mp.power(mp.sin(phi), two) - mp.power(mp.cos(phi), two)) / \
           mp.sqrt(mp.power(mp.cos(phi), two) + kcsq * mp.power(mp.sin(phi), two))

def bRhoIntegrand(kcsq, phi):
    return (mp.power(mp.sin(phi), two) - mp.power(mp.cos(phi), two)) / \
           mp.power(mp.power(mp.cos(phi), two) + kcsq * mp.power(mp.sin(phi), two), threeHalf)

def bZIntegrand(rp_, ksq_, phi_):
    return ((one - rp_) * mp.power(mp.sin(phi_), two) + (one + rp_) * mp.power(mp.cos(phi_), two)) / \
           mp.power(one - ksq_ * mp.power(mp.sin(phi_), two), threeHalf)

def A_phi(rp, zp):
    if rp == zero:
        return zero
    prefac = one / mp.sqrt(mp.power(zp, two) + mp.power(one + rp, two))
    integrand = lambda phi: aPhiIntegrand(kCSq(rp, zp), phi)
    return prefac * mp.quad(integrand, [zero, pi/two])

def B_rho(rp, zp):
    if rp == zero or zp == zero:
        return zero
    prefac = zp / mp.power(mp.power(zp, two) + mp.power(one + rp, two), threeHalf)
    integrand = lambda phi: bRhoIntegrand(kCSq(rp, zp), phi)
    return prefac * mp.quad(integrand, [zero, pi/two])

def B_z(rp_, zp_):
    if rp_ == zero:
        return pi / (two * mp.power(mp.power(zp_, two) + one, threeHalf))

    prefac_ = one / mp.power(mp.power(zp_, two) + mp.power(one + rp_, two), threeHalf)

    if zp_ == zero:
        ksq_ = four / (one/rp_ + two + rp_)
    else:
        ksq_ = kSq(rp_, zp_)

    integrand_ = lambda phi_: bZIntegrand(rp_, ksq_, phi_)
    return prefac_ * mp.quad(integrand_, [zero, pi/two])

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

    testPointsFile = "../resources/testPointsCircularWireLoop.dat"
    outFilename    = "../resources/refDataCircularWireLoop.dat"

    with open(testPointsFile, "r") as f:
        lines = f.readlines()

        with open(outFilename, "w") as outFile:

            outFile.write("# A_phi B_rho B_z of circular wire loop\n")

            numLines = len(lines)
            
            for i,line in enumerate(lines):

                # skip comment lines
                if line[0] == "#":
                    continue

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

                # compute magnetostatic quantities: A_phi, B_rho and B_z
                aPhi = A_phi(rp, zp)
                bRho = B_rho(rp, zp)
                bZ   = B_z(rp, zp)

                with mp.workdps(20):
                    outFile.write(str(aPhi) + " " + str(bRho) + " " + str(bZ) + "\n")
                    outFile.flush()

    sys.exit(0)
