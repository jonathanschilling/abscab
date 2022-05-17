#!/usr/bin/env python

import os
import sys
import mpmath as mp

# 300 digits of precision for calculations by default
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

# Implementations of A_phi, B_rho and B_z for reference calculations.
# For SI units and a given loop current I and a given loop radius a,
# A_phi needs to be multiplied by mu_0*I/pi and
# B_rho, B_z need to be multiplied by mu_0*I/(pi*a).

kSq = lambda rp, zp: \
    (four * rp) / \
    (mp.power(zp, two) + mp.power(one + rp, two))

kCSq = lambda rp, zp: \
    (mp.power(zp, two) + mp.power(one - rp, two)) / \
    (mp.power(zp, two) + mp.power(one + rp, two))

aPhiIntegrand = lambda kcsq, phi: \
    (mp.power(mp.sin(phi), two) - mp.power(mp.cos(phi), two)) / \
    mp.sqrt(mp.power(mp.cos(phi), two) + kcsq * mp.power(mp.sin(phi), two))

bRhoIntegrand = lambda kcsq, phi: \
    (mp.power(mp.sin(phi), two) - mp.power(mp.cos(phi), two)) / \
    mp.power(mp.power(mp.cos(phi), two) + kcsq * mp.power(mp.sin(phi), two), threeHalf)

bZIntegrand = lambda rp_, ksq_, phi_: \
    ((one - rp_) * mp.power(mp.sin(phi_), two) + (one + rp_) * mp.power(mp.cos(phi_), two)) / \
    mp.power(one - ksq_ * mp.power(mp.sin(phi_), two), threeHalf)

def aPhi(rp, zp):
    if rp == zero:
        return zero
    prefac = one / mp.sqrt(mp.power(zp, two) + mp.power(one + rp, two))
    integrand = lambda phi: aPhiIntegrand(kCSq(rp, zp), phi)
    return prefac * mp.quad(integrand, [zero, pi/two])

def bRho(rp, zp):
    if rp == zero or zp == zero:
        return zero
    prefac = zp / mp.power(mp.power(zp, two) + mp.power(one + rp, two), threeHalf)
    integrand = lambda phi: bRhoIntegrand(kCSq(rp, zp), phi)
    return prefac * mp.quad(integrand, [zero, pi/two])

def bZ(rp_, zp_):
    if rp_ == zero:
        return pi / (two * mp.power(mp.power(zp_, two) + one, threeHalf))

    prefac_ = one / mp.power(mp.power(zp_, two) + mp.power(one + rp_, two), threeHalf)

    if zp_ == zero:
        ksq_ = four / (one/rp_ + two + rp_)
    else:
        ksq_ = kSq(rp_, zp_)

    integrand_ = lambda phi_: bZIntegrand(rp_, ksq_, phi_)
    return prefac_ * mp.quad(integrand_, [zero, pi/two])

# convert a binary string representation of a 64-bit double variable into arbitrary precision
def binStr2Dbl(binStr):

    # check for correct header
    if not binStr[:2] == "0b":
        raise RuntimeError("probably not a binary string (does no start with 0b): '"+binStr+"'")

    # check for correct length
    bitStr = binStr[2:]
    if len(bitStr) != 64:
        raise RuntimeError("expected 64 bits, but got %d"%(len(bitStr),))

    # check that all bits are specified as either '0' or '1'
    bits = []
    for bit in bitStr:
        if bit == "0":
            bits.append(0)
        elif bit == "1":
            bits.append(1)
        else:
            raise RuntimeError("expected '0' or '1' but got '%s'"%(bit,))

    # parse sign bit
    # 0 means positive, 1 means negative, as in (-1)^signBit
    sign = None
    if bits[0] == 0:
        sign = one
    else:
        sign = -one

    # parse exponent
    expBits = bits[1:12]
    if len(expBits) != 11:
        raise RuntimeError("need exactly 11 exponent bits, but got %d"%(len(expBits),))
    exponentBias = mp.mpf("1023")
    exponentCharacteristic = zero
    for i,b in enumerate(expBits[::-1]):
        if b:
            exponentCharacteristic = exponentCharacteristic + mp.power(two, mp.mpf(str(i)))

    # parse mantissa
    mantissaBits = bits[12:]
    if len(mantissaBits) != 52:
        raise RuntimeError("need exactly 52 mantissa bits, but got %d"%(len(mantissaBits),))
    mantissaM = zero
    for i,b in enumerate(mantissaBits[::-1]):
        if b:
            mantissaM = mantissaM + mp.power(two, mp.mpf(str(i)))

    if exponentCharacteristic == zero and mantissaM == zero:
        # capture exact zero
        return zero

    # assemble resulting double/binary64 variable
    mantissa = one + mantissaM/mp.power(two, mp.mpf("52"))
    exponent = exponentCharacteristic - exponentBias
    return sign * mantissa * mp.power(two, exponent)

def parseArg(a):
    if a[:2] == "0b":
        # exact binary representation of double variable
        return binStr2Dbl(a)
    else:
        # decimal representation
        return mp.mpf(a)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: " + sys.argv[0] + " rhoP zP")
        print(" yields: k^2 k_c^2 A_phi B_rho B_z")
        sys.exit(1)

    # normalized radial position
    rp = parseArg(sys.argv[1])

    # normalized vertical position
    zp = parseArg(sys.argv[2])

    with mp.workdps(100):
        ksq   = kSq(rp, zp)
        kcsq  = kCSq(rp, zp)

        A_phi = aPhi(rp, zp)
        B_rho = bRho(rp, zp)
        B_z   = bZ(rp, zp)

    with mp.workdps(20):
        print(ksq, kcsq, A_phi, B_rho, B_z)

    sys.exit(0)
