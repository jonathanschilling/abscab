import os
import sys
import numpy as np

abscab_path = os.path.abspath("../../../main/python")
if not abscab_path in sys.path:
    sys.path.append(abscab_path)

from abscab import cel

from abscab_util import assertRelAbsEquals

def test_cel():
    tolerance = 1.0e-15

    k_c = 0.1
    p1 =  4.1
    p2 = -4.1
    a = 1.2
    b = 1.1

    cel1 =  1.5464442694017956
    cel2 = -6.7687378198360556e-1

    c1 = cel(k_c, p1, a, b)
    c2 = cel(k_c, p2, a, b)

    # ra1 = np.abs(cel1 - c1)/(1.0 + np.abs(cel1))
    # ra2 = np.abs(cel2 - c2)/(1.0 + np.abs(cel2))
    # print("case 1: rel/abs deviation = %g"%(ra1,))
    # print("case 2: rel/abs deviation = %g"%(ra2,))

    status = 0
    status |= assertRelAbsEquals(cel1, c1, tolerance)
    status |= assertRelAbsEquals(cel2, c2, tolerance)
    return status

if __name__ == "__main__":
    status = 0
    status |= test_cel()
    sys.exit(status)
