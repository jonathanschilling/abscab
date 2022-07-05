import numpy as np

def assertRelAbsEquals(expected, actual, tolerance):
    """Check if actual is equal to expected within tolerance.
    
    Check if two values are approximately equal within a prescribed tolerance.
    For values much smaller than 1, this is similar to a comparison of the
    absolute values. For values much greater than 1, this is similar to a
    comparison of the relative values.
   
    This method is described in Gill, Murray & Wright, "Practical Optimization" (1984).
   
    :param float expected:  expected result
    :param float actual:    actual result
    :param float tolerance: relative or absolute tolerance on the mismatch between the
                            expected and the actual values
    :return: 0 if the values match within the prescribed tolerance; 1 otherwise
    :rtype: int
    """
    relAbsError = np.abs(actual - expected) / (1.0 + np.abs(expected))
    if relAbsError > tolerance:
        print("expected %g, actual %g (rel/abs error %g, tolerance %g)"%
              (expected, actual, relAbsError, tolerance))
        return 1
    return 0

def errorMetric(ref, act):
    bad = 0.0
    good = -16.0

    tenToBad = 10**bad

    if np.abs(ref) > 0.0:
        if act != ref:
            relErr = np.abs((act-ref)/ref)
            return np.log10(tenToBad if tenToBad < relErr else relErr)
        else:
            return good
    else:
        return bad if np.abs(act) > 0.0 else good
