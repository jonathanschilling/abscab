import numpy as np

PI_2 = np.pi / 2.0
"""half of pi"""

SQRT_EPS = np.sqrt(np.finfo(float).eps)
"""sqrt of machine epsilon"""

def cel(k_c, p, a, b):
    r"""Generalized Complete Elliptic Integral by Bulirsch (1969).
    
    The generalized complete elliptic integral is computed using an iterative algorithm
    presented by R. Bulirsch in "Numerical Calculation of Elliptic Integrals and Elliptic Functions. III"
    in "Numerische Mathematik" 13, 305-315 (1969):
    
    .. math::
    
      cel(k_c, p, a, b) = \int\limits_0^{\pi/2} \frac{a \cos^2{\varphi} + b \sin^2{\varphi}}{  \cos^2{\varphi} + p \sin^2{\varphi}} \frac{\mathrm{d}\varphi}{\sqrt{\cos^2{\varphi} + k_c^2 \sin^2{\varphi}}}
    
    :param float k_c: parameter k_c of cel(); absolute value must not be 0
    :param float p: parameter p of cel()
    :param float a: parameter a of cel()
    :param float b: parameter b of cel()
    :return: value of cel(k_c, p, a, b)
    :rtype: float
    """
    if k_c == 0.0:
        if b != 0.0:
            # when k_c is zero and b != 0, cel diverges (?)
            return np.inf
        else:
            k_c = SQRT_EPS*SQRT_EPS
    else:
        k_c = np.abs(k_c)

    m = 1.0 # \mu
    e = k_c # \nu * \mu
    # In the iterations, \nu_i is stored in k_c.

    # initialization
    if p > 0.0:
        p = np.sqrt(p)
        b = b / p
    else: # p <= 0
        f = k_c*k_c        # f = kc^2 (re-used here; later f = a_i)
        q = 1.0 - f        # 1 - kc^2
        g = 1.0 - p
        f -= p             # kc^2 - p
        q *= b-a*p         # (1 - kc^2)*(b-a*p)
        p = np.sqrt(f/g)   # sqrt((kc^2 - p)/(1-p))              --> p0 done here
        a = (a-b)/g        # (a-b)/(1-p)                         --> a0 done here
        b = -q/(g*g*p)+a*p # -(1 - kc^2)*(b-a*p)/( (1-p)^2 * p ) --> b0 done here
    
    # iteration until convergence
    while True:
        f = a    # f = a_i
        a += b/p # a_{i+1} <-- a_i + b_i/p_i
        g = e/p  # g = (\nu_i*\mu_i) / p_i
        b += f*g # b_i + a_i * (\nu_i*\mu_i) /p_i
        b += b   # b_{i+1} <-- 2*( b_i + a_i * (\nu_i*\mu_i) /p_i )
        p += g   # p_{i+1} <-- p_i + (\nu_i*\mu_i) / p_i
        g = m    # g = mu_i
        m += k_c # mu_{i+1} <-- mu_i + \nu_i
        if np.abs(g - k_c) > g * SQRT_EPS: # |\mu_i - \nu_i| > \mu_i*EPS == |1 - \nu_i / \mu_i| > EPS, but more robust!
            # not converged yet...
            k_c = np.sqrt(e) # k_c = sqrt(\nu_i * \mu_i)
            k_c += k_c       # \nu_{i+1} <-- 2*sqrt(\nu_i * \mu_i)
            e = k_c*m        # (\nu * \mu)_{i+1} <-- \nu_{i+1} * mu_{i+1} (also update product explicitly)
        else:
            break

    # final approximation:
    # \pi/2 * (a * \mu + b)/(\mu * (\mu + p))
    return PI_2 * (a*m+b) / (m*(m+p))
