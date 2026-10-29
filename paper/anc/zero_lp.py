"""Linear program for the remark on zeros with a delay (rem:onezerofamily).

For the argument theta + tau*omega on [a, b] = [1, 2] and given zeros, the
zero-free factor must carry the phase g = theta + k*pi + tau*t - (phase of
the Blaschke factors).  Its attenuation is v_g plus the band densities of the
measure off the band (atoms at u in [0, a^2) and w > b^2 in the variable t^2),
and its linear coefficient is beta_g + L - U >= 0.  The script minimises the
mean attenuation over measures carried by a fixed grid of atoms, with the
attenuation required to be nonnegative at 399 points of the band.  It is a
computation, not a proof: with no zeros it returns 0.00658 against the exact
0.0063903 of thm:exactmean.
"""
import numpy as np
from scipy.optimize import linprog

a, b = 1.0, 2.0
f, l = (a+b)/2, (b-a)/2
NQ = 400
ph = (np.arange(NQ)+0.5)*np.pi/NQ
tq = f-l*np.cos(ph)
wq = np.pi/NQ/np.sqrt((tq+a)*(tq+b))      # int_a^b F/r dt = sum wq F(tq)


def r(x):
    return np.sqrt(np.maximum((x*x-a*a)*(b*b-x*x), 0))


def vg(g, xs):
    G = g(tq)/tq
    out = []
    for x in xs:
        Gx = g(np.array([x]))[0]/x
        out.append(2*x*r(x)/np.pi*np.sum(wq*(G-Gx)*tq/(x*x-tq*tq)))
    return np.array(out)


def phase_axis(t, e):
    return np.pi-2*np.arctan(e/t)


def phase_pair(t, w):
    x0, y0 = w.real, w.imag
    return 2*np.arctan((t-x0)/y0)+2*np.arctan((t+x0)/y0)


xs = np.linspace(a, b, 401)[1:-1]


def solve(theta, tau, zeros, kshift=0):
    def g(t):
        p = np.zeros_like(t, dtype=float)
        for w in zeros:
            p = p+(phase_axis(t, w.imag) if abs(w.real) < 1e-14 else phase_pair(t, w))
        return theta+kshift*np.pi+tau*t-p
    v = vg(g, xs)
    bg = 2/np.pi*np.sum(wq*g(tq))
    y = xs**2
    us = np.linspace(0, a*a, 60, endpoint=False)
    ws = b*b*np.exp(np.linspace(np.log(1.0001), np.log(1e4), 200))
    rx = r(xs)/xs
    Gm = np.pi/2*((a*a+b*b)/2-us-np.sqrt((a*a-us)*(b*b-us)))
    Gp = np.pi/2*(ws-(a*a+b*b)/2-np.sqrt((ws-a*a)*(ws-b*b)))
    A = np.array([rx*y/(y-u) for u in us]+[rx*y/(w-y) for w in ws]).T
    Aub = np.vstack([-A, -np.concatenate([np.ones(len(us)), -np.ones(len(ws))])[None, :]])
    res = linprog(np.concatenate([Gm, Gp]), A_ub=Aub, b_ub=np.concatenate([v, [bg]]),
                  bounds=(0, None), method='highs')
    if res.status != 0:
        return np.inf
    return (np.trapezoid(v, xs)+res.fun)/(b-a)


if __name__ == "__main__":
    print("no zeros", solve(np.pi/2, 5.0, [], kshift=-1))
    print("pair at +-1.3654+1.2608i", solve(np.pi/2, 5.0, [1.3654+1.2608j], kshift=-1))
