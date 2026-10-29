"""Figures for delay, channels and bands (sec:beyond).

Run from the paper directory: python3 anc/make_figures_beyond.py
Writes figures/fig_exactmean.png and figures/fig_delay.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.special import ellipk, ellipe
from scipy.integrate import quad
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "cm"})


def save(fig, name):
    from PIL import Image
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    im = Image.open(path).convert("RGBA")
    flat = Image.new("RGB", im.size, "white")
    flat.paste(im, mask=im.split()[3])
    flat.save(path)


a, b = 1.0, 2.0
f, l = 0.5*(a + b), 0.5*(b - a)
tau0 = 2*ellipk(1 - a*a/(b*b))/(np.pi*b)
kap = (b - a)/(b + a)
Lam = (ellipk(kap**2) - ellipe(kap**2))/kap
G0 = np.pi*l*l


def r(x):
    return np.sqrt(np.maximum((x*x - a*a)*(b*b - x*x), 0.0))


def tq(ph):
    return f - l*np.cos(ph)


def q(y):
    g = lambda ph: np.sqrt(y)/(np.sqrt(y) + tq(ph))/np.sqrt((tq(ph) + a)*(tq(ph) + b))
    return 2/np.pi*quad(g, 0, np.pi, limit=200)[0]


def dq(y):
    g = lambda ph: tq(ph)/(2*np.sqrt(y)*(np.sqrt(y) + tq(ph))**2)/np.sqrt((tq(ph) + a)*(tq(ph) + b))
    return 2/np.pi*quad(g, 0, np.pi, limit=200)[0]


def xv(v):
    return (v - np.sqrt((v - a*a)*(v - b*b)) + a*b)/(a + b)


def Gp(v):
    return np.pi/2*(v - (a*a + b*b)/2 - np.sqrt((v - a*a)*(v - b*b)))


def Gm(u):
    return np.pi/2*((a*a + b*b)/2 - u - np.sqrt((a*a - u)*(b*b - u)))


def optimum(s, tau, qf=None, dqf=None, beta=None):
    """Theorem 9.2: profile function y -> M(y) - s q(y) and the mean."""
    qf = qf or (lambda y: s*q(y))
    dqf = dqf or (lambda y: s*dq(y))
    beta = tau - s*tau0 if beta is None else beta
    if beta + qf(b*b) <= 0:
        return (lambda y: -beta - qf(y)), None
    eq = lambda v: dqf(xv(v)**2)*(v - xv(v)**2) - qf(xv(v)**2) - beta
    v = brentq(eq, b*b*(1 + 1e-12), 1e7, xtol=1e-14)
    ys = xv(v)**2
    D = dqf(ys)*(v - ys)**2/v
    C = qf(ys) - dqf(ys)*ys*(v - ys)/v
    return (lambda y: C + D*y/(v - y) - qf(y)), (v, ys, C, D)


xs = f - l*np.cos(np.linspace(0, np.pi, 81))
qx = np.array([q(x*x) for x in xs])

# 1. the least attenuation for the phase -1 + tau*omega
taus = np.linspace(-0.4, 3.0, 52)
tc = tau0 - q(b*b)
Z = np.zeros((len(taus), len(xs)))
for i, t in enumerate(taus):
    prof, _ = optimum(1.0, t)
    Z[i] = r(xs)/xs*np.array([prof(x*x) for x in xs])
X, T = np.meshgrid(xs, taus)
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(X, T, Z, cmap=cm.viridis, rstride=1, cstride=1, linewidth=0,
                antialiased=True, alpha=0.95)
prof, _ = optimum(1.0, tc)
ax.plot(xs, np.full_like(xs, tc), r(xs)/xs*np.array([prof(x*x) for x in xs]) + 0.004,
        color="crimson", lw=1.8, zorder=10)
cont = []
for t in taus[taus > tc]:
    _, par = optimum(1.0, t)
    cont.append((np.sqrt(par[1]), t))
cont = np.array(cont)
ax.plot(cont[:, 0], cont[:, 1], np.full(len(cont), 0.002), color="orange", lw=1.6, zorder=11)
ax.set_xlabel(r"$\omega/a$", labelpad=6)
ax.set_ylabel(r"slope $\tau a$", labelpad=6)
ax.set_zlabel(r"least $\alpha(\omega)$ (Np)", labelpad=6)
ax.view_init(elev=26, azim=-130)
ax.set_box_aspect((1.2, 1.1, 0.8))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_exactmean.png")

# 2. the argument pi/2 + 5 omega on the octave
s, tau = np.pi/2, 5.0
asharp = np.log(((a + b)*xs + r(xs))/(xs*xs + a*b))
prof0, par0 = optimum(s, tau)
alpha0 = r(xs)/xs*np.array([prof0(x*x) for x in xs])
eps = 0.1
Qu = lambda u: np.sqrt((1 + u*u)*(4 + u*u))
P = lambda y: 2*quad(lambda u: y/((y + u*u)*Qu(u)), 0, eps)[0]
dP = lambda y: 2*quad(lambda u: u*u/((y + u*u)**2*Qu(u)), 0, eps)[0]
be = tau - s*tau0 + 2*quad(lambda u: 1/Qu(u), 0, eps)[0]
prof1, par1 = optimum(s, tau, lambda y: s*q(y) - P(y), lambda y: s*dq(y) - dP(y), be)
alpha1 = r(xs)/xs*np.array([prof1(x*x) for x in xs])

fig = plt.figure(figsize=(7.0, 3.3))
ax = fig.add_subplot(1, 2, 1)
ax.plot(xs, asharp, color="0.35", lw=1.6, label=r"$e^{5iz}H^\sharp_{\pi/2}$")
ax.plot(xs, alpha0, color="C0", lw=1.6, label="zero-free optimum")
ax.plot(xs, alpha1, color="C3", lw=1.3, ls="--", label=r"one zero at $i/10$")
ax.set_xlim(a, b); ax.set_ylim(0, 0.4)
ax.set_xlabel(r"$\omega/a$"); ax.set_ylabel(r"attenuation (Np)")
ax.legend(loc="upper right", fontsize=7.5, frameon=False)
ins = inset_axes(ax, width="48%", height="38%", loc="center left", borderpad=1.4)
ins.plot(xs, alpha0, color="C0", lw=1.3)
ins.plot(xs, alpha1, color="C3", lw=1.1, ls="--")
ins.set_xlim(a, b); ins.set_ylim(0, 0.016)
ins.tick_params(labelsize=6.5)
ins.set_title(r"$\times 25$", fontsize=7, pad=2)
ax3 = fig.add_subplot(1, 2, 2, projection="3d")
y0s = np.linspace(1.08, 3.92, 30)
T0 = 16.0
Zd = np.zeros((len(y0s), len(xs)))
for i, y0 in enumerate(y0s):
    d = dq(y0)*(T0 - y0)**2/T0
    c = q(y0) - dq(y0)*y0*(T0 - y0)/T0
    Zd[i] = s*r(xs)/xs*(c + d*xs**2/(T0 - xs**2) - qx)
Xd, Yd = np.meshgrid(xs, np.sqrt(y0s))
ax3.plot_surface(Xd, Yd, Zd, cmap=cm.magma, rstride=1, cstride=1, linewidth=0,
                 antialiased=True, alpha=0.95)
ax3.plot(np.sqrt(y0s), np.sqrt(y0s), np.zeros_like(y0s), color="cyan", lw=1.6, zorder=10)
ax3.set_xlabel(r"$\omega/a$", labelpad=2, fontsize=8)
ax3.set_ylabel(r"$\sqrt{y_0}$", labelpad=2, fontsize=8)
ax3.set_zlabel(r"$\alpha_{y_0}$ (Np)", labelpad=4, fontsize=8)
ax3.tick_params(labelsize=6.5)
ax3.view_init(elev=28, azim=-60)
ax3.set_box_aspect((1.1, 1.1, 0.8), zoom=0.88)
fig.subplots_adjust(left=0.10, right=0.95, bottom=0.15, top=0.97, wspace=0.10)
save(fig, "fig_delay.png")
print("tau_c", tc, "zero-free", par0, "one zero", par1)
