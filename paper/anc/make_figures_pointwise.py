"""Figures for the pointwise bounds without zeros (sec:pointwise).

Run from the paper directory: python3 anc/make_figures_pointwise.py
Writes figures/fig_pointwise.png and figures/fig_leadlag.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.special import ellipk, ellipe
from scipy.integrate import quad

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
k2 = 1 - a*a/(b*b)
tau0 = 2*ellipk(k2)/(np.pi*b)
m = b*np.sqrt(ellipe(k2)/ellipk(k2))
kap = (b - a)/(b + a)
Lam = (ellipk(kap**2) - ellipe(kap**2))/kap
Lamm = (ellipe(kap**2) - (1 - kap**2)*ellipk(kap**2))/kap


def r(x):
    return np.sqrt(np.maximum((x*x - a*a)*(b*b - x*x), 0.0))


def A(x):
    if x <= a or x >= b:
        return 0.0
    def f(phi):
        t = 0.5*(a + b) - 0.5*(b - a)*np.cos(phi)
        return 1.0/((t + x)*np.sqrt((t + a)*(b + t)))
    return 2*r(x)/np.pi*quad(f, 0, np.pi, limit=200)[0]


# 1. the sharp attenuation for the phase 1 + tau*omega over the octave
xs = 0.5*(a + b) - 0.5*(b - a)*np.cos(np.linspace(0, np.pi, 91))
Ax = np.array([A(x) for x in xs])
taus = np.linspace(-2.2*tau0, 0.8*tau0, 61)
X, T = np.meshgrid(xs, taus)
C = np.maximum(0.0, -tau0 - T)
Z = np.tile(Ax, (len(taus), 1)) + C*r(X)/X
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(X, T*a, Z, cmap=cm.viridis, rstride=1, cstride=1, linewidth=0,
                antialiased=True, alpha=0.95)
ax.plot(xs, np.full_like(xs, -tau0*a), Ax + 0.004, color="crimson", lw=1.8, zorder=10)
ax.set_xlabel(r"$\omega/a$", labelpad=6)
ax.set_ylabel(r"slope $\tau a$", labelpad=6)
ax.set_zlabel(r"least $\alpha(\omega)$ (Np)", labelpad=6)
ax.view_init(elev=26, azim=-48)
ax.set_box_aspect((1.2, 1.1, 0.8))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_pointwise.png")

# 2. the sharp profiles of a unit lag and a unit lead
Lx = np.array([A(a*b/x) for x in xs])
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.plot(xs, Ax, color="C0", lw=1.6, label=r"lag $\theta=1$: $A(\omega)$")
ax.plot(xs, Lx, color="C3", lw=1.6, label=r"lead $\theta=-1$: $A(ab/\omega)$")
ax.axhline(Lam, color="C0", ls="--", lw=0.9)
ax.axhline(Lamm, color="C3", ls="--", lw=0.9)
ax.plot([m], [A(m)], "o", color="C0", ms=4)
ax.plot([a*b/m], [A(m)], "o", color="C3", ms=4)
ax.text(1.02, Lam + 0.006, r"mean $\Lambda=0.2735$", color="C0", fontsize=8.5)
ax.text(1.02, Lamm - 0.022, r"mean $\Lambda_-=0.2656$", color="C3", fontsize=8.5)
ax.text(1.5, A(m) + 0.008, r"peak $A(m)=0.3482$", fontsize=8.5, ha="center")
ax.set_xlim(a, b); ax.set_ylim(0, 0.39)
ax.set_xlabel(r"$\omega/a$"); ax.set_ylabel(r"least attenuation (Np)")
ax.legend(loc="lower center", fontsize=8.5, frameon=False)
fig.subplots_adjust(left=0.11, right=0.98, bottom=0.14, top=0.97)
save(fig, "fig_leadlag.png")
print("tau0 a", tau0*a, "A(m)", A(m), "Lam", Lam, "Lam-", Lamm)
