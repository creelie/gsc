"""Figures for the exact price on several bands (ssec:severalexact).

Run from the paper directory: python3 anc/make_figures_bands.py
Writes figures/fig_bandtest.png, figures/fig_bandloss.png and
figures/fig_bandprice.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "cm"})
NCH = 400
PH = (np.arange(NCH) + 0.5)*np.pi/NCH


def save(fig, name):
    from PIL import Image
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    im = Image.open(path).convert("RGBA")
    flat = Image.new("RGB", im.size, "white")
    flat.paste(im, mask=im.split()[3])
    flat.save(path)


class Bands:
    def __init__(self, ends):
        self.a = np.array(ends[0::2], float)
        self.b = np.array(ends[1::2], float)
        self.m = len(self.a)
        self.A, self.B = self.a**2, self.b**2
        X = lambda s: np.prod([s - u for u in self.A], axis=0)
        Y = lambda s: np.prod([s - u for u in self.B], axis=0)
        self.X, self.Y = X, Y
        self.g = np.array([brentq(lambda s: X(s) - Y(s), self.B[j], self.A[j + 1])
                           for j in range(self.m - 1)])
        self.c = np.sqrt(self.g)

    def P(self, s):
        return np.prod([s - u for u in self.g], axis=0) if self.m > 1 else np.ones_like(s)

    def quad(self, F):
        """sum_j int_{band j} F(t) dt/|r(t)| by Gauss-Chebyshev in s = t^2."""
        tot = 0.0
        for j in range(self.m):
            s = self.A[j] + (self.B[j] - self.A[j])*(1 - np.cos(PH))/2
            oth = np.ones_like(s)
            for k in range(self.m):
                if k != j:
                    oth *= np.abs((s - self.A[k])*(s - self.B[k]))
            tot += np.pi/NCH*np.sum(F(np.sqrt(s))/(2*np.sqrt(s)*np.sqrt(oth)))
        return tot

    def J(self):
        return self.quad(lambda t: np.abs(self.Y(t*t)) - np.abs(self.X(t*t)))

    def length(self):
        return float(np.sum(self.b - self.a))

    def alpha(self, x):
        """Attenuation of the optimal system for theta = 1 at x in a band."""
        rr = np.sqrt(np.abs(self.X(x*x)*self.Y(x*x)))
        rho = rr/np.abs(self.P(x*x))
        return np.array([2*rh/np.pi*self.quad(lambda t: np.abs(self.P(t*t))/(xx + t))
                         for xx, rh in zip(np.atleast_1d(x), np.atleast_1d(rho))])

    def rhoB(self, z):
        out = np.ones_like(z, dtype=complex)
        for aj, bj in zip(self.a, self.b):
            out *= np.sqrt((z - aj)/(z - bj))*np.sqrt((z + aj)/(z + bj))
        return out


def fig_bandtest():
    B = Bands((1, 2, 3, 4))
    x = np.linspace(0.02, 5.0, 340)
    y = np.linspace(0.02, 2.0, 160)
    Xg, Yg = np.meshgrid(x, y)
    rho = B.rhoB(Xg + 1j*Yg)
    W = np.real((rho - 1)**2/rho)
    W = np.clip(W, -2.2, 2.2)
    fig = plt.figure(figsize=(6.4, 4.3))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(Xg, Yg, W, cmap=cm.viridis, rstride=2, cstride=2, linewidth=0,
                    antialiased=True, alpha=0.95)
    for aj, bj in zip(B.a, B.b):
        ax.plot([aj, bj], [0, 0], [-2, -2], color="k", lw=3)
    ax.plot([B.c[0], B.c[0]], [0, 0], [0, 2.2], color="crimson", lw=1.6, zorder=10)
    ax.text(B.c[0] + 0.1, 0, 2.3, r"$c_1=\sqrt{11/2}$", color="crimson", zorder=11)
    ax.set_xlabel(r"$\mathrm{Re}\,z$")
    ax.set_ylabel(r"$\mathrm{Im}\,z$")
    ax.set_zlabel(r"$\mathrm{Re}\,W_B$")
    ax.view_init(elev=30, azim=-70)
    ax.set_title(r"$B=[1,2]\cup[3,4]$: $\mathrm{Re}\,W_B=-2$ on the bands, $W_B>0$ off them except at $c_1$",
                 fontsize=9)
    fig.tight_layout()
    save(fig, "fig_bandtest.png")


def fig_bandloss():
    fig = plt.figure(figsize=(6.4, 4.4))
    ax = fig.add_subplot(111, projection="3d")
    cols = cm.plasma(np.linspace(0.1, 0.8, 4))
    for m in range(1, 5):
        ends = []
        for j in range(m):
            ends += [2*j + 1, 2*j + 2]
        B = Bands(ends)
        for j in range(m):
            xs = np.linspace(B.a[j] + 1e-6, B.b[j] - 1e-6, 90)
            al = B.alpha(xs)
            ys = np.full_like(xs, m)
            verts = [list(zip(xs, ys, al)) + [(xs[-1], m, 0), (xs[0], m, 0)]]
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            poly = Poly3DCollection(verts, facecolor=cols[m - 1], alpha=0.75, edgecolor="k", lw=0.4)
            ax.add_collection3d(poly)
        for c in B.c:
            ax.plot([c, c], [m, m], [0, 0.25], color="crimson", lw=1.5)
        lam = B.J()/(2*B.length())
        print("m", m, "Lambda_B", lam)
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0.5, 4.5)
    ax.set_zlim(0, 0.7)
    ax.set_yticks([1, 2, 3, 4])
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"number of bands $m$")
    ax.set_zlabel(r"$\alpha_B/\theta$")
    ax.view_init(elev=24, azim=-58)
    ax.set_title(r"Least attenuation on $[1,2]\cup\dots\cup[2m-1,2m]$; resonances $c_j$ in red",
                 fontsize=9)
    fig.tight_layout()
    save(fig, "fig_bandloss.png")


def fig_bandprice():
    a2 = np.linspace(2.08, 6.0, 46)
    w = np.linspace(0.25, 3.0, 40)
    Z = np.zeros((len(w), len(a2)))
    J1 = Bands((1, 2)).J()
    for i, ww in enumerate(w):
        for k, aa in enumerate(a2):
            Jb = Bands((1, 2, aa, aa + ww)).J()
            Js = J1 + Bands((aa, aa + ww)).J()
            Z[i, k] = Jb/Js
    A2, Wg = np.meshgrid(a2, w)
    fig = plt.figure(figsize=(6.4, 4.3))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(A2, Wg, Z, cmap=cm.magma, linewidth=0, antialiased=True, alpha=0.95)
    ax.set_xlabel(r"start $a_2$ of the second band")
    ax.set_ylabel(r"width $b_2-a_2$")
    ax.set_zlabel(r"$J_B/(J_1+J_2)$")
    ax.view_init(elev=26, azim=-128)
    ax.set_title(r"Premium for one device on $[1,2]\cup[a_2,b_2]$ over two separate devices",
                 fontsize=9)
    fig.tight_layout()
    save(fig, "fig_bandprice.png")
    return Z.min(), Z.max()


if __name__ == "__main__":
    fig_bandtest()
    fig_bandloss()
    print("premium range", fig_bandprice())
    for e in [(1, 2), (1, 2, 3, 4), (1, 1.5, 4, 5), (1, 2, 3, 4, 5, 6)]:
        B = Bands(e)
        print(e, "c", B.c, "Lambda_B", B.J()/(2*B.length()))
