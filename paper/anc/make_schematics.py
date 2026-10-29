"""Schematics and further 3D plots for 'The Energy Cost of a Phase Shift'.

Run from the paper directory: python3 anc/make_schematics.py
Writes figures/fig_slab.png, fig_contour.png, fig_hstar.png, fig_cost.png,
fig_roadmap.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import FancyBboxPatch, Arc, FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.special import ellipk, ellipe

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "stix"})


def save(fig, name, tight=False):
    from PIL import Image
    path = os.path.join(OUT, name)
    if tight:
        fig.savefig(path, dpi=200, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    else:
        fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    im = Image.open(path).convert("RGBA")
    flat = Image.new("RGB", im.size, "white")
    flat.paste(im, mask=im.split()[3])
    flat.save(path)


def Lam(r):
    kap = (r - 1.0)/(r + 1.0)
    return (ellipk(kap**2) - ellipe(kap**2))/kap


# 1. slab schematic: a pulse crossing a passive slab
fig = plt.figure(figsize=(6.4, 4.0))
ax = fig.add_subplot(111, projection="3d")
L0, L1 = 0.0, 2.0
ys, zs = (-1.0, 1.0), (-1.3, 1.3)
faces = []
for x in (L0, L1):
    faces.append([(x, ys[0], zs[0]), (x, ys[1], zs[0]), (x, ys[1], zs[1]), (x, ys[0], zs[1])])
for y in ys:
    faces.append([(L0, y, zs[0]), (L1, y, zs[0]), (L1, y, zs[1]), (L0, y, zs[1])])
for z in zs:
    faces.append([(L0, ys[0], z), (L1, ys[0], z), (L1, ys[1], z), (L0, ys[1], z)])
ax.add_collection3d(Poly3DCollection(faces, facecolor=(0.55, 0.7, 0.9, 0.18),
                                     edgecolor=(0.3, 0.4, 0.6, 0.6), lw=0.6))
t = np.linspace(-3.2, 3.2, 500)
env = np.exp(-t**2/1.6)
ax.plot(t - 4.2, np.zeros_like(t), env*np.cos(5*t), color="0.2", lw=1.1)
ax.plot(t - 4.2, np.zeros_like(t), env, color="0.2", lw=0.7, ls="--")
g = 0.72
ax.plot(t + 6.2, np.zeros_like(t), g*env*np.cos(5*t - np.pi/2), color="tab:red", lw=1.1)
ax.plot(t + 6.2, np.zeros_like(t), g*env, color="tab:red", lw=0.7, ls="--")
ax.text(-4.2, 0, 1.55, "incident pulse", ha="center", fontsize=9)
ax.text(6.2, 0, 1.55, r"delayed, carrier shifted by $\theta$," + "\n" + r"attenuated by $\geq\theta\Lambda$",
        ha="center", fontsize=9)
ax.text(1.0, 0, -1.9, r"passive medium, $n_r=n_0+C/\omega$ on $[a,b]$", ha="center", fontsize=9)
ax.text(1.0, -1.0, 1.45, r"$L$", fontsize=10)
ax.set_xlim(-7.5, 9.5); ax.set_ylim(-2, 2); ax.set_zlim(-1.6, 1.6)
ax.set_axis_off()
ax.view_init(elev=18, azim=-80)
ax.set_box_aspect((3.2, 0.9, 1.0))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
save(fig, "fig_slab.png", tight=True)

# 2. contour schematic for the kernel identity
fig, ax = plt.subplots(figsize=(6.2, 3.2))
a, b, s, Rr, eps = 1.0, 2.0, 2.7, 3.6, 0.13
ax.fill_between(np.linspace(-Rr, Rr, 300), 0, np.sqrt(np.maximum(Rr**2 - np.linspace(-Rr, Rr, 300)**2, 0)),
                color=(0.85, 0.9, 1.0))
ax.add_patch(Arc((0, 0), 2*Rr, 2*Rr, theta1=0, theta2=180, color="tab:blue", lw=1.6))
pts = [-b, -a, a, b, s]
xs = [-Rr] + sum([[p - eps, p + eps] for p in sorted(pts)], []) + [Rr]
for i in range(0, len(xs), 2):
    ax.plot([xs[i], xs[i+1]], [0, 0], color="tab:blue", lw=1.6)
for p in pts:
    ax.add_patch(Arc((p, 0), 2*eps, 2*eps, theta1=0, theta2=180, color="tab:blue", lw=1.6))
for lo, hi in ((-b, -a), (a, b)):
    ax.plot([lo + eps, hi - eps], [-0.06, -0.06], color="k", lw=3.2, solid_capstyle="butt")
ax.annotate("", xy=(0.4, Rr*0.999), xytext=(1.0, Rr*0.94),
            arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1.4))
ax.annotate("", xy=(-0.2, 0.0), xytext=(-0.8, 0.0),
            arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1.4))
for p, lab in ((-b, r"$-b$"), (-a, r"$-a$"), (a, r"$a$"), (b, r"$b$"), (s, r"$s$")):
    ax.text(p, -0.38, lab, ha="center", fontsize=10)
ax.text(1.5, 0.45, r"$W=-2+ip$", ha="center", fontsize=9)
ax.text(-1.5, 0.45, r"$W=-2-ip$", ha="center", fontsize=9)
ax.text(0.0, 1.6, r"Cauchy: $\int_{\partial D}\,\frac{W(z)}{z-s}\,dz=0$", ha="center", fontsize=11)
ax.text(3.3, 2.2, r"$|z|=R$", fontsize=9)
ax.text(0.0, 0.45, r"$W>0$", ha="center", fontsize=9)
ax.text(3.2, 0.45, r"$W>0$", ha="center", fontsize=9)
ax.set_xlim(-4.2, 4.2); ax.set_ylim(-0.7, 3.9)
ax.set_aspect("equal"); ax.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
save(fig, "fig_contour.png")

# 3. Im h_star over the upper half-plane, a=1, b=2, theta=1
a, b, theta = 1.0, 2.0, 1.0
k2 = 1 - (a/b)**2
Kk = ellipk(k2)
n = 600
phi = (np.arange(n) + 0.5)*(np.pi/2)/n; dphi = (np.pi/2)/n
tq = b*np.sqrt(1 - k2*np.sin(phi)**2)


def hstar(z):
    Rz = np.sqrt(z - a)*np.sqrt(z + a)*np.sqrt(z - b)*np.sqrt(z + b)
    integ = (dphi/(tq*(z**2 - tq**2))).sum()
    return 2*theta*z*Rz/np.pi*integ


x = np.linspace(0.02, 3.2, 121); y = np.linspace(0.04, 1.4, 70)
X, Y = np.meshgrid(x, y)
Z = np.vectorize(lambda u, v: hstar(u + 1j*v).imag)(X, Y)
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(X, Y, Z, cmap=cm.magma_r, rstride=1, cstride=1, linewidth=0, antialiased=True, alpha=0.95)
tau = 2*theta*Kk/(np.pi*b)
ax.plot_surface(X, Y, tau*Y, color="0.6", alpha=0.25, linewidth=0)
xb = np.linspace(a + 1e-3, b - 1e-3, 80)
ab = []
for xx in xb:
    m = 400
    u = (np.arange(m) + 0.5)/m
    tt = a + (b - a)*u
    rr = np.sqrt((tt**2 - a*a)*(b*b - tt**2))
    ab.append(2*theta*np.sqrt((xx**2 - a*a)*(b*b - xx**2))/np.pi*((b - a)/m*(1/(rr*(xx + tt)))).sum())
ax.plot(xb, np.zeros_like(xb), ab, color="k", lw=2.0)
ax.set_xlabel(r"$\mathrm{Re}\,z$", labelpad=4)
ax.set_ylabel(r"$\mathrm{Im}\,z$", labelpad=4)
ax.set_zlabel(r"$\mathrm{Im}\,h_\star(z)$", labelpad=4)
ax.view_init(elev=26, azim=-62)
ax.set_box_aspect((1.5, 0.9, 0.75))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_hstar.png")

# 4. the cost W(t) of a lossless resonance at t, over t and the band ratio
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
lr = np.linspace(0.08, 0.7, 60)


def Wreal(t, bb):
    rho = np.sqrt((t**2 - 1.0)/(t**2 - bb**2))
    return (rho - 1)**2/rho


u = np.linspace(0.0, 1.0, 160)
for part in (0, 1):
    T = np.empty((lr.size, u.size)); Zc = np.empty_like(T)
    for i, l in enumerate(lr):
        bb = 10**l
        if part == 0:
            tt = (1.0 - 1e-4)*(1 - (1 - u)**3)
        else:
            tt = bb + 1e-4 + (6.0 - bb)*u**3
        T[i] = tt
        Zc[i] = np.minimum(np.log10(1 + Wreal(tt, bb)), 1.2)
    LR = np.repeat(lr[:, None], u.size, axis=1)
    ax.plot_surface(T, LR, Zc, cmap=cm.viridis, rstride=1, cstride=2, linewidth=0,
                    antialiased=True, vmin=0, vmax=1.2)
bl = np.linspace(0.08, 0.7, 60)
ax.plot(10**bl, bl, np.zeros_like(bl), color="k", lw=1.6)
ax.plot(np.ones_like(bl), bl, np.zeros_like(bl), color="k", lw=1.6)
ax.set_xlabel(r"resonance frequency $t/a$", labelpad=5)
ax.set_ylabel(r"$\log_{10}(b/a)$", labelpad=5)
ax.set_zlabel(r"$\log_{10}(1+W(t))$", labelpad=5)
ax.view_init(elev=28, azim=-58)
ax.set_box_aspect((1.5, 0.9, 0.75))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_cost.png")

# 5. roadmap of the paper
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.axis("off")


def box(x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec="0.3", lw=0.8))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=8.2)


def arrow(p, q):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10, color="0.25", lw=0.9))


inp = (0.92, 0.95, 1.0); core = (1.0, 0.92, 0.82); out = (0.88, 0.96, 0.88)
box(0.2, 5.6, 2.8, 1.0, "passive systems\nHerglotz form $h$", inp)
box(3.6, 5.6, 2.8, 1.0, "test function $W$\n" + r"$\mathrm{Re}\,W=-2$ on band, $>0$ off", inp)
box(7.0, 5.6, 2.8, 1.0, "elliptic integrals\nLanden, Legendre", inp)
box(2.0, 3.6, 4.2, 1.1, "sum rule\n" + r"$\int\alpha=\frac{1}{2}\int\varphi\,p+\frac{\pi}{4}\int W\,d\mu$", core)
box(7.0, 3.6, 2.8, 1.1, "the constant\n" + r"$\Lambda=(K(\kappa)-E(\kappa))/\kappa$", core)
box(0.2, 1.4, 2.2, 1.2, "main theorem\n" + r"mean $\alpha\geq\theta\Lambda$", out)
box(2.7, 1.4, 2.2, 1.2, "lead, worst\nfrequency, energy", out)
box(5.2, 1.4, 2.2, 1.2, "media: anomalous\ndispersion, stability", out)
box(7.7, 1.4, 2.1, 1.2, "extremal $H_\\star$\none zero beats it", out)
arrow((1.6, 5.6), (3.2, 4.7)); arrow((5.0, 5.6), (4.6, 4.7)); arrow((8.4, 5.6), (8.4, 4.7))
arrow((3.0, 3.6), (1.3, 2.6)); arrow((4.0, 3.6), (3.8, 2.6)); arrow((5.2, 3.6), (6.3, 2.6))
arrow((8.4, 3.6), (8.7, 2.6)); arrow((7.4, 3.6), (1.9, 2.6))
arrow((5.6, 3.6), (8.2, 2.6))
ax.text(5.0, 0.5, "inputs (blue), the identity and the constant (orange), consequences (green)",
        ha="center", fontsize=8, color="0.3")
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
save(fig, "fig_roadmap.png")
print("schematics written")
