"""Figures for 'The Energy Cost of a Phase Shift' (matplotlib, PNG).

Run from the paper directory: python3 anc/make_figures.py
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

def Lam(r):
    kap = (r - 1.0) / (r + 1.0)
    return (ellipk(kap**2) - ellipe(kap**2)) / kap

def save(fig, name):
    from PIL import Image
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    im = Image.open(path).convert("RGBA")
    flat = Image.new("RGB", im.size, "white")
    flat.paste(im, mask=im.split()[3])
    flat.save(path)

def absR(t, a, b):
    return np.sqrt(np.abs((t*t - a*a)*(t*t - b*b)))

def alpha_ext(x, a, b, theta=1.0):
    val = quad(lambda t: 1.0/((t + x)*absR(t, a, b)), a, b, limit=200)[0]
    return 2*theta*absR(x, a, b)/np.pi*val

# 1. the bound as a surface over phase and band ratio (dB)
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
deg = np.linspace(0, 90, 46)
lr = np.linspace(np.log10(1.02), 1.0, 46)
D, LR = np.meshgrid(deg, lr)
Z = 8.685889638 * np.radians(D) * Lam(10**LR)
ax.plot_surface(D, LR, Z, cmap=cm.viridis, rstride=1, cstride=1, linewidth=0, antialiased=True, alpha=0.95)
ax.set_xlabel(r"phase shift $\theta$ (degrees)", labelpad=6)
ax.set_ylabel(r"$\log_{10}(b/a)$", labelpad=6)
ax.set_zlabel("least mean loss (dB)", labelpad=6)
ax.view_init(elev=24, azim=-128)
ax.set_box_aspect((1.25, 1.0, 0.8))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_bound.png")

# 2. extremal attenuation profiles as ribbons
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
ratios = [1.25, 1.5, 2.0, 3.0, 5.0, 10.0]
cols = cm.plasma(np.linspace(0.1, 0.85, len(ratios)))
for j, r in enumerate(ratios):
    a, b = 1.0, r
    u = np.linspace(0.0, 1.0, 120)[1:-1]
    xs = a + (b - a)*u
    al = np.array([alpha_ext(x, a, b) for x in xs])
    u = np.concatenate([[0], u, [1]]); al = np.concatenate([[0], al, [0]])
    y = np.full_like(u, np.log10(r))
    verts = list(zip(u, y, al))
    ax.plot(u, y, al, color=cols[j], lw=1.8)
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    poly = Poly3DCollection([verts + [(1, y[0], 0), (0, y[0], 0)]], facecolor=cols[j], alpha=0.35, edgecolor="none")
    ax.add_collection3d(poly)
    ax.plot([0, 1], [y[0], y[0]], [Lam(r), Lam(r)], color=cols[j], lw=0.9, ls="--")
ax.set_xlabel(r"position in band $(\omega-a)/(b-a)$", labelpad=6)
ax.set_ylabel(r"$\log_{10}(b/a)$", labelpad=6)
ax.set_zlabel(r"$\alpha_\star/\theta$ (nepers/rad)", labelpad=6)
ax.view_init(elev=22, azim=-60)
ax.set_box_aspect((1.25, 1.0, 0.8))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_extremal.png")

# 3. the test function: Re W over the upper half-plane, a=1, b=2
def W(z, a=1.0, b=2.0):
    rho = np.sqrt((z - a)/(z - b))*np.sqrt((z + a)/(z + b))
    return (rho - 1)**2/rho
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
x = np.linspace(-3.2, 3.2, 161); y = np.linspace(0.02, 1.6, 80)
X, Y = np.meshgrid(x, y)
ReW = np.real(W(X + 1j*Y))
ReW = np.clip(ReW, -2.2, 2.2)
ax.plot_surface(X, Y, ReW, cmap=cm.coolwarm, rstride=1, cstride=1, linewidth=0, antialiased=True, vmin=-2.2, vmax=2.2)
for s in (-1, 1):
    xx = np.linspace(1.0, 2.0, 50)*s
    ax.plot(xx, np.zeros_like(xx), -2*np.ones_like(xx), color="k", lw=2.2)
ax.set_xlabel(r"$\mathrm{Re}\,z$", labelpad=4)
ax.set_ylabel(r"$\mathrm{Im}\,z$", labelpad=4)
ax.set_zlabel(r"$\mathrm{Re}\,W(z)$", labelpad=4)
ax.set_yticks([0, 0.8, 1.6]); ax.set_zticks([-2, -1, 0, 1, 2])
ax.view_init(elev=28, azim=-70)
ax.set_box_aspect((1.6, 0.8, 0.7))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_testfunction.png")

# 4. one zero on the imaginary axis: mean attenuation over (theta, eta), a=1, b=2
a, b = 1.0, 2.0
k2 = 1 - (a/b)**2
n = 400
phi = (np.arange(n) + 0.5)*(np.pi/2)/n; dphi = (np.pi/2)/n
t = b*np.sqrt(1 - k2*np.sin(phi)**2)
wdt = (a*a + b*b - 2*t*t)/t*dphi
J = wdt.sum(); L = J/(2*(b - a))
th = np.linspace(0.05, 3.0, 60)
fig = plt.figure(figsize=(6.2, 4.6))
ax = fig.add_subplot(111, projection="3d")
TH, U = np.meshgrid(th, np.linspace(0.0, 1.5, 60))
Zr = np.empty_like(TH)
for i in range(TH.shape[0]):
    for j in range(TH.shape[1]):
        theta = TH[i, j]; beta = b/np.tan(theta/2)*10**U[i, j]
        red = (np.arctan(t/beta)*wdt).sum()/(b - a)
        Zr[i, j] = (theta*L - red)/(theta*L)
LB = U
ax.plot_surface(TH, LB, Zr, cmap=cm.viridis, rstride=1, cstride=1, linewidth=0, antialiased=True, vmin=0, vmax=1)
ax.plot_surface(TH, LB, np.ones_like(TH), color="0.75", alpha=0.25, linewidth=0)
ax.set_xlabel(r"phase shift $\theta$ (rad)", labelpad=6)
ax.set_ylabel(r"$\log_{10}(\eta/\eta_\theta)$", labelpad=6)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r"mean loss $/\,\theta\Lambda$", labelpad=8, rotation=90)
ax.set_zlim(0, 1.05)
ax.view_init(elev=24, azim=-135)
ax.set_box_aspect((1.25, 1.0, 0.8))
fig.subplots_adjust(left=0.02, right=0.95, bottom=0.04, top=0.98)
save(fig, "fig_zero.png")

# 5. delay against constant phase: phase as a function of frequency and the pulse picture
fig = plt.figure(figsize=(6.4, 4.4))
ax = fig.add_subplot(111, projection="3d")
tt = np.linspace(-6, 6, 600)
env = np.exp(-tt**2/6)
for k, (lab, shift, cep, col) in enumerate([("input", 0.0, 0.0, "0.25"),
                                             ("delay", 1.5, 0.0, "tab:blue"),
                                             ("constant phase", 0.0, np.pi/2, "tab:red")]):
    yv = np.full_like(tt, k)
    ax.plot(tt, yv, np.exp(-(tt - shift)**2/6)*np.cos(3*(tt - shift) - cep), color=col, lw=1.2)
    ax.plot(tt, yv, np.exp(-(tt - shift)**2/6), color=col, lw=0.8, ls="--")
ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["input", "delay", r"phase $\pi/2$"])
ax.set_xlabel("time", labelpad=4); ax.set_zlabel("field", labelpad=2)
ax.view_init(elev=26, azim=-62)
ax.set_box_aspect((1.8, 1.0, 0.6))
fig.subplots_adjust(left=0.0, right=0.96, bottom=0.04, top=0.98)
save(fig, "fig_pulses.png")
print("figures written; Lambda(2) =", Lam(2.0), " J/2 =", J/2)
