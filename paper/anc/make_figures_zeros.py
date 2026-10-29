"""Figures of Section 8 (the least loss with zeros) of 'The Energy Cost of a Phase Shift'.

Run from the paper directory: python3 anc/make_figures_zeros.py
Writes figures/fig_allzeros.png, fig_lambdaall.png, fig_slitmap.png.

The optimal profile is computed by solving the boundary problem of
Lemma 8.3 directly (the band condition e^u cos(theta) = cos(v) with the
off-band data of (8.11)) by a Chebyshev collocation and a damped fixed
point, independently of the formula theta*Lambda - delta(eta_star).
"""
import os
import numpy as np
from matplotlib.patches import FancyArrowPatch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
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


def Lam(a, b):
    kap = (b - a)/(b + a)
    return (ellipk(kap**2) - ellipe(kap**2))/kap


def optimal_profile(th, a=1.0, b=2.0, N=60, iters=500):
    """Nodes x_j and the least attenuation there, for 0 < th < pi.

    v = Im zeta on the band is written as sqrt(1-s^2) sum c_n U_n(s) plus the
    tight off-band data; u = Re zeta is its Hilbert transform (Lemma 8.3)."""
    m, h = (a + b)/2, (b - a)/2
    j = np.arange(1, N + 1); y = np.cos(j*np.pi/(N + 1)); x = m + h*y
    n = np.arange(N); phi = np.arccos(y)
    U = np.sin(np.outer(phi, n + 1))/np.sin(phi)[:, None]
    T1 = np.cos(np.outer(phi, n + 1))
    M = 400; k = np.arange(1, M + 1); sq = np.cos(k*np.pi/(M + 1))
    wq = np.pi/(M + 1)*np.sin(k*np.pi/(M + 1))**2
    tq = m + h*sq
    Uq = np.sin(np.outer(np.arccos(sq), n + 1))/np.sqrt(1 - sq**2)[:, None]
    R = (h/(tq[None, :] + x[:, None])) @ (wq[:, None]*Uq)
    Tmat = (-np.pi*T1 - R)/np.pi
    u0 = 0.5*np.log((x - a)*(b - x)/((x + a)*(x + b)))   # data pi on (0,a), pi/2 on band, 0 beyond
    sw = np.sqrt(1 - y**2)
    c = np.zeros(N)
    for _ in range(iters):
        u = u0 - Tmat @ c
        d = np.arcsin(np.clip(np.cos(th)*np.exp(u), -1, 1))   # pi/2 - v
        c = 0.5*c + 0.5*np.linalg.solve(U*sw[:, None], d)
    u = u0 - Tmat @ c
    e = np.exp(u)
    alpha = 2*np.arctanh(np.sin(th)*e/np.sqrt(1 - (np.cos(th)*e)**2))
    return x, alpha


def mean_profile(x, alpha, a=1.0, b=2.0):
    xs = np.concatenate(([a], x[::-1], [b])); al = np.concatenate(([0.0], alpha[::-1], [0.0]))
    return np.trapezoid(al, xs)/(b - a) if hasattr(np, "trapezoid") else np.trapz(al, xs)/(b - a)


def alpha_right(w, a=1.0, b=2.0):
    q = np.clip((w - a)*(b - w)/((w + a)*(w + b)), 0, None)
    return 2*np.arctanh(np.sqrt(q))


def r_(t, a, b):
    return np.sqrt((t*t - a*a)*(b*b - t*t))


def eta_star(th, a=1.0, b=2.0):
    I0 = quad(lambda t: 1/r_(t, a, b), a, b, limit=200)[0]
    f = lambda e: quad(lambda t: 2*np.arctan(t/e)/r_(t, a, b), a, b, limit=200)[0] - th*I0
    return brentq(f, 1e-6, 1e6)


def mean_onezero(th, eta, a=1.0, b=2.0):
    p = lambda t: (a*a + b*b - 2*t*t)/r_(t, a, b)
    delta = quad(lambda t: np.arctan(t/eta)*p(t), a, b, limit=200)[0]/(b - a)
    return th*Lam(a, b) - delta


a, b = 1.0, 2.0

# 1. the least attenuation over (theta, omega): a 3D surface
ths = np.linspace(0.04, np.pi - 0.04, 64)
X, A = [], []
for th in ths:
    x, al = optimal_profile(th)
    X.append(np.concatenate(([a], x[::-1], [b]))); A.append(np.concatenate(([0.0], al[::-1], [0.0])))
X = np.array(X); A = np.array(A); TH = np.repeat(ths[:, None], X.shape[1], axis=1)
fig = plt.figure(figsize=(6.2, 4.7))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(TH, X, A, cmap=cm.viridis, rstride=1, cstride=1, linewidth=0, antialiased=True, alpha=0.93)
ax.set_xlabel(r"phase $\theta$", labelpad=6)
ax.set_ylabel(r"frequency $\omega/a$", labelpad=6)
ax.set_zlabel(r"least attenuation $\alpha^\sharp_\theta$ (Np)", labelpad=6)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
ax.view_init(elev=26, azim=-132)
ax.set_box_aspect((1.3, 1.0, 0.75))
fig.subplots_adjust(left=0.0, right=0.94, bottom=0.04, top=0.99)
save(fig, "fig_allzeros.png")

# 2. mean attenuation against theta
th2 = np.linspace(0.02, np.pi - 0.02, 120)
lamb = Lam(a, b)
kap = (b - a)/(b + a)
lead = (ellipe(kap**2) - (1 - kap**2)*ellipk(kap**2))/kap
opt = np.array([mean_onezero(t, eta_star(t)) for t in th2])
onez = np.array([mean_onezero(t, b/np.tan(t/2)) for t in th2])
fig, ax = plt.subplots(figsize=(5.0, 3.4))
ax.plot(th2, th2*lamb, color="tab:blue", lw=1.4, label=r"zero-free bound $\theta\Lambda$")
ax.plot(th2, (np.pi - th2)*lead, color="tab:orange", lw=1.2, ls="--", label="lead branch")
ax.plot(th2, onez, color="tab:green", lw=1.2, ls="-.", label=r"one zero at $i\eta_\theta$")
ax.plot(th2, opt, color="k", lw=1.8, label=r"least over all, $\Lambda_{\rm all}(\theta)$")
ax.plot([np.pi/2], [np.pi/2*(3 - 2*np.sqrt(2))], "o", color="tab:red", ms=5)
ax.annotate(r"$\frac{\pi}{2}(3-2\sqrt{2})$", (np.pi/2, np.pi/2*(3 - 2*np.sqrt(2))), (1.75, 0.08),
            arrowprops=dict(arrowstyle="->", lw=0.7), fontsize=9)
ax.set_xlim(0, np.pi); ax.set_ylim(0, 0.75)
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels(["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
ax.set_xlabel(r"phase $\theta$"); ax.set_ylabel("mean attenuation (Np)")
ax.legend(fontsize=8, loc="upper left", frameon=False)
ax.grid(alpha=0.25)
fig.tight_layout()
save(fig, "fig_lambdaall.png")

# 3. schematic of the conformal map of Proposition 8.6
th = 1.0; T = 1.35
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 3.0), gridspec_kw={"width_ratios": [1.05, 1]})
cols = {"circ": "tab:purple", "neg": "tab:red", "slit": "tab:green", "pos": "tab:blue"}
ax1.fill_between([-2.2, 4.6], 0, 2.2, color="0.93")
ax1.plot([-2.2, 0], [0, 0], color=cols["circ"], lw=3)
ax1.plot([0, 1], [0, 0], color=cols["neg"], lw=3)
ax1.plot([1, 4], [0, 0], color=cols["slit"], lw=3)
ax1.plot([4, 4.6], [0, 0], color=cols["pos"], lw=3)
for xv, lab in [(0, r"$0$"), (1, r"$a^2$"), (4, r"$b^2$")]:
    ax1.plot([xv], [0], "ko", ms=3.5); ax1.text(xv, -0.32, lab, ha="center", fontsize=9)
ax1.plot([-1.3], [0], "o", color="k", ms=3.5, mfc="white")
ax1.text(-1.3, -0.32, r"$-\eta_\star^2$", ha="center", fontsize=9)
ax1.text(1.2, 1.2, r"$\sigma=z^2\in\mathbb{C}_+$", fontsize=10)
ax1.text(-1.1, 0.18, r"$i\mathbb{R}_+$", fontsize=8, color=cols["circ"])
ax1.text(2.5, 0.18, "band", fontsize=8, color=cols["slit"], ha="center")
ax1.set_xlim(-2.2, 4.6); ax1.set_ylim(-0.6, 2.2); ax1.axis("off")
# right panel: G-plane
ax2.fill_between([-2.6, 2.6], -3.2, 0, color="0.93")
tt = np.linspace(np.pi, 2*np.pi, 200)
ax2.fill(np.cos(tt), np.sin(tt), color="white")
ax2.plot(np.cos(tt), np.sin(tt), color=cols["circ"], lw=3)
ax2.plot([-2.6, -1], [0, 0], color=cols["neg"], lw=3)
ax2.plot([1, 2.6], [0, 0], color=cols["pos"], lw=3)
ax2.plot([np.cos(th)]*2, [-T, -3.2], color=cols["slit"], lw=3)
ax2.plot([np.cos(th)], [-T], "ko", ms=3)
ax2.text(np.cos(th) + 0.1, -T + 0.05, r"$\cos\theta-iT$", fontsize=8)
ax2.plot([np.cos(th)], [-np.sin(th)], "o", color="k", ms=3.5, mfc="white")
ax2.text(np.cos(th) + 0.14, -np.sin(th) - 0.12, r"$e^{-i\theta}$", fontsize=8)
for xv, lab in [(-1, r"$-1$"), (1, r"$1$")]:
    ax2.plot([xv], [0], "ko", ms=3.5); ax2.text(xv, 0.14, lab, ha="center", fontsize=9)
ax2.text(-2.4, -2.6, r"$D_T$", fontsize=11)
ax2.text(np.cos(th) - 0.65, -2.9, r"$P_L$", fontsize=9); ax2.text(np.cos(th) + 0.15, -2.9, r"$P_R$", fontsize=9)
ax2.text(-0.2, 0.55, r"$G=1/F$", fontsize=10)
ax2.set_xlim(-2.6, 2.6); ax2.set_ylim(-3.2, 0.9); ax2.set_aspect("equal"); ax2.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.03, top=0.97, wspace=0.14)
fig.patches.append(FancyArrowPatch((0.475, 0.62), (0.535, 0.62), transform=fig.transFigure,
                                   arrowstyle="->", mutation_scale=12, lw=1.2, color="k"))
fig.text(0.505, 0.66, r"$\Phi$", fontsize=11, ha="center")
save(fig, "fig_slitmap.png")

x, al = optimal_profile(np.pi/2)
print("max |alpha - alpha_right| at pi/2:", np.max(np.abs(al - alpha_right(x))))
print("figures written; octave Lambda_all(pi/2) =", mean_profile(x, al), "closed form", np.pi/2*(3 - 2*np.sqrt(2)))
