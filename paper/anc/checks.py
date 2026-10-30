"""Numerical cross-checks for the paper 'The Exact Price of a Phase Shift'.

Every statement of the paper is proved analytically; this script only
re-evaluates the identities and constants at sample bands.  mpmath's
ellipk(m), ellipe(m) take the parameter m = k^2.
"""
import mpmath as mp
mp.mp.dps = 30
passed = failed = 0
def check(name, ok):
    global passed, failed
    if ok: passed += 1
    else: failed += 1
    print(("ok   " if ok else "FAIL ") + name)

def K(k): return mp.ellipk(k**2)
def E(k): return mp.ellipe(k**2)

def band(a, b):
    R = lambda t: mp.sqrt(abs((t*t-a*a)*(t*t-b*b)))
    w = lambda t: (a*a+b*b-2*t*t)/R(t)
    return R, w

def rho(z, a, b):
    return mp.sqrt((z-a)/(z-b))*mp.sqrt((z+a)/(z+b))
def W(z, a, b):
    r = rho(z, a, b); return (r-1)**2/r

for (a, b) in [(1, 2), (1, mp.mpf('1.1')), (mp.mpf('0.3'), 5)]:
    a, b = mp.mpf(a), mp.mpf(b)
    R, w = band(a, b)
    tag = " [a=%s, b=%s]" % (mp.nstr(a, 3), mp.nstr(b, 3))
    # boundary values of W on the band and off it
    t = (a+b)/2; z = t + 1j*mp.mpf('1e-25')
    check("Re W = -2, Im W = w on the band"+tag,
          abs(W(z, a, b) - (-2 + 1j*w(t))) < 1e-10)
    for s in [a/3, 3*b]:
        val = W(s + 1j*mp.mpf('1e-25'), a, b)
        check("W real and positive off the band"+tag, abs(mp.im(val)) < 1e-10 and mp.re(val) > 0)
    # kernel identity K(s) = Re W(s)
    for s in [a/2, a+(b-a)/3, (a+3*b)/4, 2*b]:
        f = lambda x: w(x)*2*x/(x*x-s*s)
        if a < s < b:
            g = lambda x: (w(x)*2*x - w(s)*2*s)/(x*x-s*s)
            # PV of 1/(x^2-s^2) over (a,b) in closed form
            pv = (mp.log(abs((b-s)/(b+s))) - mp.log(abs((a-s)/(a+s))))/(2*s)
            val = (mp.quad(g, [a, s, b]) + w(s)*2*s*pv)/mp.pi
            target = -2
        else:
            val = mp.quad(f, [a, b])/mp.pi
            target = mp.re(W(s + 1j*mp.mpf('1e-25'), a, b))
        check("kernel identity at s=%s" % mp.nstr(s, 4) + tag, abs(val-target) < 1e-8)
    # first moment of w vanishes, J closed forms, Landen
    check("int t w = 0"+tag, abs(mp.quad(lambda x: x*w(x), [a, b])) < 1e-12)
    J = mp.quad(w, [a, b])
    k = mp.sqrt(1-(a/b)**2); kp = a/b; kap = (b-a)/(b+a)
    check("J = b[(1+k'^2)K-2E]"+tag, abs(J - b*((1+kp**2)*K(k)-2*E(k))) < 1e-12)
    check("J = 2(a+b)(K-E)(kappa)"+tag, abs(J - 2*(a+b)*(K(kap)-E(kap))) < 1e-12)
    check("Landen E(k) = (1+k')E(kappa) - k'K(k)"+tag, abs(E(k)-((1+kp)*E(kap)-kp*K(k))) < 1e-20)
    Lam = J/(2*(b-a))
    check("pi kappa/4 < Lambda <= pi kappa/(4(1-kappa^2))"+tag,
          mp.pi*kap/4 < Lam <= mp.pi*kap/(4*(1-kap**2)))
    I0 = mp.quad(lambda x: 1/R(x), [a, b])
    check("I0 = K(k)/b"+tag, abs(I0 - K(k)/b) < 1e-12)
    check("int t^2/|R| = b E(k)"+tag, abs(mp.quad(lambda x: x*x/R(x), [a, b]) - b*E(k)) < 1e-12)
    check("int w/t = pi(b-a)^2/(2ab)"+tag, abs(mp.quad(lambda x: w(x)/x, [a, b]) - mp.pi*(b-a)**2/(2*a*b)) < 1e-12)
    check("M3 = -pi(b^2-a^2)^2/8"+tag, abs(mp.quad(lambda x: x**3*w(x), [a, b]) + mp.pi*(b*b-a*a)**2/8) < 1e-10)
    # extremal system for theta = 1
    theta = mp.mpf(1)
    m2 = b*b*E(k)/K(k)
    check("m^2 = b^2 E/K balances the band"+tag, abs(mp.quad(lambda s: (s*s-m2)/R(s), [a, b])) < 1e-12)
    D = mp.quad(lambda s: (m2-s*s)/R(s), [0, a])
    check("Legendre: int_0^a (m^2-s^2)/|R| = pi/(2 I0)"+tag, abs(D - mp.pi/(2*I0)) < 1e-12)
    tau = theta/D
    check("tau = 2 theta K(k)/(pi b) = 4 theta K(kappa)/(pi(a+b))"+tag,
          abs(tau - 2*theta*K(k)/(mp.pi*b)) < 1e-12 and abs(tau - 4*theta*K(kap)/(mp.pi*(a+b))) < 1e-12)
    alpha1 = lambda x: tau*mp.quad(lambda s: (m2-s*s)/R(s), [a, x])
    def alpha2(x):
        # Cauchy form 2 theta x |R(x)|/pi PV int dt/(|R(t)|(x^2-t^2))
        pv = -(mp.log(abs((b-x)/(b+x))) - mp.log(abs((a-x)/(a+x))))/(2*x)
        return 2*theta*x*R(x)/mp.pi*(mp.quad(lambda t: 1/R(t)/(x*x-t*t) - 1/R(x)/(x*x-t*t), [a, x, b]) + pv/R(x))
    alpha3 = lambda x: 2*theta*R(x)/mp.pi*mp.quad(lambda t: 1/((t+x)*R(t)), [a, b])
    for x in [a+(b-a)/5, (a+b)/2, b-(b-a)/7]:
        v1, v2, v3 = alpha1(x), alpha2(x), alpha3(x)
        check("extremal attenuation: three formulas agree, positive at x=%s" % mp.nstr(x, 4)+tag,
              abs(v1-v2) < 1e-8 and abs(v1-v3) < 1e-8 and v1 > 0)
    tot = mp.quad(alpha3, [a, b])
    check("extremal attains theta J/2"+tag, abs(tot - theta*J/2) < 1e-8)
    # phase of the extremal on the band equals theta: h(x) = tau int_0^x (s^2-m^2)/R(s) ds, R<0 on (0,a)
    check("extremal phase = theta"+tag, abs(tau*mp.quad(lambda s: (m2-s*s)/R(s), [0, a]) - theta) < 1e-12)
    # plain mean of the inverted (phase-lead) extremal: theta (E - kappa'^2 K)/kappa per unit bandwidth
    lead = mp.quad(lambda om: alpha3(a*b/om), [a, b])
    check("lead extremal mean = theta(E-kappa'^2 K)/kappa"+tag,
          abs(lead/(b-a) - theta*(E(kap)-(1-kap**2)*K(kap))/kap) < 1e-8)
    # monotone-weight lemma: int t w f(t^2) > 0 for decreasing f
    for beta in [mp.mpf('0.5'), 3, 30]:
        val = mp.quad(lambda x: mp.atan(x/beta)*w(x), [a, b])
        check("int atan(t/beta) w > 0 (beta=%s)" % mp.nstr(beta, 3)+tag, val > 0)

# table of Lambda
for r, ref in [(mp.mpf('1.01'), '0.0039075'), (mp.mpf('1.1'), '0.037432'), (mp.mpf('1.5'), '0.15950'),
               (2, '0.27353'), (3, '0.43658'), (10, '0.94478'), (100, '2.01598'), (1000, '3.15018')]:
    kap = (r-1)/(r+1); Lam = (K(kap)-E(kap))/kap
    check("Lambda(b/a=%s) = %s" % (mp.nstr(r, 4), ref), abs(Lam - mp.mpf(ref)) < 6e-6*max(1, Lam))
    print("     Lambda =", mp.nstr(Lam, 8), " wide-band approx", mp.nstr((mp.log(4*r)-2)/2, 6))

# single-zero reduction asymptotics, a = 1, b = 2
a, b = mp.mpf(1), mp.mpf(2)
R, w = band(a, b)
for beta, ref in [(3, '0.03263'), (10, '0.001144'), (30, '4.349e-5')]:
    red = mp.quad(lambda x: mp.atan(x/beta)*w(x), [a, b])/(b-a)
    check("reduction at beta=%s is %s" % (beta, ref), abs(red/mp.mpf(ref) - 1) < 2e-3)
    print("     reduction", mp.nstr(red, 6), " asymptotic", mp.nstr(mp.pi*(b*b-a*a)**2/(24*beta**3), 6))

# numbers of the source example
h = mp.mpf('6.62607015e-34'); f = mp.mpf(5e6)
check("1 degree at 5 MHz is 555.6 ps", abs(1/f/360 - mp.mpf('555.5556e-12')) < 1e-15)
check("h f dt has value 1.83e-36 J s", abs(h*(f-f/360)*mp.mpf('555e-12') - mp.mpf('1.8337e-36')) < 1e-39)
check("h f/360 = 9.20e-30 J", abs(h*f/360 - mp.mpf('9.2029e-30')) < 1e-33)
check("pi^2/720 = 0.013708", abs(mp.pi**2/720 - mp.mpf('0.0137078')) < 1e-6)

# octave extremal numbers (rem:extremal), theta = 1
a, b = mp.mpf(1), mp.mpf(2)
R, w = band(a, b)
k = mp.sqrt(1-(a/b)**2); kap = (b-a)/(b+a)
I0 = K(k)/b; m2 = b*b*E(k)/K(k); m = mp.sqrt(m2)
tau = 2*K(k)/(mp.pi*b)
Lam = (K(kap)-E(kap))/kap
def astar(x):
    # alpha* = tau int_a^x (m^2 - s^2)/|R(s)| ds, with s = a + (x - a) sin^2(u) to remove the edge singularity
    if x <= a: return mp.mpf(0)
    def g(u):
        s = a + (x-a)*mp.sin(u)**2
        d = mp.sqrt((s+a)*(b*b-s*s))
        return (m2-s*s)*2*mp.sqrt(x-a)*mp.cos(u)/d if d > 0 else mp.mpf(0)
    return tau*mp.quad(g, [0, mp.pi/2])
check("octave: tau* a = 0.6864", abs(tau - mp.mpf('0.6864')) < 5e-5)
check("octave: m = 1.4988 a", abs(m - mp.mpf('1.4988')) < 5e-5)
check("octave: alpha*(m) = 0.3482 Np", abs(astar(m) - mp.mpf('0.3482')) < 5e-5)
check("octave: mean alpha* = Lambda = 0.27353",
      abs(mp.quad(lambda u: astar(a+(b-a)*mp.sin(u)**2)*2*(b-a)*mp.sin(u)*mp.cos(u), [0, mp.pi/4, mp.pi/2])/(b-a) - Lam) < 1e-8 and abs(Lam - mp.mpf('0.27353')) < 5e-6)
check("octave: max alpha*/(theta Lambda) = 1.273", abs(astar(m)/Lam - mp.mpf('1.273')) < 5e-4)
check("W(0) = a/b + b/a - 2", abs(W(mp.mpf('1e-30')+1j*mp.mpf('1e-30'), a, b) - (a/b+b/a-2)) < 1e-12)
# lead extremal constant of tab:zeros
lead = (E(kap)-(1-kap**2)*K(kap))/kap
check("octave lead constant (E - kappa'^2 K)/kappa = 0.26560", abs(lead - mp.mpf('0.26560')) < 5e-6)
# tab:zeros: one zero at eta_theta = b cot(theta/2)
for th, c1, c2, c3 in [('0.5', '0.13677', '0.70160', '0.13443'), ('1', '0.27353', '0.56880', '0.25394'),
                       ('pi/2', '0.42966', '0.41720', '0.34659'), ('2', '0.54707', '0.30320', '0.35977'),
                       ('2.5', '0.68383', '0.17040', '0.27677')]:
    t = mp.pi/2 if th == 'pi/2' else mp.mpf(th)
    eta = b*mp.cot(t/2)
    delta = mp.quad(lambda x: mp.atan(x/eta)*w(x), [a, b])/(b-a)
    ok = (abs(t*Lam - mp.mpf(c1)) < 6e-6 and abs((mp.pi-t)*mp.mpf('0.26560') - mp.mpf(c2)) < 6e-5
          and abs(t*Lam - delta - mp.mpf(c3)) < 6e-6)
    check("tab:zeros row theta=%s" % th, ok)
    if th in ('0.5', '1'):
        check("one zero lowers the mean below theta Lambda and theta < a pi/(a+b) [theta=%s]" % th,
              t < a*mp.pi/(a+b) and delta > 0)
# cor:stable ratios J/P
for (a1, b1, ref) in [(1, 2, '0.2584'), (1, 10, '0.617')]:
    a1, b1 = mp.mpf(a1), mp.mpf(b1)
    R1, w1 = band(a1, b1)
    c = mp.sqrt((a1*a1+b1*b1)/2)
    J1 = mp.quad(w1, [a1, b1]); P1 = mp.quad(w1, [a1, c]) - mp.quad(w1, [c, b1])
    check("J/P(b/a=%s) = %s" % (mp.nstr(b1/a1, 3), ref), abs(J1/P1 - mp.mpf(ref)) < 6*10**(-len(ref)+1))
# 5 MHz band 4.5-5.5 MHz, one degree
kap = mp.mpf('0.1'); Lam5 = (K(kap)-E(kap))/kap; thl = mp.pi/180*Lam5
check("Lambda(kappa=0.1) = 0.078836", abs(Lam5 - mp.mpf('0.078836')) < 5e-7)
check("theta Lambda = 1.3760e-3 Np = 0.011951 dB",
      abs(thl - mp.mpf('1.3760e-3')) < 5e-8 and abs(thl*20/mp.log(10) - mp.mpf('0.011951')) < 5e-7)
loss = 1 - mp.e**(-2*thl)
check("worst-frequency loss 2.748e-3, formula pi^2/720 * 0.2 = 2.742e-3",
      abs(loss - mp.mpf('2.748e-3')) < 5e-7 and abs(mp.pi**2/720*mp.mpf('0.2') - mp.mpf('2.742e-3')) < 5e-7)
check("hf = 3.313e-27 J, deposit >= 9.10e-30 J", abs(h*f - mp.mpf('3.313e-27')) < 5e-31
      and abs(loss*h*f - mp.mpf('9.10e-30')) < 5e-33)
check("ratio to hf/360 is 0.989, against (pi^2/2)(Delta f/f)",
      abs(loss*360 - mp.mpf('0.989')) < 5e-4 and abs(mp.pi**2/2*mp.mpf('0.2') - 0.987) < 2e-3)
# sec:allzeros, the least loss with zeros allowed
mp.mp.dps = 20
def zsetup(a, b):
    r = lambda t: mp.sqrt((t*t-a*a)*(b*b-t*t))
    I0 = mp.quad(lambda t: 1/r(t), [a, b])
    avg = lambda fn: mp.quad(lambda t: fn(t)/r(t), [a, b])/I0
    kp = (b-a)/(b+a)
    p = lambda t: (a*a+b*b-2*t*t)/r(t)
    return r, avg, (K(kp)-E(kp))/kp, p
def etastar(th, a, b):
    r, avg, Lz, p = zsetup(a, b)
    lo, hi = mp.mpf('1e-8')*a, mp.mpf('1e8')*b
    for _ in range(160):
        mid = mp.sqrt(lo*hi)
        if 2*avg(lambda t: mp.atan(t/mid)) > th: lo = mid
        else: hi = mid
    return mp.sqrt(lo*hi)
def Lall(th, a, b):
    r, avg, Lz, p = zsetup(a, b); e = etastar(th, a, b)
    return e, th*Lz - mp.quad(lambda t: mp.atan(t/e)*p(t), [a, b])/(b-a), Lz
def vg(x, th, eta, a, b):
    r = lambda t: mp.sqrt((t*t-a*a)*(b*b-t*t))
    G = lambda t: (th-2*mp.atan(t/eta))/t
    if x <= a or x >= b: return mp.mpf(0)
    def fn(t):
        d = r(t)*(x*x-t*t)
        return (G(t)-G(x))/d*t if d != 0 else mp.mpf(0)
    return 2*x*r(x)/mp.pi*mp.quad(fn, [a, x, b])
def aright(x, a, b):
    return mp.log(((a+b)*x + mp.sqrt((x*x-a*a)*(b*b-x*x)))/(x*x+a*b))
for rb in [2, 3, 10]:
    a, b = mp.mpf(1), mp.mpf(rb)
    tag = " [b/a=%s]" % rb
    r, avg, Lz, p = zsetup(a, b)
    check("eta* = sqrt(ab) at a right angle"+tag, abs(2*avg(lambda t: mp.atan(t/mp.sqrt(a*b))) - mp.pi/2) < 1e-12)
    xs = [a+(b-a)*mp.mpf(j)/7 for j in (1, 3, 6)]
    check("v_g at pi/2 equals both closed forms of alpha#"+tag,
          all(abs(vg(x, mp.pi/2, mp.sqrt(a*b), a, b) - aright(x, a, b)) < 1e-9 and
              abs(aright(x, a, b) - 2*mp.atanh(mp.sqrt((x-a)*(b-x)/((x+a)*(x+b))))) < 1e-15 for x in xs))
    tot = mp.quad(lambda u: aright(a+(b-a)*mp.sin(u)**2, a, b)*2*(b-a)*mp.sin(u)*mp.cos(u), [0, mp.pi/2])
    check("int alpha# = (pi/2)(sqrt b - sqrt a)^2"+tag, abs(tot - mp.pi/2*(mp.sqrt(b)-mp.sqrt(a))**2) < 1e-12)
    check("Lambda_all(pi/2) = theta Lambda - delta(sqrt(ab))"+tag,
          abs(Lall(mp.pi/2, a, b)[1] - mp.pi/2*(mp.sqrt(b)-mp.sqrt(a))/(mp.sqrt(b)+mp.sqrt(a))) < 1e-10)
    w0 = mp.sqrt(a*b)
    check("alpha#(sqrt(ab)) = log(b/a)/2, symmetric under w -> ab/w"+tag,
          abs(aright(w0, a, b) - mp.log(b/a)/2) < 1e-15 and abs(aright(xs[0], a, b) - aright(a*b/xs[0], a, b)) < 1e-15
          and all(aright(x, a, b) < aright(w0, a, b) for x in xs if x != w0))
check("octave mean (pi/2)(3-2 sqrt2) = 0.26951, 62.7% of 0.42966",
      abs(mp.pi/2*(3-2*mp.sqrt(2)) - mp.mpf('0.26951')) < 5e-6 and abs((3-2*mp.sqrt(2))/mp.mpf('0.27353')/1 - mp.mpf('0.6272')) < 5e-4)
check("10 log10(2) = 3.0103 dB", abs(10*mp.log10(2) - mp.mpf('3.0103')) < 5e-5)
# tab:allzeros
a, b = mp.mpf(1), mp.mpf(2)
for th, c1, c2, c3, c4, c0 in [('0.5', '0.72576', '0.13090', '0.9571', '0.16677', '0.13443'),
                               ('1', '0.71853', '0.22860', '0.8357', '0.29204', '0.25394'),
                               ('pi/2', '0.70711', '0.26951', '0.6272', '0.34657', '0.34659'),
                               ('2', '0.69843', '0.24357', '0.4452', '0.31540', '0.35977'),
                               ('2.5', '0.69050', '0.15941', '0.2331', '0.20805', '0.27677')]:
    t = mp.pi/2 if th == 'pi/2' else mp.mpf(th)
    e, L, Lz = Lall(t, a, b)
    xm = mp.findroot(lambda x: mp.diff(lambda y: vg(y, t, e, a, b), x), mp.mpf('1.42'))
    ok = (abs(e/(b/mp.tan(t/2)) - mp.mpf(c1)) < 6e-6 and abs(L - mp.mpf(c2)) < 6e-6
          and abs(L/(t*Lz) - mp.mpf(c3)) < 6e-5 and abs(vg(xm, t, e, a, b) - mp.mpf(c4)) < 6e-6
          and e < b/mp.tan(t/2) and L < mp.mpf(c0))
    check("tab:allzeros row theta=%s (and below the one-zero value %s)" % (th, c0), ok)
check("worst frequency at theta=1: 0.29204 = 1.068 theta Lambda", abs(mp.mpf('0.29204')/mp.mpf('0.27353') - mp.mpf('1.068')) < 5e-4)
# the optimal profile is a mean: mean of v_g = theta Lambda - delta(eta*), and v_g > 0 on the band
for rb in [2, 10, 30]:
    a, b = mp.mpf(1), mp.mpf(rb)
    for th in [mp.mpf('0.3'), mp.mpf(1), mp.mpf('2.9')]:
        e = etastar(th, a, b)
        xs = [a+(b-a)*mp.mpf(j)/12 for j in range(1, 12)]
        check("alpha# > 0 on the band [b/a=%s, theta=%s]" % (rb, mp.nstr(th, 2)),
              all(vg(x, th, e, a, b) > 0 for x in xs))
a, b = mp.mpf(1), mp.mpf(2)
e, L, Lz = Lall(mp.mpf(1), a, b)
mean_vg = mp.quad(lambda ph: vg((a+b)/2-(b-a)/2*mp.cos(ph), 1, e, a, b)*(b-a)/2*mp.sin(ph), [0, mp.pi/2, mp.pi])/(b-a)
check("octave theta=1: mean of alpha# = theta Lambda - delta(eta*)", abs(mean_vg - L) < 1e-8)
check("octave theta=1: zero-free extremal alpha* >= alpha# pointwise",
      all(astar(x) > vg(x, 1, e, a, b) for x in [a+(b-a)*mp.mpf(j)/10 for j in range(1, 10)]))
# cor:narrowzeros
for rb, ref in [(mp.mpf('1.2'), '1.0021'), (mp.mpf('1.05'), '1.000149')]:
    a, b = mp.mpf(1), rb
    eps = 2*(b-a)/(a+b)
    val = mp.pi/2*(mp.sqrt(b)-1)/(mp.sqrt(b)+1)/(mp.pi/8*eps)
    check("narrow band ratio %s at b/a=%s" % (ref, mp.nstr(rb, 3)), abs(val - mp.mpf(ref)) < 6*10**(-len(ref)+2))
a, b = mp.mpf('4.5'), mp.mpf('5.5')
e, L, Lz = Lall(mp.pi/180, a, b)
check("5 MHz, 1 degree: Lambda_all = 1.37588e-3 < theta Lambda = 1.37595e-3",
      abs(L - mp.mpf('1.37588e-3')) < 5e-9 and abs(mp.pi/180*Lz - mp.mpf('1.37595e-3')) < 5e-9)
# cor:nome
for rb, ref in [(2, '0.085796'), (10, '0.26220')]:
    bb = mp.mpf(rb)
    q = mp.exp(-mp.pi*K(1/bb)/K(mp.sqrt(1-1/bb**2)))
    check("nome q(b/a=%s) = %s" % (rb, ref), abs(q - mp.mpf(ref)) < 6e-6)
for rb in [2, 10]:
    a, b = mp.mpf(1), mp.mpf(rb)
    I0 = K(mp.sqrt(1-(a/b)**2))/b
    xph = lambda ph: (a+b)/2-(b-a)/2*mp.cos(ph)
    def ltanh(ph):
        X = xph(ph)
        return (mp.log((b-a)*mp.sin(ph/2)**2) + mp.log((b-a)*mp.cos(ph/2)**2) - mp.log((X+a)*(X+b)))/2
    val = mp.quad(lambda ph: ltanh(ph)/mp.sqrt((xph(ph)+a)*(xph(ph)+b)), [0, mp.pi/2, mp.pi])/I0
    check("nome bound is an equality for H#_{pi/2} [b/a=%s]" % rb,
          abs(val + mp.pi*K(a/b)/K(mp.sqrt(1-(a/b)**2))) < 1e-12)

# sec:pointwise, pointwise bounds without zeros
def Rc(z, a, b):
    z = mp.mpc(z)
    return z*z*mp.sqrt(1-a*a/(z*z))*mp.sqrt(1-b*b/(z*z))

def pw_setup(a, b):
    k = mp.sqrt(1-a*a/(b*b))
    t0 = 2*K(k)/(mp.pi*b)
    m2 = b*b*E(k)/K(k)
    rr = lambda x: mp.sqrt((x*x-a*a)*(b*b-x*x))
    def A(x):
        if x <= a or x >= b: return mp.mpf(0)
        f = lambda ph: 1/(((a+b)/2-(b-a)/2*mp.cos(ph)+x)*mp.sqrt(((a+b)/2-(b-a)/2*mp.cos(ph)+a)*(b+(a+b)/2-(b-a)/2*mp.cos(ph))))
        return 2*rr(x)/mp.pi*mp.quad(f, [0, mp.pi])
    def hc(z):
        z = mp.mpc(z)
        Y = abs(z) + b
        return t0*mp.quad(lambda s: (s*s-m2)/Rc(s, a, b), [0, mp.mpc(0, Y), mp.mpc(z.real, Y), z])
    def cauchy(kf, z):
        f = lambda ph: (lambda s: 2*z*kf(s)/(s*s-z*z)*(b-a)/2*mp.sin(ph))((a+b)/2-(b-a)/2*mp.cos(ph))
        return mp.quad(f, [0, mp.pi])/mp.pi
    return k, t0, m2, rr, A, hc, cauchy

for (a, b) in [(1, 2), (1, mp.mpf('1.1')), (mp.mpf('0.3'), 5)]:
    a, b = mp.mpf(a), mp.mpf(b)
    k, t0, m2, rr, A, hc, cauchy = pw_setup(a, b)
    tag = " [a=%s, b=%s]" % (mp.nstr(a, 3), mp.nstr(b, 3))
    kap = (b-a)/(b+a)
    Lam = (K(kap)-E(kap))/kap
    Lm = (E(kap)-(1-kap**2)*K(kap))/kap
    check("R(0) = -ab and R ~ z^2"+tag, abs(Rc(mp.mpc(0, 1e-30), a, b)+a*b) < 1e-20
          and abs(Rc(mp.mpc(0, 1e6), a, b)/mp.mpc(0, 1e6)**2-1) < 1e-10)
    check("R(x+i0) = i r(x) on the band"+tag,
          abs(Rc(mp.mpc((a+b)/2, 1e-30), a, b) - 1j*rr((a+b)/2)) < 1e-15)
    for z in [mp.mpc((a+b)/2, (b-a)/3), mp.mpc(a/2, a/4), mp.mpc(2*b, b)]:
        k0 = lambda s: rr(s)/(a*b*s)
        lhs = Rc(z, a, b)/(a*b*z)
        rhs = z/(a*b) - 1/z + cauchy(k0, z)
        check("lem:atoms identity for K_0"+tag, abs(lhs-rhs) < 1e-15)
        for t in [a/2, (a+3*b)/2]:
            Rt = Rc(mp.mpf(t)+mp.mpc(0, 1e-40), a, b).real
            kt = lambda s: 2*s*rr(s)/(Rt*(t*t-s*s))
            lhs = 2*z*Rc(z, a, b)/(Rt*(t*t-z*z))
            rhs = -2*z/Rt + 2*z/(t*t-z*z) + cauchy(kt, z)
            check("lem:atoms identity for K_t, t=%s" % mp.nstr(t, 3)+tag, abs(lhs-rhs) < 1e-15)
    # kernels positive on the band
    xs = [a+(b-a)*j/7 for j in range(1, 7)]
    okk = True
    for t in [a/3, a*mp.mpf('0.99'), b*mp.mpf('1.01'), 4*b]:
        Rt = Rc(mp.mpf(t)+mp.mpc(0, 1e-40), a, b).real
        okk = okk and all(2*x*rr(x)/(Rt*(t*t-x*x)) > 0 for x in xs)
    check("k_t > 0 on the band for t off the band"+tag, okk)
    # beta identity: (2/pi) int 2t/(r (t0^2-t^2)) = 2/R(t0)
    for t in [a/2, 3*b]:
        Rt = Rc(mp.mpf(t)+mp.mpc(0, 1e-40), a, b).real
        val = 2/mp.pi*mp.quad(lambda ph: (lambda s: 2*s/((t*t-s*s)*mp.sqrt((s+a)*(b+s))))((a+b)/2-(b-a)/2*mp.cos(ph)), [0, mp.pi])
        check("eq:betarec for one resonance, t=%s" % mp.nstr(t, 3)+tag, abs(val-2/Rt) < 1e-15)
    # h_c on the band and the lead identity
    for z in [mp.mpc((a+b)/2, (b-a)/5), mp.mpc(a/2, a/3), mp.mpc(3*b, b/2)]:
        lhs = hc(-a*b/z); rhs = -hc(z) + t0*Rc(z, a, b)/z
        check("lem:leadprofile h(-ab/z) = -h(z) + tau R/z"+tag, abs(lhs-rhs) < 1e-12)
    for x in [a+(b-a)/5, (a+b)/2, b-(b-a)/7]:
        check("A(ab/x) = tau r/x - A(x) > 0"+tag, abs(A(a*b/x) - (t0*rr(x)/x - A(x))) < 1e-15 and A(a*b/x) > 0)
    check("int_a^b r/x = pi (b-a)^2/4"+tag, abs(mp.quad(lambda x: rr(x)/x, [a, b]) - mp.pi*(b-a)**2/4) < 1e-15)
    check("mean of A(ab/x) is Lambda_-"+tag, abs(mp.quad(lambda x: A(a*b/x), [a, mp.sqrt(m2), b])/(b-a) - Lm) < 1e-12)
    check("mean of A is Lambda"+tag, abs(mp.quad(A, [a, mp.sqrt(m2), b])/(b-a) - Lam) < 1e-12)
    # affine optimum below the threshold: Herglotz, linear coefficient, mean
    tau = -mp.mpf(3)/2*t0; C = -tau - t0
    hopt = lambda z: hc(z) + tau*z + C*Rc(z, a, b)/z
    pts = [mp.mpc(x, y) for x in [a/2, (a+b)/2, 2*b, -(a+b)/2] for y in [(b-a)/50, (b-a)/3, 2*b]]
    check("thm:affine optimum has Im h >= 0 in C+"+tag, all(mp.im(hopt(z)) > 0 for z in pts))
    y = mp.mpf(10)**6
    check("thm:affine optimum has linear coefficient 0 at the threshold"+tag,
          abs(hopt(mp.mpc(0, y))/mp.mpc(0, y)) < 1e-8)
    check("cor:pointwisecons(iv) mean"+tag,
          abs(mp.quad(lambda x: A(x) + C*rr(x)/x, [a, mp.sqrt(m2), b])/(b-a) - (Lam + mp.pi/4*(b-a)*C)) < 1e-12)

a, b = mp.mpf(1), mp.mpf(2)
k, t0, m2, rr, A, hc, cauchy = pw_setup(a, b)
kap = mp.mpf(1)/3
Lam = (K(kap)-E(kap))/kap; Lm = (E(kap)-(1-kap**2)*K(kap))/kap
check("octave A(m) = 0.34824", abs(A(mp.sqrt(m2)) - mp.mpf('0.348236277699756')) < 1e-12)
check("octave A(m)/Lambda = 1.2731", abs(A(mp.sqrt(m2))/Lam - mp.mpf('1.2731')) < 5e-5)
check("octave Lambda_- = 0.26560", abs(Lm - mp.mpf('0.265596407637276')) < 1e-12)
check("octave theta_c = 1.5477", abs(mp.pi*Lm/(Lam+Lm) - mp.mpf('1.5477')) < 5e-5)
check("octave tau_circ a = 0.6864", abs(t0 - mp.mpf('0.6864')) < 5e-5)


import random
# ssec:closedform, the optimum in closed form
mp.mp.dps = 20
def ell_setup(th, a, b):
    k = mp.sqrt(1-a*a/(b*b)); u = th*K(k)/mp.pi
    eta = b*mp.ellipfun('cn', u, k**2)/mp.ellipfun('sn', u, k**2)
    Q = mp.sqrt((eta**2+a*a)*(eta**2+b*b)); n = (b*b-a*a)/(b*b+eta**2)
    c = mp.ellippi(n, k**2)/((b*b+eta**2)*K(k))
    def al(x):
        psi = mp.asin(mp.sqrt((b*b-x*x)/(b*b-a*a)))
        return 2*eta*Q/b*(c*mp.ellipf(psi, k**2) - mp.ellippi(n, psi, k**2)/(b*b+eta**2))
    return k, eta, Q, c, al
for (a, b) in [(mp.mpf(1), mp.mpf(2)), (mp.mpf(1), mp.mpf(5))]:
    tag = " [b/a=%s]" % mp.nstr(b/a, 2)
    r = lambda t: mp.sqrt((t*t-a*a)*(b*b-t*t))
    for th in [mp.mpf('0.5'), mp.pi/2, mp.mpf('2.5')]:
        k, eta, Q, c, al = ell_setup(th, a, b)
        t2 = tag + " theta=%s" % mp.nstr(th, 3)
        e1 = etastar(th, a, b)
        check("thm:etaexplicit eta* = b cn/sn"+t2, abs(e1-eta) < 1e-9*eta)
        xs = [a+(b-a)*mp.mpf(j)/5 for j in (1, 2, 4)]
        check("thm:ellipticprofile alpha# equals v_g"+t2, max(abs(al(x)-vg(x, th, eta, a, b)) for x in xs) < 1e-10)
        Pint = lambda x: 2*x*r(x)*mp.quad(lambda u: 1/((x*x+u*u)*mp.sqrt((u*u+a*a)*(u*u+b*b))), [0, eta])
        k0 = mp.sqrt(1-a*a/(b*b)); t0 = 2*K(k0)/(mp.pi*b)
        Aq = lambda x: r(x)/x*2/mp.pi*mp.quad(lambda t: x/((x+t)*r(t)), [a, b])
        check("eq:asharpP alpha# = 2xr int - (pi-theta) A"+t2, max(abs(al(x)-(Pint(x)-(mp.pi-th)*Aq(x))) for x in xs) < 1e-10)
        e, L, Lz = Lall(th, a, b)
        check("eq:ellipticmean pi eta (Qc-1)/(b-a) = Lambda_all"+t2, abs(mp.pi*eta*(Q*c-1)/(b-a) - L) < 1e-10)
        xp = mp.sqrt(1/c-eta**2)
        check("peak x_theta: derivative vanishes, a < x < b"+t2, a < xp < b and abs(mp.diff(al, xp)) < 1e-8)
        x = (a+2*b)/3
        check("eq:alphaellipticprime"+t2, abs(mp.diff(al, x) - 2*eta*Q/r(x)*(1/(x*x+eta**2)-c)) < 1e-10)
        # transfer function: compare with -B e^{i h} at a point near 0, via the log derivative integral
        z = mp.mpc(a/3, a/5); ph = mp.asin(z/a)
        Hel = mp.exp(2j*Q/(b*eta)*(mp.ellippi(-a*a/eta**2, ph, (a/b)**2) - eta**2*c*mp.ellipf(ph, (a/b)**2)))
        Rz = lambda w: -mp.sqrt(a*a-w*w)*mp.sqrt(b*b-w*w)
        Hint = mp.exp(mp.quad(lambda w: 2j*eta*Q/Rz(w*z)*(c-1/((w*z)**2+eta**2))*z, [0, 1]))
        check("eq:transferelliptic"+t2, abs(Hel-Hint) < 1e-12)
    k, eta, Q, c, al = ell_setup(mp.pi/2, a, b)
    check("closed form at pi/2: eta = sqrt(ab), c = 1/(2ab)"+tag, abs(eta-mp.sqrt(a*b)) < 1e-15 and abs(c-1/(2*a*b)) < 1e-15)
a, b = mp.mpf(1), mp.mpf(2)
for th, mean in [(mp.mpf('0.5'), '0.13090'), (mp.mpf(1), '0.22860'), (mp.pi/2, '0.26951'), (mp.mpf(2), '0.24357'), (mp.mpf('2.5'), '0.15941')]:
    k, eta, Q, c, al = ell_setup(th, a, b)
    check("tab:allzeros mean %s from the closed form" % mean, abs(mp.pi*eta*(Q*c-1) - mp.mpf(mean)) < 5e-6)
# lem:zeroatom
r = lambda t: mp.sqrt((t*t-a*a)*(b*b-t*t))
Rf = lambda z: z*z*mp.sqrt(1-a*a/(z*z))*mp.sqrt(1-b*b/(z*z))
for u in [mp.mpf('0.3'), mp.mpf(2)]:
    Qu = mp.sqrt((u*u+a*a)*(u*u+b*b))
    phi = lambda z: z/(z*z+u*u)*(1+Rf(z)/Qu)
    cau = lambda z: 2*z*Rf(z)/mp.pi*mp.quad(lambda t: t/(t*t+u*u)/(r(t)*(z*z-t*t)), [a, b])
    zz = [mp.mpc(3, 1), mp.mpc(mp.mpf('0.5'), 2), mp.mpc(-1, mp.mpf('0.4'))]
    check("lem:zeroatom h_{j_u} formula [u=%s]" % u, max(abs(phi(z)-cau(z)) for z in zz) < 1e-10)
    check("lem:zeroatom beta = 1/Q_u [u=%s]" % u, abs(2/mp.pi*mp.quad(lambda t: t/((t*t+u*u)*r(t)), [a, b]) - 1/Qu) < 1e-10)
    check("R(iu) = -Q_u [u=%s]" % u, abs(Rf(mp.mpc(0, u))+Qu) < 1e-14)

# sec:beyond, a lead with a delay
mp.mp.dps = 20
def ex_setup(a, b):
    r = lambda t: mp.sqrt((t*t-a*a)*(b*b-t*t))
    k = mp.sqrt(1-a*a/(b*b)); t0 = 2*K(k)/(mp.pi*b)
    kap = (b-a)/(b+a); Lam = (K(kap)-E(kap))/kap; Lm = (E(kap)-(1-kap**2)*K(kap))/kap
    tq = lambda ph: (a+b)/2-(b-a)/2*mp.cos(ph)
    def q(y): return 2/mp.pi*mp.quad(lambda ph: mp.sqrt(y)/(mp.sqrt(y)+tq(ph))/mp.sqrt((tq(ph)+a)*(tq(ph)+b)), [0, mp.pi])
    def dq(y): return 2/mp.pi*mp.quad(lambda ph: tq(ph)/(2*mp.sqrt(y)*(mp.sqrt(y)+tq(ph))**2)/mp.sqrt((tq(ph)+a)*(tq(ph)+b)), [0, mp.pi])
    def d2q(y):
        x = mp.sqrt(y)
        return -2/mp.pi*mp.quad(lambda ph: tq(ph)*(3*x+tq(ph))/(4*x**3*(x+tq(ph))**3)/mp.sqrt((tq(ph)+a)*(tq(ph)+b)), [0, mp.pi])
    f, l = (a+b)/2, (b-a)/2; G0 = mp.pi*l*l
    Dl = lambda v: mp.sqrt((v-a*a)*(v-b*b))
    xv = lambda v: (v-Dl(v)+a*b)/(a+b)
    Gp = lambda v: mp.pi/2*(v-(a*a+b*b)/2-Dl(v))
    Gm = lambda u: mp.pi/2*((a*a+b*b)/2-u-mp.sqrt((a*a-u)*(b*b-u)))
    return r, t0, Lam, Lm, q, dq, d2q, f, l, G0, xv, Gp, Gm
def ex_solve(s, tau, a, b, S):
    r, t0, Lam, Lm, q, dq, d2q, f, l, G0, xv, Gp, Gm = S
    Phi = lambda v: t0-q(xv(v)**2)+dq(xv(v)**2)*(v-xv(v)**2)
    lo, hi = b*b*(1+mp.mpf('1e-12')), b*b*mp.mpf(10)**6
    for _ in range(200):
        mid = (lo+hi)/2
        if s*Phi(mid) < tau: lo = mid
        else: hi = mid
    v = (lo+hi)/2; xs = xv(v); ys = xs**2
    D = s*dq(ys)*(v-ys)**2/v; C = s*(q(ys)-dq(ys)*ys*(v-ys)/v)
    return v, xs, ys, C, D, -s*Lam+(C*G0+D*Gp(v))/(b-a)
for (a, b) in [(mp.mpf(1), mp.mpf(2)), (mp.mpf(1), mp.mpf(4))]:
    tag = " [b/a=%s]" % mp.nstr(b/a, 2)
    S = ex_setup(a, b); r, t0, Lam, Lm, q, dq, d2q, f, l, G0, xv, Gp, Gm = S
    xs3 = [a+(b-a)*mp.mpf(j)/6 for j in (1, 3, 5)]
    pw = pw_setup(a, b); A = pw[4]; hc = pw[5]
    check("lem:q q(x^2) = x A(x)/r(x)"+tag, max(abs(q(x*x)-x*A(x)/r(x)) for x in xs3) < 1e-12)
    check("lem:q q(y) -> tau_circ, q' > 0, q'' < 0"+tag, abs(q(mp.mpf(10)**16)-t0) < 1e-7
          and all(dq(y) > 0 and d2q(y) < 0 for y in [a*a/2, a*b, b*b, 3*b*b]))
    check("lem:q q' and q'' agree with numerical derivatives"+tag,
          abs(dq(a*b)-mp.diff(q, a*b)) < 1e-12 and abs(d2q(a*b)-mp.diff(q, a*b, 2)) < 1e-9)
    s = mp.mpf(1)
    tc = s*(t0-q(b*b))
    check("tau_c = (2s/pi) int t/((b+t) r) > 0"+tag, abs(tc - 2*s/mp.pi*mp.quad(lambda t: t/((b+t)*r(t)), [a, b])) < 1e-10 and tc > 0)
    check("Gamma_- and Gamma_+ closed forms"+tag,
          abs(Gm(a*a/3) - mp.quad(lambda x: x*r(x)/(x*x-a*a/3), [a, b])) < 1e-14
          and abs(Gp(3*b*b) - mp.quad(lambda x: x*r(x)/(3*b*b-x*x), [a, b])) < 1e-14 and abs(Gm(0)-G0) < 1e-20)
    check("M_s(0) = s Lambda_-"+tag, abs(-s*Lam+mp.pi/4*(b-a)*s*t0 - s*Lm) < 1e-14)
    # (i): below the threshold the origin-only optimum is Herglotz and attains the bound
    Rcz = lambda z: Rc(z, a, b)
    for tau in [-s*t0/2, tc/2, tc]:
        ho = lambda z: -s*hc(z)+tau*z+(s*t0-tau)*Rcz(z)/z
        pts = [mp.mpc(x, y) for x in [a/2, (a+b)/2, 2*b, -(a+b)/2] for y in [(b-a)/50, (b-a)/3, 2*b]]
        check("thm:exactmean(i) optimum Herglotz, band density = bound [tau=%s]" % mp.nstr(tau, 4) + tag,
              all(mp.im(ho(z)) > -1e-15 for z in pts) and
              max(abs(mp.im(ho(mp.mpc(x, mp.mpf('1e-12'))))-r(x)/x*(s*t0-tau-s*q(x*x))) for x in xs3) < 1e-9)
    for tau in [2*tc, mp.mpf(3)]:
        v, xs, ys, C, D, M = ex_solve(s, tau, a, b, S)
        t2 = tag + " tau=%s" % mp.nstr(tau, 4)
        al = lambda x: r(x)/x*(C+D*x*x/(v-x*x)-s*q(x*x))
        grid = [a+(b-a)*mp.mpf(j)/40 for j in range(1, 40)]
        check("thm:exactmean(ii) C, D > 0, f < x* < b"+t2, C > 0 and D > 0 and f < xs < b)
        check("thm:exactmean(ii) band density >= 0, zero at x*"+t2, min(al(x) for x in grid) >= 0 and abs(al(xs)) < 1e-12)
        check("thm:exactmean(ii) linear coefficient 0"+t2, abs(tau-s*t0+C-D) < 1e-12)
        ho = lambda z: -s*hc(z)+tau*z+Rcz(z)*(C/z+D*z/(v-z*z))
        pts = [mp.mpc(x, y) for x in [a/2, (a+b)/2, 2*b, mp.sqrt(v), -(a+b)/2] for y in [(b-a)/50, (b-a)/3, 2*b]]
        check("thm:exactmean(ii) optimum is Herglotz"+t2, all(mp.im(ho(z)) > 0 for z in pts))
        check("eq:exactmean equals the integral of the profile"+t2, abs(mp.quad(al, [a, xs, b])/(b-a) - M) < 1e-12)
        l1 = mp.pi*(xs-a)*(b-xs); l2 = mp.pi*(xs-f)**2
        check("dual form of the mean"+t2, abs((b-a)*M - (-s*(b-a)*Lam+l1*s*q(ys)+l2*(s*t0-tau))) < 1e-12)
        check("lem:meandual lambda1+lambda2 = Gamma0, low and high slacks"+t2,
              abs(l1+l2-G0) < 1e-15
              and min(Gm(a*a*mp.mpf(j)/30)-l1*ys/(ys-a*a*mp.mpf(j)/30)-l2 for j in range(1, 30)) > 0
              and min(Gp(b*b*(1+mp.mpf(j)**2/20))+l2-l1*ys/(b*b*(1+mp.mpf(j)**2/20)-ys) for j in range(1, 60)) > -1e-15
              and abs(Gp(v)+l2-l1*ys/(v-ys)) < 1e-12)
        # random admissible measures: masses at 0, below and above the band
        random.seed(7)
        qg = {x: q(x*x) for x in grid}
        cnt = 0; worst = mp.inf
        for trial in range(300):
            c0 = 3*mp.mpf(random.random())
            lows = [(a*mp.mpf(random.random()), mp.mpf(random.random())) for _ in range(random.randint(0, 2))]
            highs = [(b*(1+6*mp.mpf(random.random())), 6*mp.mpf(random.random())) for _ in range(random.randint(0, 3))]
            Mf = lambda y: c0+sum(m*y/(y-t*t) for t, m in lows)+sum(m*y/(t*t-y) for t, m in highs)
            beta = tau-s*t0+c0+sum(m for t, m in lows)-sum(m for t, m in highs)
            if beta < 0 or any(Mf(x*x) < s*qg[x] for x in grid): continue
            cost = (-s*(b-a)*Lam+c0*G0+sum(m*Gm(t*t) for t, m in lows)+sum(m*Gp(t*t) for t, m in highs))/(b-a)
            cnt += 1; worst = min(worst, cost-M)
        check("random admissible measures cost at least the optimum (%d)" % cnt + t2, cnt > 20 and worst > 0)
# the square identity, also continued to negative w
a, b = mp.mpf(1), mp.mpf(3); f, l = 2, 1
for xs in [mp.mpf('2.2'), mp.mpf('2.8')]:
    ys = xs*xs; l1 = mp.pi*(xs-a)*(b-xs); l2 = mp.pi*(xs-f)**2; ep = xs-f
    for w in [mp.mpf(10), mp.mpf(40), mp.mpf(-1), mp.mpf(-7)]:
        sq = mp.sqrt((a*a-w)*(b*b-w)) if w < 0 else mp.sqrt((w-a*a)*(w-b*b))
        xw = (w+sq+a*b)/(a+b) if w < 0 else (w-sq+a*b)/(a+b)
        rh = xw-f
        Gw = -mp.pi/2*((a*a+b*b)/2-w-sq) if w < 0 else mp.pi/2*(w-(a*a+b*b)/2-sq)
        ok = abs(w-(f*f+l*l+f*(rh+l*l/rh))) < 1e-15 and abs(Gw-mp.pi*f*rh) < 1e-15
        ok = ok and abs((w-ys)*(Gw+l2-l1*ys/(w-ys)) - mp.pi*f*(f+l*l/rh)*(rh-ep)**2) < 1e-13
        if w < 0: ok = ok and -l*l/f < rh < 0
        check("eq:square at x*=%s, w=%s" % (xs, w), ok)
# cor:delaybranches monotonicity in s
a, b = mp.mpf(1), mp.mpf(2); S = ex_setup(a, b)
Ms = [ex_solve(s, mp.mpf(5), a, b, S)[5] for s in [mp.mpf('1.2'), mp.pi/2, mp.mpf(2)]]
check("cor:delaybranches M_s(5) increases with s", Ms[0] < Ms[1] < Ms[2])

# ssec:delay, the counterexample on [1,2]
r, t0, Lam, Lm, q, dq, d2q, f, l, G0, xv, Gp, Gm = S
s = mp.pi/2; T0 = 16
pw = pw_setup(a, b); hc = pw[5]
check("tau_circ <= 1 and >= 1/2 on [1,2]", mp.mpf(1)/2 <= t0 <= 1)
grid = [1+mp.mpf(j)/20 for j in range(0, 21)]
check("q' <= tau/8, -q'' <= tau/8 on [1,4]", all(dq(x*x) <= t0/8 and -d2q(x*x) <= t0/8 for x in grid))
check("q' >= 1/128, -q'' >= 1/1024 on [1,4]", all(dq(x*x) >= mp.mpf(1)/128 and -d2q(x*x) >= mp.mpf(1)/1024 for x in grid))
check("q - y q' >= 3/32 on [9/4,4], q(4) >= tau/2", all(q(x*x)-x*x*dq(x*x) >= mp.mpf(3)/32 for x in grid if x >= 1.5) and q(4) >= t0/2)
for y0 in [mp.mpf('1.3'), mp.mpf('2.5'), mp.mpf('3.7')]:
    d = dq(y0)*(T0-y0)**2/T0; c = q(y0)-dq(y0)*y0*(T0-y0)/T0
    al = lambda x: s*r(x)/x*(c+d*x*x/(T0-x*x)-q(x*x))
    hy = lambda z: -s*hc(z)+5*z+s*Rc(z, a, b)*(c/z+d*z/(T0-z*z))
    pts = [mp.mpc(x, yy) for x in [mp.mpf('0.5'), mp.mpf('1.5'), 4, 7, -2] for yy in [mp.mpf('0.02'), mp.mpf('0.4'), 3]]
    check("thm:delay witness y0=%s: c, d > 0, Herglotz, alpha(sqrt y0) = 0, alpha >= 0" % y0,
          c > 0 and d > 0 and all(mp.im(hy(z)) > 0 for z in pts) and abs(al(mp.sqrt(y0))) < 1e-14
          and min(al(x) for x in grid) >= -1e-20 and 5-s*t0+s*(c-d) > 0)
y0 = mp.mpf(5)/2
d = dq(y0)*(T0-y0)**2/T0; c = q(y0)-dq(y0)*y0*(T0-y0)/T0
al = lambda x: s*r(x)/x*(c+d*x*x/(T0-x*x)-q(x*x))
M52 = mp.quad(al, [1, mp.sqrt(y0), 2])
check("mean of H_{5/2} = 0.0084888 and the exact formula", abs(M52-mp.mpf('0.0084888')) < 5e-8
      and abs(M52 - s*(-Lam+mp.pi*c/4+mp.pi*d/2*(16-mp.mpf(5)/2-mp.sqrt(15*12)))) < 1e-14)
check("0 <= f0 <= (7/96)(y-y0)^2 on [1,4]", all(0 <= c+d*x**2/(T0-x**2)-q(x**2) <= mp.mpf(7)/96*(x*x-y0)**2+1e-25 for x in grid))
check("eq:delaybound 21 pi^2/1024 < (pi/2)(3-2sqrt2)", M52 <= 21*mp.pi**2/1024 < mp.pi/2*(3-2*mp.sqrt(2)))
check("Lambda_all(pi/2) on [1,2] = 0.26951", abs(mp.pi/2*(3-2*mp.sqrt(2)) - mp.mpf('0.26951')) < 5e-6)
v, xs, ys, C, D, M = ex_solve(s, mp.mpf(5), a, b, S)
check("M_{pi/2}(5) parameters sqrt v, x*, C, D",
      abs(mp.sqrt(v)-mp.mpf('8.870977')) < 5e-7 and abs(xs-mp.mpf('1.504922')) < 5e-7
      and abs(C-mp.mpf('0.426778')) < 5e-7 and abs(D-mp.mpf('4.348520')) < 5e-7)
check("M_{pi/2}(5) = 0.00639032, ratio about 42", abs(M-mp.mpf('0.00639032')) < 5e-9 and abs(mp.pi/2*(3-2*mp.sqrt(2))/M-42) < 0.5)
check("zero-free branch: M <= 21pi^2/1024 < pi^2/24 <= (pi/2)Lambda", M < 21*mp.pi**2/1024 < mp.pi**2/24 <= mp.pi/2*Lam)
check("tau_c on the octave, s=1, is 0.2852", abs(t0-q(4) - mp.mpf('0.2852')) < 5e-5)

# ssec:delayzeros, one zero
Qu = lambda u: mp.sqrt((1+u*u)*(4+u*u))
def onezero(eps):
    P = lambda y: 2*mp.quad(lambda u: y/((y+u*u)*Qu(u)), [0, eps])
    dP = lambda y: 2*mp.quad(lambda u: u*u/((y+u*u)**2*Qu(u)), [0, eps])
    qe = lambda y: s*q(y)-P(y); dqe = lambda y: s*dq(y)-dP(y)
    be = 5-s*t0+2*mp.quad(lambda u: 1/Qu(u), [0, eps])
    lo, hi = 4*(1+mp.mpf('1e-12')), mp.mpf(10)**6
    for _ in range(200):
        mid = (lo+hi)/2
        if dqe(xv(mid)**2)*(mid-xv(mid)**2)-qe(xv(mid)**2) < be: lo = mid
        else: hi = mid
    v1 = (lo+hi)/2; y1 = xv(v1)**2
    D1 = dqe(y1)*(v1-y1)**2/v1; C1 = qe(y1)-dqe(y1)*y1*(v1-y1)/v1
    M1 = -s*Lam+2*mp.quad(lambda u: Gm(-u*u)/Qu(u), [0, eps])+C1*G0+D1*Gp(v1)
    al1 = lambda x: r(x)/x*(C1+D1*x*x/(v1-x*x)-qe(x*x))
    return M1, C1, D1, v1, y1, al1, be, qe, dqe, P, dP
eps = mp.mpf(1)/10
M1, C1, D1, v1, y1, al1, be, qe, dqe, P, dP = onezero(eps)
check("prop:onezerodelay P' <= eps^3/3, -P'' <= 2eps^3/3, P <= eps on [1,4]",
      all(0 <= dP(x*x) <= eps**3/3 and P(x*x) <= eps for x in grid)
      and all(-mp.diff(P, x*x, 2) <= 2*eps**3/3 for x in grid[::5]))
check("prop:onezerodelay C, D > 0, beta > 0, x* in (3/2,2)", C1 > 0 and D1 > 0 and be > 0 and 1.5 < mp.sqrt(y1) < 2)
check("prop:onezerodelay band density >= 0", min(al1(x) for x in grid[1:-1]) >= 0)
check("prop:onezerodelay linear coefficient 0", abs(be+C1-D1) < 1e-12)
check("prop:onezerodelay mean = 0.0063752 < M_{pi/2}(5)", abs(M1-mp.mpf('0.0063752')) < 5e-8 and M1 < M)
check("one-zero mean formula equals the integral of the profile",
      abs(M1 - mp.quad(al1, [1, mp.sqrt(y1), 2])) < 1e-12)
xs1 = mp.sqrt(y1); l1 = mp.pi*(xs1-1)*(2-xs1); l2 = mp.pi*(xs1-f)**2
Lj = lambda u: (Gm(-u*u)-l1*y1/(y1+u*u)-l2)/Qu(u)
check("L(j_u) < 0 for u in (0, eps]", all(Lj(eps*mp.mpf(j)/10) < 0 for j in range(1, 11)))
check("arg(-B_eps) = pi - p_eps on the band", abs(mp.arg(-(mp.mpf('1.3')-1j*eps)/(mp.mpf('1.3')+1j*eps)) - (mp.pi-2*mp.atan(eps/mp.mpf('1.3')))) < 1e-15)
# rem:onezerofamily
lo, hi = mp.mpf('0.3'), mp.mpf('0.6')
for _ in range(40):
    mid = (lo+hi)/2
    if onezero(mid)[1] > 0: lo = mid
    else: hi = mid
e0 = (lo+hi)/2; M0 = onezero(e0)[0]
check("rem:onezerofamily eps0 = 0.47247, mean 0.0051493", abs(e0-mp.mpf('0.47247')) < 5e-5 and abs(M0-mp.mpf('0.0051493')) < 5e-7)
check("rem:onezerofamily mean decreases in eps", onezero(mp.mpf('0.2'))[0] < M1 and onezero(mp.mpf('0.4'))[0] < onezero(mp.mpf('0.2'))[0])
# thm:poisson
best = mp.mpf(0)
for j in range(1, 60):
    y = mp.mpf(j)/100
    om = (mp.atan(2/y)-mp.atan(1/y))/mp.pi
    ep = (mp.exp(-5*y)/2**(1-om))**(1/om)
    best = max(best, om/(y/(mp.pi*(1+y*y)))*(-mp.log(1-ep)))
check("thm:poisson floor for the example is positive", 1e-22 < best < 1e-18)

def delay_checks():
    dps0 = mp.mp.dps
    # prop:floor, the zero-free floor
    mp.mp.dps = 20
    for (a, b) in [(mp.mpf(1), mp.mpf(2)), (mp.mpf(1), mp.mpf(4))]:
        tag = " [b/a=%s]" % mp.nstr(b/a, 2)
        S = ex_setup(a, b); r, t0, Lam, Lm, q, dq, d2q, f, l, G0, xv, Gp, Gm = S
        check("prop:floor int x r = f^2 Gamma_0"+tag, abs(mp.quad(lambda x: x*r(x), [a, b]) - f*f*G0) < 1e-15)
        minf = mp.quad(lambda x: r(x)/x*(q(f*f)-q(x*x)), [a, f, b])/(b-a)
        check("prop:floor m_inf > 0 and (b-a) m_inf = Gamma_0 q(f^2) - (b-a) Lambda"+tag,
              minf > 0 and abs((b-a)*minf - (G0*q(f*f)-(b-a)*Lam)) < 1e-12)
        check("prop:floor Gamma_-(u) >= Gamma_0 f^2/(f^2-u), Gamma_+(w) >= Gamma_0 f^2/(w-f^2)"+tag,
              all(Gm(u) >= G0*f*f/(f*f-u)-1e-18 for u in [0, a*a/3, a*a*mp.mpf('0.99')])
              and all(Gp(w) >= G0*f*f/(w-f*f)-1e-18 for w in [b*b*mp.mpf('1.01'), 3*b*b, 100*b*b]))
    a, b = mp.mpf(1), mp.mpf(2)
    S = ex_setup(a, b); r, t0, Lam, Lm, q, dq, d2q, f, l, G0, xv, Gp, Gm = S
    s = mp.pi/2
    minf = mp.quad(lambda x: r(x)/x*(q(f*f)-q(x*x)), [a, f, b])/(b-a)
    check("prop:floor (pi/2) m_inf = 0.0060483 on the octave", abs(s*minf - mp.mpf('0.0060483')) < 5e-8)
    Ms = [ex_solve(s, mp.mpf(t), a, b, S)[5] for t in [5, 10, 50, 10**4]]
    check("prop:floor M_{pi/2}(tau) decreases to (pi/2) m_inf (tau = 5, 10, 50, 1e4)",
          all(Ms[i] > Ms[i+1] for i in range(3)) and all(M >= s*minf for M in Ms) and Ms[3]-s*minf < 1e-6)
    check("prop:floor M_{pi/2}(5) is less than six per cent above the floor", Ms[0]/(s*minf) < mp.mpf('1.06'))

    # thm:expdecay, the constants
    check("thm:expdecay e^{3-12/e} < 1/4", mp.exp(3-12/mp.e) < mp.mpf(1)/4)
    check("thm:expdecay 2 sqrt m <= 3(m-1) for m >= 4", all(2*mp.sqrt(m) <= 3*(m-1) for m in range(4, 400, 2)))
    Imf = lambda m: 2*mp.quadosc(lambda u: (mp.sin(u)/u)**m if u != 0 else mp.mpf(1), [0, mp.inf], period=mp.pi)
    check("thm:expdecay I_m >= 4 sqrt6/(3 sqrt m) > 3/sqrt m (m = 4, ..., 30)",
          all(Imf(m) >= 4*mp.sqrt(6)/(3*mp.sqrt(m)) > 3/mp.sqrt(m) for m in range(4, 31, 2)))
    check("thm:expdecay (1-v^2)^m >= 1 - m v^2 and sin u/u >= 1 - u^2/6",
          all((1-v*v)**m >= 1-m*v*v for m in [4, 10, 50] for v in [mp.mpf(j)/20 for j in range(21)])
          and all(mp.sin(u)/u >= 1-u*u/6 for u in [mp.mpf(j)/10 for j in range(1, 25)]))
    check("thm:expdecay |sin w / w| <= e^{|Im w|}",
          all(abs(mp.sin(w)/w) <= mp.exp(abs(mp.im(w))) for w in [mp.mpc(x, y) for x in [0.3, 2, 7] for y in [-3, -0.5, 0.2, 4]]))
    check("thm:expdecay (1/2)(1-eps)^{-2} <= 8/9 at eps = 1/4", mp.mpf(1)/2/(mp.mpf(3)/4)**2 <= mp.mpf(8)/9)
    for (a, b) in [(mp.mpf(1), mp.mpf(2)), (mp.mpf(1), mp.mpf(4))]:
        tag = " [b/a=%s]" % mp.nstr(b/a, 2)
        r = lambda t: mp.sqrt((t*t-a*a)*(b*b-t*t))
        check("thm:expdecay int t/r = pi/2 and x r(x) <= b(b^2-a^2)/2"+tag,
              abs(mp.quad(lambda t: t/r(t), [a, b]) - mp.pi/2) < 1e-10
              and all(x*r(x) <= b*(b*b-a*a)/2 for x in [a+(b-a)*mp.mpf(j)/40 for j in range(41)]))
        for tau in [12, 20, 40]:
            tau = mp.mpf(tau)/a
            m = 2*int(mp.floor(a*tau/(2*mp.e))); kap = tau/m
            eps = 2*(kap*a)**(1-m)/((m-1)*Imf(m))
            check("thm:expdecay eps <= e^{3-a tau/e}, m in [4, a tau/e], kappa a >= e [a tau=%s]%s" % (mp.nstr(a*tau, 3), tag),
                  eps <= mp.exp(3-a*tau/mp.e) and 4 <= m <= a*tau/mp.e and kap*a >= mp.e)
    # the right-angle example on the octave
    mp.mp.dps = 15
    a, b = mp.mpf(1), mp.mpf(2)
    k = mp.sqrt(1-a*a/(b*b)); t0 = 2*K(k)/(mp.pi*b)
    def expo_mean(tau):
        m = 2*int(mp.floor(a*tau/(2*mp.e))); kap = tau/m; Im = Imf(m)
        tail = lambda x: 2/Im*mp.quadosc(lambda u: (mp.sin(u)/u)**m, [kap*x, mp.inf], period=mp.pi)
        xs = [a+(b-a)*mp.mpf(j)/20 for j in range(21)]
        vals = [-mp.log(1-tail(x)) for x in xs]
        h = (b-a)/20
        mean = h/3*(vals[0]+vals[-1]+4*sum(vals[1:-1:2])+2*sum(vals[2:-1:2]))/(b-a)
        bound = mp.exp(3-a*tau/mp.e)*(2+mp.pi/4*(b-a)*t0+b*(b*b-a*a)*tau/(2*mp.e*a*a))
        return mean, bound
    for tau, val in [(12, '1.7e-3'), (20, '4.9e-5'), (30, '1.4e-7')]:
        mean, bound = expo_mean(mp.mpf(tau))
        check("thm:expdecay theta = pi/2 on the octave: mean %s at tau = %d, below the bound" % (val, tau),
              abs(mean/mp.mpf(val)-1) < 0.04 and mean <= bound)
    mp.mp.dps = 20

    # cor:rate
    for (a, b) in [(mp.mpf(1), mp.mpf(2)), (mp.mpf(1), mp.mpf(4))]:
        tag = " [b/a=%s]" % mp.nstr(b/a, 2)
        y = mp.mpf('1e-4')
        om = (mp.atan(b/y)-mp.atan(a/y))/mp.pi
        check("cor:rate omega_y/y -> (b-a)/(pi a b), omega_y/P_y(a) -> a(b-a)/b"+tag,
              abs(om/y-(b-a)/(mp.pi*a*b)) < 1e-6 and abs(om/(y/(mp.pi*(a*a+y*y)))-a*(b-a)/b) < 1e-6)
    check("cor:rate on the octave pi ab/(b-a) = 2 pi", abs(mp.pi*2/(2-1)-2*mp.pi) < 1e-30)

    # rem:onezerofamily, the linear program
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import zero_lp
    import numpy as np
    check("rem:onezerofamily LP: 0.00658 without zeros, 4.3e-4 with the pair",
          abs(zero_lp.solve(np.pi/2, 5.0, [], kshift=-1)-0.00658) < 5e-6
          and abs(zero_lp.solve(np.pi/2, 5.0, [1.3654+1.2608j], kshift=-1)-4.3e-4) < 5e-6)

    mp.mp.dps = dps0
delay_checks()

def construction_checks():
    from fractions import Fraction as Fr
    dps0 = mp.mp.dps
    mp.mp.dps = 30
    # thm:convex (iii), the strong convexity of -log on (0,1]
    random.seed(7); okc = True
    for _ in range(200):
        r0, r1, t = [mp.mpf(random.uniform(1e-3, 1)) for _ in range(3)]
        okc = okc and -mp.log((1-t)*r0+t*r1) <= -(1-t)*mp.log(r0)-t*mp.log(r1)-t*(1-t)*(r1-r0)**2/2+mp.mpf(10)**-25
    check("thm:convex(iii) -log((1-t)r0+t r1) <= (1-t)(-log r0)+t(-log r1)-t(1-t)(r1-r0)^2/2 (random)", okc)
    # lem:toeplitz
    def toep(c):
        N = len(c); return mp.matrix([[c[i-j] if i >= j else 0 for j in range(N)] for i in range(N)])
    def taylor_blaschke(zs, N):
        # coefficients of prod (w - z)/(1 - conj(z) w) with real symmetric zeros
        co = [mp.mpf(1)]+[mp.mpf(0)]*N
        for z in zs:
            num = [-z, 1]; den = [mp.conj(z)**k for k in range(N+1)]   # 1/(1-conj z w)
            f = [sum(num[i]*den[k-i] for i in range(2) if 0 <= k-i <= N) for k in range(N+1)]
            co = [sum(co[i]*f[k-i] for i in range(k+1)) for k in range(N+1)]
        return [mp.re(x) for x in co]
    c = taylor_blaschke([mp.mpc('0.3', '0.5'), mp.mpc('0.3', '-0.5'), mp.mpf('-0.6')], 12)
    c = [x*mp.mpf('0.9') for x in c]
    nrm = max(mp.sqrt(abs(e)) for e in mp.eighe(toep(c).T*toep(c))[0])
    check("lem:toeplitz ||T_N(c)|| <= 1 and sum c_k^2 <= 1 for 0.9 x a Blaschke product (N = 12)", nrm <= 1 and sum(x*x for x in c) <= 1)
    N = len(c)-1; okf = True
    for th in [0.3, 1.1, 2.5, -1.7]:
        w = mp.expjpi(th/mp.pi); u = mp.matrix([mp.conj(w)**k for k in range(N+1)])
        lhs = (u.H*toep(c)*u)[0]/(N+1); rhs = sum((1-mp.mpf(k)/(N+1))*c[k]*w**k for k in range(N+1))
        okf = okf and abs(lhs-rhs) < 1e-25 and abs(rhs) <= 1
    check("lem:toeplitz u* T_N u/(N+1) = p_N(w) on |w| = 1, and |p_N| <= 1", okf)
    random.seed(11); okr = True
    for _ in range(20):
        cc = [mp.mpf(random.uniform(-1, 1)) for _ in range(8)]
        nr = max(mp.sqrt(abs(e)) for e in mp.eighe(toep(cc).T*toep(cc))[0]); cc = [x/nr for x in cc]
        okr = okr and max(abs(sum((1-mp.mpf(k)/8)*cc[k]*mp.expjpi(2*j/mp.mpf(97))**k for k in range(8))) for j in range(97)) <= 1+mp.mpf(10)**-20
    check("lem:toeplitz random contractions give Fejer polynomials bounded by 1 on the circle", okr)
    # lem:tail, |A_{I,k}| <= C_I/k for affine phases
    mp.mp.dps = 20; okt = True
    for (u, v, th, tau, sg) in [(1, 2, 1.0, 5, 1), (1, mp.mpf('1.5'), 0.4, 12, 2), (mp.mpf('1.5'), 2, 2.0, 0, mp.mpf('0.5'))]:
        u, v, sg = mp.mpf(u), mp.mpf(v), mp.mpf(sg)
        CI = (v*v+sg*sg)/sg+abs(tau)*((v**3-u**3)/(6*sg)+sg*(v-u)/2)
        for k in [1, 3, 10, 40]:
            A = mp.quad(lambda x: mp.expj(-(th+tau*x))*((x-1j*sg)/(x+1j*sg))**k, mp.linspace(u, v, 9))
            okt = okt and abs(A) <= CI/k
        jf = lambda x: (x*x+sg*sg)/(2*sg)
        okt = okt and abs(CI-(2*jf(v)+mp.quad(lambda x: jf(x)*abs(tau), [u, v]))) < 1e-15
    check("lem:tail |A_{I,k}| <= C_I/k, and the affine formula for C_I", okt)
    check("lem:tail sum_{k>N} k^-2 <= 1/N (N = 1, 5, 50)", all(mp.nsum(lambda k: 1/k**2, [N+1, mp.inf]) <= mp.mpf(1)/N for N in [1, 5, 50]))
    # thm:sharperate, the constants
    mp.mp.dps = 30
    def Cab(a, b):
        d = mp.pi/(a+b); s0 = mp.sin(d*a); return mp.atanh(s0)/d, d, s0
    okC = True
    for (a, b) in [(1, 2), (1, mp.mpf('1.1')), (mp.mpf('0.3'), 5), (1, 10)]:
        a, b = mp.mpf(a), mp.mpf(b); C, d, s0 = Cab(a, b)
        okC = okC and abs(C-(a+b)/mp.pi*mp.log(mp.cot(mp.pi*(b-a)/(4*(a+b))))) < 1e-25 and a/mp.e < a < C < mp.pi*a*b/(b-a)
    check("thm:sharperate C_{a,b} = gamma/d = (a+b)/pi log cot(pi(b-a)/(4(a+b))), a < C_{a,b} < pi ab/(b-a)", okC)
    C12 = Cab(mp.mpf(1), mp.mpf(2))[0]
    check("thm:sharperate C_{1,2} = (3/pi) log(2+sqrt3) = 1.2576021548", abs(C12-3/mp.pi*mp.log(2+mp.sqrt(3))) < 1e-25 and abs(C12-mp.mpf('1.2576021548')) < 1e-10)
    check("thm:sharperate artanh(sin t) = int_0^t sec > t", all(abs(mp.atanh(mp.sin(t))-mp.quad(mp.sec, [0, t])) < 1e-20 and mp.atanh(mp.sin(t)) > t for t in [0.1, 0.7, 1.3, 1.5]))
    # the kernels K_m, the normaliser I_m and Pi_m, as polynomials in s (high precision: their coefficients are large)
    mp.mp.dps = 160
    a, b = mp.mpf(1), mp.mpf(2); C, d, s0 = Cab(a, b); gam = mp.atanh(s0)
    def kernel(m):
        X = [(-1-s0**2)/(1-s0**2), 0, 2/(1-s0**2)]   # coefficients in s
        Tp, Tc = [mp.mpf(1)], X[:]
        mul = lambda p, q: [sum(p[i]*q[k-i] for i in range(len(p)) if 0 <= k-i < len(q)) for k in range(len(p)+len(q)-1)]
        if m == 0: T = Tp
        else:
            for _ in range(m-1):
                nx = mul([2*x for x in X], Tc); pad = Tp+[0]*(len(nx)-len(Tp))
                Tp, Tc = Tc, [nx[i]-pad[i] for i in range(len(nx))]
            T = Tc
        Km = mul(T, T); Prim = [0]+[Km[k]/(k+1) for k in range(len(Km))]
        Im = sum(Prim); return Km, Prim, Im
    ev = lambda p, s: sum(p[k]*s**k for k in range(len(p)))
    okk = True; logs = []
    for m in [1, 2, 3, 5, 8]:
        Km, Prim, Im = kernel(m); eps = (1-s0)/Im
        Pi = lambda s: ev(Prim, s)/Im
        okk = okk and len(Prim)-1 == 4*m+1 and abs(Pi(1)-1) < 1e-20 and abs(Pi(-1)+1) < 1e-20
        okk = okk and all(ev(Km, s) <= 1+mp.mpf(10)**-20 and -1e-20 <= 1-Pi(s) <= eps+1e-20 for s in mp.linspace(s0, 1, 25))
        okk = okk and Im <= mp.e**(4*m*gam)
        okk = okk and all(abs(ev(Km, s)-mp.cosh(m*mp.acosh((1+s0**2-2*s*s)/(1-s0**2)))**2) < 1e-15*ev(Km, s) for s in [0, s0/3, s0/2])
        logs.append(mp.log(Im)/m)
    check("thm:sharperate K_m <= 1 and 0 <= 1-Pi_m <= eps_m on [s0, 1], Pi_m(+-1) = +-1, degree 4m+1", okk)
    check("lem:chebnorm chi(0) = 2 gamma, I_m <= e^{4m gamma}, log I_m/m increases towards 4 gamma",
          abs(mp.acosh((1+s0**2)/(1-s0**2))-2*gam) < 1e-25 and all(logs[i] < logs[i+1] < 4*gam for i in range(len(logs)-1)))
    I80 = mp.quad(lambda s: mp.chebyt(80, (2*s*s-1-s0**2)/(1-s0**2))**2, [0, s0/4, s0/2, s0, 1])
    check("lem:chebnorm log I_80/80 (by quadrature) is within 0.06 of 4 gamma", 0 < 4*gam-mp.log(I80)/80 < mp.mpf('0.06'))
    # passivity of F_m: a polynomial in q = e^{idz} with real coefficients, bounded by 1 on the circle
    okp = True
    for m in [1, 2]:
        Km, Prim, Im = kernel(m); N = 4*m+1; th = mp.mpf('0.9')
        coef = [mp.mpc(0)]*(2*N+1)
        for k, pk in enumerate(Prim):
            if pk == 0: continue
            for j in range(k+1):   # ((q - 1/q)/(2i))^k
                cbin = mp.binomial(k, j)*(-1)**(k-j)/(2j)**k
                coef[N+j-(k-j)] += 1j*mp.sin(th)*pk/Im*cbin
        coef[N] += mp.cos(th)
        okp = okp and max(abs(mp.im(x)) for x in coef) < 1e-20
        okp = okp and max(abs(sum(coef[n]*mp.expjpi(2*j/mp.mpf(301))**n for n in range(2*N+1))) for j in range(301)) <= 1+mp.mpf(10)**-20
        z = mp.mpc('1.3', '0.2'); qq = mp.exp(1j*d*z)
        okp = okp and abs(sum(coef[n]*qq**n for n in range(2*N+1))-qq**N*(mp.cos(th)+1j*mp.sin(th)*ev(Prim, mp.sin(d*z))/Im)) < 1e-20
    check("lem:chebnorm q^N(cos th + i sin th Pi_m) is a real polynomial in q of degree 2N, bounded by 1 on |q| = 1", okp)
    # the mean bound of thm:sharperate at a right angle, and in general
    mp.mp.dps = 20
    t0 = ex_setup(a, b)[1]
    okm = True
    for m in [2, 3]:
        Km, Prim, Im = kernel(m); eps = (1-s0)/Im
        mean = mp.quad(lambda x: -mp.log(ev(Prim, mp.sin(d*x))/Im), [a, b])/(b-a)
        Dab = b*(b*b-a*a)/(4*a)*(d/a+(1-s0)/a**2)
        okm = okm and eps <= mp.mpf(1)/4 and mean <= 2*eps and 2*eps <= (2*Dab+(2+mp.pi/4*(b-a)*t0)*(1-s0))/Im
        okm = okm and max(d*mp.cos(d*x)*ev(Km, mp.sin(d*x))/Im for x in mp.linspace(a, b, 41)) <= d/Im
    check("thm:sharperate octave, m = 2, 3: eps_m <= 1/4, mean of -log U_m <= 2 eps_m <= the general bound, |U_m'| <= d/I_m", okm)
    # prop:sixtap
    mp.mp.dps = 30
    P = lambda s: (192*s**5-560*s**3+735*s)/367
    Hs = lambda z: (-12-80*mp.exp(2j*z)-435*mp.exp(4j*z)+435*mp.exp(6j*z)+80*mp.exp(8j*z)+12*mp.exp(10j*z))/734
    pts = [mp.mpc('1.3', '0.4'), mp.mpc('-0.7', '2'), mp.mpc(3, '0.01')]
    check("prop:sixtap the taps equal i e^{5iz} P(sin z), and 367 P(sin x) = 435 sin x + 80 sin 3x + 12 sin 5x",
          all(abs(Hs(z)-1j*mp.exp(5j*z)*P(mp.sin(z))) < 1e-25 for z in pts)
          and all(abs(367*P(mp.sin(x))-(435*mp.sin(x)+80*mp.sin(3*x)+12*mp.sin(5*x))) < 1e-25 for x in [mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('-2.2')]))
    check("prop:sixtap P(+-1) = +-1, P' = 15(8s^2-7)^2/367, and P is Pi_1 for s0 = sqrt3/2",
          P(1) == 1 and P(-1) == -1 and all(abs(mp.diff(P, s)-mp.mpf(15)/367*(8*s*s-7)**2) < 1e-20 for s in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.93')])
          and all(abs(ev(kernel(1)[1], s)/kernel(1)[2]-P(s)) < 1e-25 for s in [mp.mpf('0.2'), mp.mpf('0.7'), 1]))
    check("prop:sixtap |H| <= 1 on the real axis and in the upper half-plane (samples)",
          max(abs(Hs(x)) for x in mp.linspace(0, mp.pi, 200)) <= 1+mp.mpf(10)**-25 and all(abs(Hs(z)) <= 1 for z in pts))
    check("prop:sixtap sin 1 > 1 - 1/6 + 1/120 - 1/5040 > 5/6, sin 2 > sin 1, P(5/6) = 59225/59454",
          mp.sin(1) > 1-mp.mpf(1)/6+mp.mpf(1)/120-mp.mpf(1)/5040 > mp.mpf(5)/6 and mp.sin(2) > mp.sin(1)
          and Fr(192*5**5-560*5**3*36+735*5*6**4, 367*6**5) == Fr(59225, 59454))
    check("prop:sixtap P(sin x) in (59225/59454, 1] and arg H = pi/2 + 5x on [1,2]",
          all(mp.mpf(59225)/59454 < P(mp.sin(x)) <= 1 and abs(Hs(x)/abs(Hs(x))-mp.expj(mp.pi/2+5*x)) < 1e-25 for x in mp.linspace(1, 2, 51)))
    Ff = lambda t: t+mp.mpf(176)/1101*t**3+mp.mpf(192)/1835*t**5
    check("prop:sixtap d/dx(-F(cos x)) = P(sin x), int_1^2 P(sin x) = F(cos 1) - F(cos 2)",
          all(abs(-mp.diff(lambda y: Ff(mp.cos(y)), x)-P(mp.sin(x))) < 1e-20 for x in [mp.mpf('0.4'), mp.mpf('1.2'), mp.mpf('1.9')])
          and abs(mp.quad(lambda x: P(mp.sin(x)), [1, 2])-(Ff(mp.cos(1))-Ff(mp.cos(2)))) < 1e-25)
    fact = lambda n: 1 if n == 0 else n*fact(n-1)
    L = sum(Fr((-1)**k, fact(2*k)) for k in range(12)); U = sum(Fr((-1)**k*4**k, fact(2*k)) for k in range(13))
    Fq = lambda t: t+Fr(176, 1101)*t**3+Fr(192, 1835)*t**5
    cert = (1-Fq(L)+Fq(U))/Fr(59225, 59454)
    check("prop:sixtap L < cos 1, cos 2 < U, and (1-F(L)+F(U))/M = 0.00069587... < 7/10000 exactly",
          mp.mpf(L.numerator)/L.denominator < mp.cos(1) and mp.cos(2) < mp.mpf(U.numerator)/U.denominator
          and cert < Fr(7, 10000) and abs(mp.mpf(cert.numerator)/cert.denominator-mp.mpf('0.0006958765626')) < 1e-12)
    loss = mp.quad(lambda x: -mp.log(abs(Hs(x))), [1, 2])
    check("prop:sixtap mean loss 0.00069359 nepers, below the floor 0.0060483 by more than a factor eight",
          abs(loss-mp.mpf('0.00069359')) < 5e-9 and 8*loss < mp.mpf('0.0060483') and loss < mp.mpf(cert.numerator)/cert.denominator)
    roots = mp.polyroots([192, 0, -560, 0, 735])
    check("prop:sixtap the nonzero roots of P satisfy s^2 = (35 +- 14 i sqrt5)/24, and H vanishes at a point of C_+ above them",
          all(min(abs(r*r-(35+sg*14j*mp.sqrt(5))/24) for sg in [1, -1]) < 1e-20 for r in roots)
          and any(mp.im(mp.asin(r)) > 0 and abs(Hs(mp.asin(r))) < 1e-20 for r in roots))
    mp.mp.dps = dps0
construction_checks()

# ssec:channels
A0m = mp.matrix([[1, 0], [0, -1]]); C1m = mp.matrix([[1, 0], [0, 0]]); C2m = mp.matrix([[5, 2], [2, 1]])/4
check("thm:noloewner det C2 = det(C2-A0) = 1/16", abs(mp.det(C2m)-mp.mpf(1)/16) < 1e-18 and abs(mp.det(C2m-A0m)-mp.mpf(1)/16) < 1e-18)
ev = sorted(mp.eig(C2m-C1m)[0], key=lambda z: mp.re(z))
check("thm:noloewner eigenvalues of C2-C1 are -1/4 and 3/4", abs(ev[0]+mp.mpf(1)/4) < 1e-18 and abs(ev[1]-mp.mpf(3)/4) < 1e-18)
random.seed(11)
def herm(n):
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(i, n):
            x = mp.mpf(random.random()-0.5); yv = mp.mpf(random.random()-0.5) if i != j else 0
            M[i, j] = mp.mpc(x, yv); M[j, i] = mp.mpc(x, -yv)
    return M
def posp(M):
    E, Q = mp.eighe(M)
    return Q*mp.diag([max(e, 0) for e in E])*Q.H
pw = pw_setup(a, b); A = pw[4]
okt = True
for trial in range(6):
    n = 3; X = herm(n); Th = X*X.H; T = herm(n)*4
    B = t0*Th+T; Cp = posp(-B)
    Fo = lambda z: pw[5](z)*Th+z*T+Cp*Rc(z, a, b)/z
    for z in [mp.mpc(1.5, 0.3), mp.mpc(0.4, 2), mp.mpc(3, 0.1)]:
        F = Fo(z); ImF = (F-F.H)/(2j)
        okt = okt and min(mp.re(e) for e in mp.eighe(ImF)[0]) > -1e-15
    x = mp.mpf('1.4'); F = Fo(mp.mpc(x, mp.mpf('1e-14'))); ImF = (F-F.H)/(2j)
    trb = A(x)*sum(Th[i, i] for i in range(n))+r(x)/x*sum(Cp[i, i] for i in range(n))
    okt = okt and abs(sum(ImF[i, i] for i in range(n))-trb) < 1e-9
    # every orthonormal basis gives a trace bound at most the eigenbasis one
    Ub = mp.eighe(herm(n))[1]
    alt = A(x)*sum(Th[i, i] for i in range(n))+r(x)/x*sum(max(0, -mp.re((Ub[:, j].H*B*Ub[:, j])[0])) for j in range(n))
    okt = okt and mp.re(alt) <= mp.re(trb)+1e-15
check("thm:trace optimum is a generator and attains the trace bound (random)", okt)

# ssec:bands, two bands
bands = [(mp.mpf(1), mp.mpf(2)), (mp.mpf(3), mp.mpf(4))]; cg = mp.mpf('2.5')
def Rm(z):
    z = mp.mpc(z); out = mp.mpc(1)
    for (aa, bb) in bands: out *= z*z*mp.sqrt(1-aa*aa/(z*z))*mp.sqrt(1-bb*bb/(z*z))
    return out
Dm = lambda z: Rm(z)/(z*z-cg*cg); Qm = lambda z: Dm(z)/z
pts = [mp.mpc(x, y) for x in [-5, mp.mpf('-2.2'), mp.mpf('0.3'), mp.mpf('1.5'), mp.mpf('2.5'), mp.mpf('3.5'), 6] for y in [mp.mpf('0.01'), mp.mpf('0.3'), 2, 10]]
check("lem:multicomp Q_m Herglotz on two bands", min(mp.im(Qm(z)) for z in pts) > 0)
check("lem:multicomp D_m = i rho, rho > 0, on both bands", all(abs(mp.re(Dm(mp.mpc(x, mp.mpf('1e-15'))))) < 1e-10 and mp.im(Dm(mp.mpc(x, mp.mpf('1e-15')))) > 0 for x in [mp.mpf('1.5'), mp.mpf('3.5')]))
check("lem:multicomp nu_0 = a1 b1 a2 b2 / c^2", abs(-Dm(mp.mpc(0, mp.mpf('1e-18'))) - 24/cg**2) < 1e-12)
hh = mp.mpf('1e-8')
check("lem:multicomp gap residue negative, linear coefficient 1", mp.re(Qm(cg+hh))*hh < 0 and abs(Qm(mp.mpc(0, 10**6))/mp.mpc(0, 10**6)-1) < 1e-10)

# ssec:severalexact, the exact price on several bands
def sev_setup(ends):
    a = [mp.mpf(v) for v in ends[0::2]]; b = [mp.mpf(v) for v in ends[1::2]]; m = len(a)
    A = [x*x for x in a]; B = [x*x for x in b]
    X = lambda s: mp.fprod([s-u for u in A]); Y = lambda s: mp.fprod([s-u for u in B])
    g = [mp.findroot(lambda s: X(s)-Y(s), (B[j], A[j+1]), solver='bisect') for j in range(m-1)]
    P = lambda s: mp.fprod([s-x for x in g])
    def br(F):
        tot = 0
        for j in range(m):
            oth = [k for k in range(m) if k != j]
            def hfun(ph, j=j, oth=oth):
                s = A[j]+(B[j]-A[j])*(1-mp.cos(ph))/2
                o = mp.fprod([abs((s-A[k])*(s-B[k])) for k in oth])
                return F(mp.sqrt(s))/(2*mp.sqrt(s)*mp.sqrt(o))
            tot += mp.quad(hfun, [0, mp.pi])
        return tot
    def rhoB(z):
        out = mp.mpc(1)
        for aa, bb in zip(a, b):
            out *= mp.sqrt((z-aa)/(z-bb))*mp.sqrt((z+aa)/(z+bb))
        return out
    WB = lambda z: (rhoB(z)-1)**2/rhoB(z)
    return a, b, m, A, B, X, Y, g, P, br, WB
sev_ok = {k: True for k in ["roots", "W", "band", "mom", "glauert", "equal", "alpha", "super"]}
for ends in [(1, 2, 3, 4), (1, '1.5', 4, 5), (1, 2, 3, 4, 5, 6), (1, '1.2', 2, 3, '3.1', 7)]:
    a, b, m, A, B, X, Y, g, P, br, WB = sev_setup(ends)
    c = [mp.sqrt(x) for x in g]
    sev_ok["roots"] &= all(B[j] < g[j] < A[j+1] for j in range(m-1))
    offs = [mp.mpf('0.3')*a[0], (b[0]+a[1])/2] + [b[-1]+1, 3*b[-1]] + [(b[j]+c[j])/2 for j in range(m-1)]
    sev_ok["W"] &= all(mp.re(WB(mp.mpc(t, mp.mpf('1e-30')))) > 0 for t in offs) and all(abs(WB(mp.mpc(cc, mp.mpf('1e-30')))) < 1e-15 for cc in c)
    pB = lambda t: (abs(Y(t*t))-abs(X(t*t)))
    for j in range(m):
        x = (a[j]+2*b[j])/3; w = WB(mp.mpc(x, mp.mpf('1e-30')))
        sev_ok["band"] &= abs(mp.re(w)+2) < 1e-12 and abs(mp.im(w)-pB(x)/mp.sqrt(abs(X(x*x)*Y(x*x)))) < 1e-10
    d = sum(B)-sum(A)
    sev_ok["mom"] &= abs(br(lambda t: t*pB(t))) < 1e-15 and abs(br(lambda t: t**3*pB(t))+mp.pi*d*d/8) < 1e-12
    W0 = WB(mp.mpc(mp.mpf('1e-30'), mp.mpf('1e-30')))
    sev_ok["mom"] &= abs(br(lambda t: pB(t)/t)-mp.pi/2*mp.re(W0)) < 1e-12
    ip = lambda t: abs(P(t*t))
    sev_ok["glauert"] &= all(abs(br(lambda t: ip(t)/(cc-t))-br(lambda t: ip(t)/(cc+t))) < 1e-14 for cc in c)
    rho = lambda x: mp.sqrt(abs(X(x*x)*Y(x*x)))/abs(P(x*x))
    alpha = lambda x: 2*rho(x)/mp.pi*br(lambda t: ip(t)/(x+t))
    loss = sum(mp.quad(alpha, [a[j], (a[j]+b[j])/2, b[j]]) for j in range(m))
    JB = br(pB)
    sev_ok["equal"] &= abs(loss-JB/2) < 1e-12
    sev_ok["alpha"] &= all(alpha(a[j]+(b[j]-a[j])*k/7) > 0 for j in range(m) for k in range(1, 7))
    Js = sum(sev_setup((a[j], b[j]))[9](lambda t, j=j: (b[j]**2-t*t)-(t*t-a[j]**2)) for j in range(m))
    sev_ok["super"] &= JB > Js
check("lem:WB(iii) one root of X-Y in each gap", sev_ok["roots"])
check("lem:WB(iii) W_B > 0 off the bands, = 0 at c_j", sev_ok["W"])
check("lem:WB(iv) W_B = -2 + i p_B on the bands", sev_ok["band"])
check("lem:kernelB moments: first 0, third -pi delta^2/8, p_B/x = pi W_B(0)/2", sev_ok["mom"])
check("lem:glauertB identity at the gap points", sev_ok["glauert"])
check("thm:severalbands(iii) loss of H_B equals theta J_B/2", sev_ok["equal"])
check("thm:severalbands(i) alpha > 0 on the bands", sev_ok["alpha"])
check("thm:severalbands(iv) J_B > sum of single-band J", sev_ok["super"])
a, b, m, A, B, X, Y, g, P, br, WB = sev_setup((1, 2, 3, 4))
JB = br(lambda t: (abs(Y(t*t))-abs(X(t*t))))
check("two bands: c_1^2 = 11/2 on [1,2]u[3,4]", abs(g[0]-mp.mpf(11)/2) < 1e-20)
check("two bands: Lambda_B = 0.36445, total 0.72890", abs(JB/4-mp.mpf('0.36445')) < 5e-6 and abs(JB/2-mp.mpf('0.72890')) < 5e-6)
L34 = (K(mp.mpf(1)/7)-E(mp.mpf(1)/7))*7; L12 = (K(mp.mpf(1)/3)-E(mp.mpf(1)/3))*3
check("single bands: 0.27353 + 0.11307 = 0.38660", abs(L12-mp.mpf('0.27353')) < 5e-6 and abs(L34-mp.mpf('0.11307')) < 5e-6 and abs(L12+L34-mp.mpf('0.38660')) < 5e-6)
lams = []
for mm in range(1, 5):
    ee = []
    for j in range(mm): ee += [2*j+1, 2*j+2]
    _a, _b, _m, _A, _B, _X, _Y, _g, _P, _br, _W = sev_setup(ee)
    lams.append(_br(lambda t: (abs(_Y(t*t))-abs(_X(t*t))))/(2*mm))
check("fig:bandloss constants 0.27353, 0.36445, 0.41420, 0.44683", all(abs(l-mp.mpf(v)) < 5e-6 for l, v in zip(lams, ['0.27353', '0.36445', '0.41420', '0.44683'])))

print("%d checks passed, %d failed" % (passed, failed))
