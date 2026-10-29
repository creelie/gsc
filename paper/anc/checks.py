"""Numerical cross-checks for the paper 'The Energy Cost of a Phase Shift'.

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
print("%d checks passed, %d failed" % (passed, failed))
