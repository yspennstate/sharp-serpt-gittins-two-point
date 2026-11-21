import sympy as sp
p, t, z, rho = sp.symbols('p t z rho', positive=True)
H = t**5 + 4*t**4 - 20*t**3 + 12*t**2 + 16*t - 8
Hp = 16*(1-t)**4 + 88*t*(1-t)**3 + 108*t**2*(1-t)**2 + 32*t**3*(1-t) + t**4
print("H' identity:", sp.expand(sp.diff(H, t) - Hp) == 0)
K = (t**2-2)*p**2 + 2*t*(1-t)*p + t
L = (p**3-2*p**2-1)*t**2 + 4*p*t - 2*p**2
N = t*(1-p)*(p-t); D = (1-t)*(2*p-(1-p+p**2)*t)
print("Np:", sp.expand(sp.diff(N,p)*D - N*sp.diff(D,p) + t*(t-1)*K) == 0)
print("Nt:", sp.expand(sp.diff(N,t)*D - N*sp.diff(D,t) - (p-1)*L) == 0)
print("K(t,t):", sp.factor(K.subs(p, t)), " K(1,t):", sp.factor(K.subs(p, 1)))
print("Res:", sp.factor(sp.resultant(K, L, p)))
# f from FCFS boundary
lam = t; q = 1-p
EX = 1/p; EX2 = (1+p-p**2)/p**2
fcfs = EX + lam*EX2/(2*(1-lam*EX))
gap = lam*q/(2*p*(1-lam))
print("f:", sp.simplify(gap/fcfs - N/D) == 0)
# SERPT: conditional means, a=1
def ET(w, p, lam, b):
    q = 1-p; A = 1-lam; EX = p + q*b; EX2 = p + q*b**2; r = lam*EX
    B = lam*EX2/(2*(1-r)*A)
    return p*(1 + lam*(1 + q*w**2)/(2*A)) + q*(B + (b-w)/A + w)
q = 1-p
b = 1 + 1/p + z; mu = 1/p + q*z; lamr = rho/mu
G = ET(1/p, p, lamr, b); S = ET(mu, p, lamr, b)
Dd = 2*(1-lamr)*(1-rho)*mu*G
J = p + q/p
Dproof = (1-rho)*((2*mu-rho)*J + 2*mu*q*(1+z)) + rho*q*(1 + 2*q*(1/p+z) + q*(1/p+z)**2)
print("D:", sp.simplify(Dd - Dproof) == 0)
print("ratio-1:", sp.simplify((S-G)/G - rho*p*q**3*(1-rho)*z**2/Dproof) == 0)
P = sp.Poly(sp.expand(Dproof*p**2), z)
print("D coeffs*p^2:", [sp.factor(c) for c in P.all_coeffs()])
# Gittins mean at the boundary d = 1/p (a = 1, b = 1 + 1/p) equals the FCFS (Pollaczek-Khinchine) mean
print("boundary:", sp.simplify(ET(1/p, p, t, 1 + 1/p) - fcfs) == 0)
# FCFS minus Gittins for pd > a (general a, d, p, lambda): the remark after Proposition 3.3
a, d, lm = sp.symbols('a d lm', positive=True)
q = 1-p; b = a + d; mu = a + q*d; EX2 = p*a**2 + q*b**2; r = lm*mu; A = 1 - lm*a
B = lm*EX2/(2*(1-r)*A)
Gm = p*(a + lm*(a**2 + q*(a/p)**2)/(2*A)) + q*(B + (b-a/p)/A + a/p)
Ff = mu + lm*EX2/(2*(1-r))
Y = p*b - q*a + lm*a*(q*b - p*d)
print("FCFS-Gittins:", sp.simplify(Ff - Gm - lm*q*(p*d-a)*Y/(2*p*A*(1-r))) == 0)
print("mu*Y bound:", sp.expand(mu*(p*b-q*a) + a*(q*b-p*d) - p*(a**2 + 2*q*a*d + q*d**2)) == 0)
