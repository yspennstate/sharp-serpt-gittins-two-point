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
