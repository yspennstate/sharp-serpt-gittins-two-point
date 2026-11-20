#!/usr/bin/env python3
"""Serial checks for the binary M-SERPT calculation. No network or external writes."""
import argparse
import json
import math
import random
from fractions import Fraction as F


def means(a, b, p, lam):
    q = 1 - p
    mu = p * a + q * b
    rho = lam * mu
    assert a > 0 and b >= a and 0 <= p <= 1 and 0 <= rho < 1
    if p in (0, 1) or b == a or p * (b - a) <= a:
        fcfs = mu + lam * (p * a*a + q * b*b) / (2 * (1-rho))
        return fcfs, fcfs
    A, h = 1-lam*a, a/p
    B = lam * (p*a*a+q*b*b) / (2*(1-rho)*A)
    ms = p*(a+lam*a*a/(2*A)) + q*(B+b/A)
    gi = p*(a+lam*(a*a+q*h*h)/(2*A)) + q*(B+(b-h)/A+h)
    return ms, gi


def limiting_gap(p, t):
    return t*(1-p)*(p-t)/((1-t)*(2*p-(1-p+p*p)*t))


def cutoff_mean(a,b,p,lam,h):
    """High initial stage, low middle stage, protected final h units."""
    q=1-p;mu=p*a+q*b;A=1-lam*a
    assert 0<=h<b-a and lam*mu<1
    B=lam*(p*a*a+q*b*b)/(2*(1-lam*mu)*A)
    return p*(a+lam*(a*a+q*h*h)/(2*A))+q*(B+(b-h)/A+h)


def serpt_mean(a,b,p,lam):
    mu=p*a+(1-p)*b
    if p in (0,1) or p*(b-a)<=a:
        return means(a,b,p,lam)[1]
    return cutoff_mean(a,b,p,lam,mu)


def H(t):
    return t**5+4*t**4-20*t**3+12*t*t+16*t-8


def K(p, t):
    return (t*t-2)*p*p+2*t*(1-t)*p+t


class Interval:
    """Exact rational interval arithmetic, used only to certify decimal bounds."""
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def coerce(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, x):
        x = self.coerce(x)
        return Interval(self.lo+x.lo, self.hi+x.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, x):
        return self + -self.coerce(x)

    def __rsub__(self, x):
        return self.coerce(x) + -self

    def __mul__(self, x):
        x = self.coerce(x)
        z = [a*b for a in (self.lo,self.hi) for b in (x.lo,x.hi)]
        return Interval(min(z), max(z))
    __rmul__ = __mul__

    def __truediv__(self, x):
        x = self.coerce(x)
        assert not x.lo <= 0 <= x.hi
        return self * Interval(1/x.hi, 1/x.lo)


def exact_checks():
    import sympy as s
    import mpmath as mp
    p,t,lam,d = s.symbols('p t lam d')
    q=1-p; b=1+d; A=1-lam; rho=lam*(1+q*d); h=1/p
    B=lam*(p+q*b*b)/(2*(1-rho)*A)
    M=p*(1+lam/(2*A))+q*(B+b/A)
    G=p*(1+lam*(1+q*h*h)/(2*A))+q*(B+(b-h)/A+h)
    assert s.cancel(M-G-lam*q/(2*p*A)) == 0
    fcfs = 1/p + lam*(1+p-p*p)/(2*p*(p-lam))
    assert s.cancel(G.subs(d,1/p)-fcfs) == 0
    f=limiting_gap(p,t)
    assert s.cancel((M-G).subs({d:1/p,lam:t})/G.subs({d:1/p,lam:t})-f)==0
    kval = K(p,t)
    lval=(p**3-2*p*p-1)*t*t+4*p*t-2*p*p
    N=t*(1-p)*(p-t);D=(1-t)*(2*p-(1-p+p*p)*t)
    assert s.expand(s.diff(N,p)*D-N*s.diff(D,p)+t*(t-1)*kval)==0
    assert s.expand(s.diff(N,t)*D-N*s.diff(D,t)-(p-1)*lval)==0
    resultant=s.factor(s.resultant(kval,lval,p))
    assert s.expand(resultant-t*t*(t-1)*(t*t-t-1)*H(t))==0
    bernstein=16*(1-t)**4+88*t*(1-t)**3+108*t*t*(1-t)**2+32*t**3*(1-t)+t**4
    assert s.expand(s.diff(H(t),t)-bernstein)==0
    print('symbolic identities: PASS; resultant =',resultant)
    hh=s.symbols('hh')
    Th=p*(1+lam*(1+q*hh*hh)/(2*A))+q*(B+(b-hh)/A+hh)
    assert s.cancel(Th-G-lam*p*q*(hh-1/p)**2/(2*A))==0
    r,z=s.symbols('r z',positive=True)
    muz=1/p+q*z;J=p+q/p;EX2=1+2*q*(1/p+z)+q*(1/p+z)**2
    Dz=(1-r)*((2*muz-r)*J+2*muz*q*(1+z))+r*q*EX2
    assert s.cancel((Dz-2*(1-lam)*(1-rho)*(1+q*d)*G).subs({d:1/p+z,lam:r/muz}))==0
    assert s.simplify(s.Poly(Dz,z).coeff_monomial(z*z)-q*q*(2-r))==0
    rr=(r*(1-r)/(2-r))
    assert s.factor(s.together(s.diff(rr,r)))==(r*r-4*r+2)/(r-2)**2
    print('cutoff quadratic regret, SERPT denominator and derivative identities: PASS')

    # Exact signs and interval enclosures; no floating-point proof obligation.
    tlo=F('0.450614580799471903808270806938')
    thi=F('0.450614580799471903808270806940')
    plo=F('0.6571391408919625221961487')
    phi=F('0.6571391408919625221961489')
    assert H(tlo)<0<H(thi)
    ti=Interval(tlo,thi);pi=Interval(plo,phi)
    assert K(Interval(plo),ti).lo>0 and K(Interval(phi),ti).hi<0
    fi=limiting_gap(pi,ti)+1
    clo,chi=F('1.060173421912767960604515'),F('1.060173421912767960604516')
    assert clo<fi.lo<=fi.hi<chi
    print('certified C interval:',str(clo),str(chi))
    mp.mp.dps=45
    ts=mp.findroot(H,.45)
    ps=(ts*(1-ts)+mp.sqrt(ts*ts*(1-ts)**2+ts*(2-ts*ts)))/(2-ts*ts)
    print('stationary point:',json.dumps({k:str(v) for k,v in {
        't':ts,'p':ps,'rho':ts/ps,'d_over_a':1/ps,'C':1+limiting_gap(ps,ts)}.items()}))
    checks=0; max_grid=(0,None)
    for i in range(1,30):
        pp=F(i,30)
        for j in range(1,30):
            rr=F(j,30)
            for factor in (F(0),F(1,2),F(1),F(1001,1000),F(2),F(10)):
                aa=F(1);bb=aa+factor/pp;ll=rr/(pp*aa+(1-pp)*bb)
                mm,gg=means(aa,bb,pp,ll)
                ss=serpt_mean(aa,bb,pp,ll)
                assert gg<=mm and mm/gg<chi
                assert gg<=ss and ss/gg<F('1.042893218813453')
                if factor>1:
                    assert mm-gg==ll*(1-pp)/(2*pp*(1-ll))
                    assert mm/gg<=1+limiting_gap(pp,ll)
                    mu=pp*aa+(1-pp)*bb
                    assert ss-gg==ll*pp*(1-pp)*(mu-aa/pp)**2/(2*(1-ll*aa))
                    bound=1+rr*(1-rr)/(2-rr)*pp*(1-pp)
                    assert ss/gg<bound
                else:
                    assert mm==gg
                    assert ss==gg
                if mm/gg>max_grid[0]: max_grid=(mm/gg,(pp,rr,factor))
                checks+=1
    # Degenerate, zero-load, scale, published SOAP example, and unstable inputs.
    for pp in (F(0),F(1),F(1,2)):
        assert means(F(1),F(3),pp,F(0))[0]==3-2*pp
        if pp in (0,1):
            assert means(F(1),F(3),pp,F(1,10))[0]==means(F(1),F(3),pp,F(1,10))[1]
    mm,gg=means(F(2),F(14),F(1,2),F(1,10))
    published=(F(6,10)/(1-F(2,10))+2+F(5)/((1-F(8,10))*(1-F(2,10)))+10/(1-F(2,10))+4)/2
    assert gg==published
    assert means(F(6),F(42),F(1,2),F(1,30))==(3*mm,3*gg)
    try:
        means(F(1),F(3),F(1,2),F(1,2))
        raise RuntimeError('unstable input accepted')
    except AssertionError:
        pass
    print('exact rational cases:',checks,'grid max:',float(max_grid[0]),'parameters:',list(map(str,max_grid[1])))
    print('SOAP CoinFlip{2,14}, lambda=1/10: M =',mm,'G =',gg,'ratio =',mm/gg)
    print('degenerate, zero-load, scale, stability checks: PASS')
    for eps in (F(1,10),F(1,1000),F(1,10**6),F(1,10**9)):
        pp=F(str(ps));tt=F(str(ts));bb=1+(1+eps)/pp
        mm,gg=means(F(1),bb,pp,tt)
        print('sharpness epsilon:',str(eps),'ratio:',format(float(mm/gg),'.15f'))
    pp=F(2,3);tt=F(9,20);bb=1+1/pp
    print('boundary p=2/3, lambda=9/20: actual ratio =',means(F(1),bb,pp,tt)[0]/means(F(1),bb,pp,tt)[1],
          'right limit =',float(1+limiting_gap(pp,tt)))
    print('SERPT sharp constant:',mp.nstr((7-2*mp.sqrt(2))/4,45))
    for dd in (10,1000,10**6,10**9):
        pp=F(1,2);bb=F(1+dd);rr=F(str(2-mp.sqrt(2)));ll=rr/(pp+(1-pp)*bb)
        ss=serpt_mean(F(1),bb,pp,ll);gg=means(F(1),bb,pp,ll)[1]
        print('SERPT sharpness d:',dd,'ratio:',format(float(ss/gg),'.15f'))


def simulate(a,b,p,lam,policy,seed,count=100000,warmup=2000):
    """Event-driven sample paths from ranks, without using response-time formulas."""
    rng=random.Random(seed)
    time=0.;next_arrival=rng.expovariate(lam);arrived=0;queue=[]
    total=0.;short=0.;long=0.;ns=0;nl=0
    mu=p*a+(1-p)*b
    def rank(job):
        age,stage=job[2],job[3]
        if policy=='M': return mu if stage==0 else max(mu,b-a)
        if policy=='S': return mu-age if stage==0 else b-age
        return min((a-age)/p,mu-age) if stage==0 else b-age
    while arrived<count or queue:
        if not queue:
            time=next_arrival
            large=rng.random()>=p
            queue.append([arrived,time,0.,0,b if large else a])
            arrived+=1
            next_arrival=time+rng.expovariate(lam) if arrived<count else math.inf
        index=min(range(len(queue)),key=lambda i:(rank(queue[i]),queue[i][0]))
        job=queue[index]
        checkpoint=a if job[3]==0 else b
        end=time+checkpoint-job[2]
        if next_arrival<end:
            job[2]+=next_arrival-time
            time=next_arrival
            large=rng.random()>=p
            queue.append([arrived,time,0.,0,b if large else a])
            arrived+=1
            next_arrival=time+rng.expovariate(lam) if arrived<count else math.inf
        else:
            time=end;job[2]=checkpoint
            if job[4]==checkpoint:
                if job[0]>=warmup:
                    sojourn=time-job[1];total+=sojourn
                    if job[4]==a:short+=sojourn;ns+=1
                    else:long+=sojourn;nl+=1
                queue.pop(index)
            else: job[3]=1
    return {'mean':total/(ns+nl),'short':short/ns,'long':long/nl,'jobs':ns+nl}


def simulation_checks():
    cases=[('near-bound',1.,1.+1.000001/.6571391408919625,.6571391408919625,.4506145807994719),
           ('SOAP-example',2.,14.,.5,.1),('below-transition',1.,2.,.5,.4),
           ('at-transition',1.,2.5,2/3,.45),
           ('SERPT-near-bound',1.,1001.,.5,(2-math.sqrt(2))/501)]
    for name,a,b,p,lam in cases:
        ms,gi=means(a,b,p,lam)
        observed={policy:simulate(a,b,p,lam,policy,7115) for policy in ('M','G','S')}
        print(json.dumps({'case':name,'parameters':[a,b,p,lam],
            'theory':{'M':ms,'G':gi,'S':serpt_mean(a,b,p,lam),'ratio':ms/gi},'simulation':observed,
            'observed_ratio':observed['M']['mean']/observed['G']['mean']}))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['exact','simulate','all'],default='all',nargs='?')
    mode=ap.parse_args().mode
    if mode in ('exact','all'):exact_checks()
    if mode in ('simulate','all'):simulation_checks()
