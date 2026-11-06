"""Independent EOS evaluation for verification only. No stellar ODE here.
The unmodified R17 TOV/GL4 solvers consume these exact native Hermite cells.
SI source -> geometric cm units: P_geom=G_CGS*P_SI*10/C_CGS**4.
"""
from pathlib import Path
import numpy as np
G=6.67430e-11;C=299792458.; CONV=G/C**4/1e4
class NativeEOS:
 def __init__(self,path,alpha=0.,split=69):
  a=np.genfromtxt(path,delimiter=',',names=True)
  if a.ndim!=1 or len(a)<3 or not {'mu_b_j','pressure_pa','n_b_m3'}.issubset(a.dtype.names or ()):
   raise ValueError('Expected normalized SI table with at least three rows')
  if type(split) is not int or not 2<=split<=len(a) or not np.isfinite(alpha):raise ValueError('Invalid perturbation control')
  self.mu=np.array(a['mu_b_j']);self.p=np.array(a['pressure_pa']);self.n=np.array(a['n_b_m3'])
  if not all(np.all(np.isfinite(x)) and np.all(x>0) for x in [self.mu,self.p,self.n]):
   raise ValueError('Nonfinite or nonpositive EOS input')
  if not np.all(np.diff(self.mu)>0):raise ValueError('Material jumps or unordered mu are not supported')
  j=split-1
  if alpha:
   mu=self.mu[:j];x=(mu-self.mu[0])/(self.mu[j]-self.mu[0]);b=x*x*(1-x)**2
   db=2*x*(1-x)*(1-2*x)/(self.mu[j]-self.mu[0]);p=self.p[:j].copy();n=self.n[:j].copy()
   self.p[:j]=p*(1+alpha*b);self.n[:j]=n*(1+alpha*b)+alpha*p*db
  self.h=np.diff(self.mu);n0=self.n[:-1];n1=self.n[1:];p0=self.p[:-1];p1=self.p[1:]
  self.coeff=np.array([p0,self.h*n0,3*(p1-p0)-self.h*(2*n0+n1),2*(p0-p1)+self.h*(n0+n1)]).T
  if not np.all(self.h>0) or not np.all(np.diff(self.p)>0):raise ValueError('unordered EOS')
  # On a cubic P(mu), P'' is linear and d(mu P''-P')/dmu=mu P''' has fixed sign.
  dn0=2*self.coeff[:,2]/self.h**2;dn1=(2*self.coeff[:,2]+6*self.coeff[:,3])/self.h**2
  if np.any(dn0<=0) or np.any(dn1<=0) or np.any(self.mu[:-1]*dn0<self.n[:-1]) or np.any(self.mu[1:]*dn1<self.n[1:]):
   raise ValueError('Native cubic is unstable or acausal; no smoothing is applied')
 def at_t(self,i,t):
  if not 0<=i<len(self.h):raise ValueError('Unknown native cell')
  t=np.asarray(t)
  if not np.all(np.isfinite(t)) or np.any(t < -128*np.finfo(float).eps) or np.any(t>1+128*np.finfo(float).eps):
   raise ValueError('Local coordinate outside the cell')
  t=np.clip(t,0.,1.)
  a,b,c,d=self.coeff[i];h=self.h[i];mu=self.mu[i]+h*t
  p=a+t*(b+t*(c+d*t));n=(b+t*(2*c+3*d*t))/h;dn=(2*c+6*d*t)/h**2
  e=mu*n-p;chi=n/(mu*dn)
  if not all(np.all(np.isfinite(x)) and np.all(x>0) for x in [p,e,n,dn,chi]) or np.any(chi>1):raise ValueError('noncausal EOS')
  return p,e,chi,n,mu
 def inverse(self,i,p):
  p=np.asarray(p);low,high=self.p[i:i+2]
  tol=128*np.finfo(float).eps*max(low,high)
  if np.any(p<low-tol) or np.any(p>high+tol):raise ValueError('EOS domain')
  p=np.clip(p,low,high) # endpoint roundoff only; prior inequality excludes extrapolation
  a,b,c,d=self.coeff[i];delta=p-a;t=np.clip(delta/(high-low),0,1);lo=np.zeros_like(t);hi=np.ones_like(t)
  for _ in range(100):
   f=t*(b+t*(c+d*t))-delta
   good=abs(f)<=8*np.finfo(float).eps*np.maximum(abs(p),low)
   if np.all(good):return t
   lo=np.where(f<0,t,lo);hi=np.where(f>0,t,hi)
   tn=t-f/(b+t*(2*c+3*d*t));bad=(tn<=lo)|(tn>=hi)|~np.isfinite(tn)
   tn=np.where(bad,.5*(lo+hi),tn);t=np.where(good,t,tn)
  raise ArithmeticError('EOS inversion did not close')
 def parts(self,pc,cut):
  limits=np.r_[pc,self.p[(self.p>cut)&(self.p<pc)][::-1],cut]
  if not self.p[0]<=cut<pc<=self.p[-1]:raise ValueError('outside table')
  parts=[]
  for top,bot in zip(limits[:-1],limits[1:]):
   i=int(np.searchsorted(self.p,(top+bot)*.5)-1)
   def eos(pgeom,index=i):
    t=self.inverse(index,np.asarray(pgeom)/CONV);p,e,chi,n,mu=self.at_t(index,t)
    return e*CONV,chi
   parts.append(dict(top=float(top*CONV),bottom=float(bot*CONV),eos=eos,side=2,region='core',native_cell=i))
  return parts
