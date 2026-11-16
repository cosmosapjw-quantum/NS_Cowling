"""Source-defined nucleonic RMF field inputs for the eight GPPVA crusts.

Raw fields phi, omega, rho and masses are in MeV; positions in fm and
number/scalar densities in fm^-3. The rho vertex is g_rho*tau3/2.
M_i*=M_i-g_sigma*phi. Potentials U and W have units MeV^4.
The density-dependent couplings use LOCAL baryon density in a static cell.
This module supplies the field functional and uniform stationary references;
finite-range geometry and profiles are solved in spatial_tf.py.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq
from .cold_leptons import HBARC, ELECTRON_MASS, free_fermion, leptons
from .gm1_uniform import scalar_density

MODELS = ('DD2','DDME2','TW','FSU2','FSU2R','FSU2H','NL3wrL55','TM1e')
_DD = {
 'DD2': dict(mn=939.56536,mp=938.27203,n0=.149065,ms=546.212459,mw=783.,mr=763.,gs=10.686681,gw=13.342362,gr=2*3.626940,sa=1.357630,sb=.634442,sc=1.005358,sd=.575810,wa=1.369718,wb=.496475,wc=.817753,wd=.638452,ra=.518903),
 'DDME2': dict(mn=939.,mp=939.,n0=.152,ms=550.1238,mw=783.,mr=763.,gs=10.5396,gw=13.0189,gr=2*3.6836,sa=1.3881,sb=1.0943,sc=1.7057,sd=.4421,wa=1.3892,wb=.9240,wc=1.4620,wd=.4775,ra=.5647),
 'TW': dict(mn=939.,mp=939.,n0=.153,ms=550.,mw=783.,mr=763.,gs=10.72854,gw=13.29015,gr=2*3.66098,sa=1.365469,sb=.226061,sc=.409704,sd=.901995,wa=1.402488,wb=.172577,wc=.344293,wd=.983955,ra=.515),
}
_NL = {
 'FSU2': dict(mn=939.,mp=939.,ms=497.479,mw=782.5,mr=763.,gs=math.sqrt(108.0943),gw=math.sqrt(183.7893),gr=math.sqrt(80.4656),kappa=3.0029,lam=-.000533,zeta=.0256,Lambda=.000823),
 'FSU2R': dict(mn=939.,mp=939.,ms=497.479,mw=782.5,mr=763.,gs=math.sqrt(107.5751),gw=math.sqrt(182.3949),gr=math.sqrt(206.4260),kappa=3.0911,lam=-.001680,zeta=.024,Lambda=.045),
 'FSU2H': dict(mn=939.,mp=939.,ms=497.479,mw=782.5,mr=763.,gs=math.sqrt(102.7200),gw=math.sqrt(169.5315),gr=math.sqrt(197.2692),kappa=4.0014,lam=-.013298,zeta=.008,Lambda=.045),
 'NL3wrL55': dict(mn=939.,mp=939.,ms=508.194,mw=782.501,mr=763.,gs=10.217,gw=12.868,gr=11.2766,kappa_unscaled=4.384*939.,lambda_unscaled=-173.31,Lambda=.03),
 'TM1e': dict(mn=938.,mp=938.,ms=511.198,mw=783.,mr=770.,gs=10.0289,gw=12.6139,gr=13.9714,g2=-7.2325,g3=.6183,c3=71.3075,Lambda=.0429),
}
SOURCES = {
 'DD2':'https://doi.org/10.1103/PhysRevC.81.015803',
 'DDME2':'https://doi.org/10.1103/PhysRevC.71.024312',
 'TW':'https://doi.org/10.1016/S0375-9474(99)00310-3',
 'FSU2':'https://iopscience.iop.org/0004-637X/863/1/104/suppdata/apjaad049t1_ascii.txt',
 'FSU2R':'https://arxiv.org/abs/1708.08681',
 'FSU2H':'https://arxiv.org/abs/1708.08681',
 'NL3wrL55':'https://arxiv.org/html/1203.3052',
 'TM1e':'https://user.numazu-ct.ac.jp/~sumi/eos/table4/guide_EOS4.pdf',
}

def model(name: str) -> dict:
    name=name.removeprefix('GPPVA_')
    if name not in MODELS:raise ValueError('unknown GPPVA model: '+name)
    p=dict((_DD if name in _DD else _NL)[name]);p.update(name=name,kind='density_dependent' if name in _DD else 'nonlinear',source=SOURCES[name],rho_vertex='g_rho*tau3/2')
    # U=.5ms²phi²+a3 phi³/3+a4 phi⁴/4;
    # W=.5mw²w²+.5mr²r²+b4 w⁴/4+cross w²r².
    if name in _DD:p.update(a3=0.,a4=0.,b4=0.,cross=0.)
    elif name.startswith('FSU'):
        p.update(a3=p['kappa']*p['gs']**3/2,a4=p['lam']*p['gs']**4/6,b4=p['zeta']*p['gw']**4/6,cross=p['Lambda']*p['gw']**2*p['gr']**2)
    elif name=='NL3wrL55':
        p.update(a3=p['kappa_unscaled']/2,a4=p['lambda_unscaled']/6,b4=0.,cross=p['Lambda']*p['gw']**2*p['gr']**2)
    else:
        p.update(a3=-p['g2']*HBARC,a4=p['g3'],b4=p['c3'],cross=p['Lambda']*p['gw']**2*p['gr']**2)
    return p

def couplings(name: str, density) -> dict:
    """Vectorized vertices and derivatives dg/dn in fm^3, including n=0."""
    p=model(name);n=np.asarray(density,dtype=float)
    if not np.all(np.isfinite(n)) or np.any(n<0):raise ValueError('finite nonnegative local density required')
    if p['kind']=='nonlinear':return {**{g:np.full_like(n,p[g]) for g in ('gs','gw','gr')},**{g:np.zeros_like(n) for g in ('dgs','dgw','dgr')}}
    x=n/p['n0'];result={}
    for vertex,prefix in [('gs','s'),('gw','w')]:
        a,b,c,d=(p[prefix+k] for k in 'abcd');u=x+d
        result[vertex]=p[vertex]*a*(1+b*u*u)/(1+c*u*u)
        result['d'+vertex]=p[vertex]*2*a*(b-c)*u/(1+c*u*u)**2/p['n0']
    result['gr']=p['gr']*np.exp(-p['ra']*(x-1))
    result['dgr']=-p['ra']*result['gr']/p['n0']
    return result

def scalar_potential(name: str, phi) -> dict:
    p=model(name);s=np.asarray(phi,dtype=float)
    if not np.all(np.isfinite(s)):raise ValueError('finite scalar field required')
    return {'U':p['ms']**2*s*s/2+p['a3']*s**3/3+p['a4']*s**4/4,
            'dU':p['ms']**2*s+p['a3']*s*s+p['a4']*s**3,
            'ddU':p['ms']**2+2*p['a3']*s+3*p['a4']*s*s}

def vector_potential(name: str, omega, rho) -> dict:
    p=model(name);w,r=np.broadcast_arrays(np.asarray(omega,float),np.asarray(rho,float))
    if not np.all(np.isfinite(w)) or not np.all(np.isfinite(r)):raise ValueError('finite vector fields required')
    return {'W':p['mw']**2*w*w/2+p['mr']**2*r*r/2+p['b4']*w**4/4+p['cross']*w*w*r*r,
            'dw':p['mw']**2*w+p['b4']*w**3+2*p['cross']*w*r*r,
            'dr':p['mr']**2*r+2*p['cross']*w*w*r,
            'dww':p['mw']**2+3*p['b4']*w*w+2*p['cross']*r*r,
            'drr':p['mr']**2+2*p['cross']*w*w,
            'dwr':4*p['cross']*w*r}

def rearrangement(name: str,density,scalar_density_value,phi,omega,rho,proton_density):
    """Local DD rearrangement in MeV; identical addition for n and p."""
    c=couplings(name,density)
    return c['dgw']*omega*density+c['dgr']*rho*(2*np.asarray(proton_density)-density)/2-c['dgs']*phi*scalar_density_value

def baryons(name: str,density: float,proton_fraction: float) -> dict:
    if not math.isfinite(density) or density<0 or not math.isfinite(proton_fraction) or not 0<=proton_fraction<=1:
        raise ValueError('finite nonnegative density and 0<=Yp<=1 required')
    p=model(name);c={k:float(v) for k,v in couplings(name,density).items()}
    nn,np_=density*(1-proton_fraction),density*proton_fraction;n3=np_-nn;hc3=HBARC**3
    def kinetic(phi):
        mn,mp=p['mn']-c['gs']*phi,p['mp']-c['gs']*phi
        sn,sp=scalar_density(nn,mn/HBARC),scalar_density(np_,mp/HBARC)
        return free_fermion(nn,mn),free_fermion(np_,mp),sn+sp
    def scalar_residual(phi):return float(scalar_potential(name,phi)['dU'])/hc3-c['gs']*kinetic(phi)[2]
    if density==0:phi=omega=rho=0.;roots=[0.]
    else:
        grid=np.linspace(0,min(p['mn'],p['mp'])/c['gs']*(1-1e-11),65)
        values=[scalar_residual(s) for s in grid];roots=[]
        for i in range(len(grid)-1):
            if values[i]*values[i+1]<=0:
                x=brentq(scalar_residual,grid[i],grid[i+1],xtol=1e-12)
                step=max(1e-5,abs(x)*1e-6)
                if (scalar_residual(x+step)-scalar_residual(x-step))/(2*step)>0:roots.append(x)
        if not roots:raise RuntimeError('no locally stable positive-mass scalar root')
        phi=min(roots,key=lambda x:kinetic(x)[0][0]+kinetic(x)[1][0]+float(scalar_potential(name,x)['U'])/hc3)
        def r_of_w(w):return c['gr']*n3*hc3/(2*(p['mr']**2+2*p['cross']*w*w))
        omega=brentq(lambda w:float(vector_potential(name,w,r_of_w(w))['dw'])-c['gw']*density*hc3,0,c['gw']*density*hc3/p['mw']**2*(1+1e-10),xtol=1e-12)
        rho=r_of_w(omega)
    kn,kp,ns=kinetic(phi);us=scalar_potential(name,phi);vw=vector_potential(name,omega,rho)
    sr=float(rearrangement(name,density,ns,phi,omega,rho,np_))
    mun=kn[2]+c['gw']*omega-c['gr']*rho/2+sr
    mup=kp[2]+c['gw']*omega+c['gr']*rho/2+sr
    eps=kn[0]+kp[0]+(float(us['U'])-float(vw['W']))/hc3+c['gw']*omega*density+c['gr']*rho*n3/2
    pressure=kn[1]+kp[1]+(float(vw['W'])-float(us['U']))/hc3+density*sr
    return {'model':p['name'],'n_b':density,'Y_p':proton_fraction,'epsilon_baryons':eps,'pressure_baryons':pressure,'mu_n':mun,'mu_p':mup,
            'effective_mass_n':p['mn']-c['gs']*phi,'effective_mass_p':p['mp']-c['gs']*phi,'phi':phi,'omega':omega,'rho':rho,'scalar_density':ns,'rearrangement':sr,
            'scalar_root_count':len(roots),'scalar_residual_fm3':scalar_residual(phi),
            'omega_residual_fm3':float(vw['dw'])/hc3-c['gw']*density,'rho_residual_fm3':float(vw['dr'])/hc3-c['gr']*n3/2,
            'euler_residual':eps+pressure-nn*mun-np_*mup}

def beta_state(name: str,density: float,include_muons: bool=True) -> dict:
    """Constrained neutral npe(mu) minimum, including absent-species boundaries."""
    if not math.isfinite(density) or density<=0:raise ValueError('positive finite density required')
    def lepton_at_charge(charge):
        if charge==0:return ELECTRON_MASS,leptons(ELECTRON_MASS,include_muons)
        upper=free_fermion(charge,ELECTRON_MASS)[2]
        if include_muons and upper>105.6583755:
            mu=brentq(lambda q:leptons(q,True)['n_e']+leptons(q,True)['n_mu']-charge,
                       ELECTRON_MASS,upper*(1+1e-13),xtol=2e-12)
        else:mu=upper
        return mu,leptons(mu,include_muons)
    def beta_residual(y):
        s=baryons(name,density,y);mu,_=lepton_at_charge(density*y)
        return s['mu_n']-s['mu_p']-mu
    f0,f1=beta_residual(0.),beta_residual(1.)
    if f0<=0:yp,branch=0.,'NEUTRON_ONLY_KKT'
    elif f1>=0:yp,branch=1.,'PROTON_ONLY_KKT'
    else:yp,branch=brentq(beta_residual,0.,1.,xtol=1e-13),'INTERIOR_BETA'
    s=baryons(name,density,yp);mu,l=lepton_at_charge(density*yp)
    eps=s['epsilon_baryons']+l['epsilon_leptons'];pressure=s['pressure_baryons']+l['pressure_leptons']
    # The baryon multiplier at y=1 is mu_p+mu_e; the neutron chemical
    # potential is a one-sided insertion threshold when that species is absent.
    mu_b=s['mu_p']+mu if yp==1 else s['mu_n']
    return dict(s,**l,Y_n=1-yp,Y_e=l['n_e']/density,Y_mu=l['n_mu']/density,mu_e=mu,mu_B=mu_b,
                epsilon=eps,pressure=pressure,composition_branch=branch,
                beta_residual_mev=s['mu_n']-s['mu_p']-mu,
                charge_residual=density*yp-l['n_e']-l['n_mu'],
                total_euler_residual=eps+pressure-density*mu_b)


def symmetry_energy(name: str,density: float,step: float=.001) -> float:
    """Quadratic isospin coefficient: S=(1/4)d(mu_n-mu_p)/d(delta).

    A five-point derivative avoids subtracting nucleon rest energies. The
    split bare masses of DD2 are retained; odd isospin terms are not assumed zero.
    """
    if not math.isfinite(step) or not 0<step<.25:raise ValueError('finite 0<isospin step<.25 required')
    def difference(delta):
        s=baryons(name,density,(1-delta)/2)
        return s['mu_n']-s['mu_p']
    return (difference(-2*step)-8*difference(-step)+8*difference(step)-difference(2*step))/(48*step)
