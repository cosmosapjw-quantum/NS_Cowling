"""Cold catalysed, discrete-nucleus outer crust from nuclear masses.

Pearson et al. MNRAS481,2994(2018), Sec2.1/AppB; PRC83,065810(2011).
MeV, fm, c retained through rest energies and hbar*c constants. Not inner crust.
The mass catalogue, not a provider EOS table, is the microscopic input.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import csv
import math
import numpy as np
from scipy.optimize import brentq
from .eos_data import EOSData

HBARC = 197.3269804
ME = 0.51099895
MN = 939.56542052
MP = 938.27208816
U = 931.49410242
ALPHA = 1 / 137.035999084
E2 = ALPHA * HBARC
BCC = 0.895929255682
U1 = 0.5113875
SOURCE = "Pearson2018 Sec2.1, AppendixB; Pearson2011 Eqs16-22; AME2016/HFB22/Welker2017"
MASS_PATH = Path(__file__).resolve().parents[2] / "data/outer_crust/hfb22_ame2016_cu.csv"

@dataclass(frozen=True)
class MassCatalogue:
    z: np.ndarray
    a: np.ndarray
    excess: np.ndarray
    radius: np.ndarray
    source: tuple[str, ...]

    def __post_init__(self):
        values = (self.z, self.a, self.excess, self.radius)
        if any(np.ndim(x) != 1 for x in values) or not len(self.z):
            raise ValueError("nonempty one-dimensional mass arrays required")
        if any(len(x) != len(self.z) for x in values) or len(self.source) != len(self.z):
            raise ValueError("mass catalogue lengths differ")
        if not all(np.all(np.isfinite(x)) for x in values):
            raise ValueError("mass catalogue must be finite")
        if np.any(self.z <= 0) or np.any(self.a < self.z) or np.any(self.radius < 0):
            raise ValueError("invalid nuclear mass catalogue")
        if np.any(self.z != np.rint(self.z)) or np.any(self.a != np.rint(self.a)):
            raise ValueError("integer nuclear charges/mass numbers required")
        if len(set(zip(self.z, self.a))) != len(self.z):
            raise ValueError("duplicate nuclear mass entry")

    @property
    def atomic_corrected(self):
        return self.a*U + self.excess + 1.44381e-5*self.z**2.39 + 1.55468e-12*self.z**5.35

@lru_cache(maxsize=1)
def masses() -> MassCatalogue:
    with MASS_PATH.open() as stream:
        rows = list(csv.DictReader(stream))
    array = lambda key: np.array([float(r[key]) for r in rows])
    return MassCatalogue(array("Z"), array("A"), array("mass_excess_mev"),
                         array("rch_fm"), tuple(r["source"] for r in rows))

def _energy_pressure(ne, z, a, excess, radius):
    """Vectorized energy/pressure. Differentiation holds isotope A,Z fixed."""
    ne, z, a, excess, radius = np.broadcast_arrays(ne, z, a, excess, radius)
    x = HBARC*(3*np.pi**2*ne)**(1/3)/ME
    gamma = np.sqrt(1+x*x)
    inv = np.arcsinh(x)
    # In the supported domain x >= 0.1, these expressions have no severe cancellation.
    factor = ME**4/(8*np.pi**2*HBARC**3)
    ee = factor*(x*(1+2*x*x)*gamma-inv)
    pe = factor/3*(x*(2*x*x-3)*gamma+3*inv)
    aa = x*gamma-inv
    k = ALPHA*ME**4/(4*np.pi**3*HBARC**3)
    ex = k*(1.5*aa*aa-x**4)
    px = k*(2*aa*x**3/gamma-x**4/3-1.5*aa*aa)
    # Pearson2018 B11, with derivative of log10(x) explicitly retained.
    corr_pref = ALPHA**2*ME**4/(12*np.pi**4*HBARC**3)
    corr = corr_pref*x**4*(-12.51+np.log10(x))
    pcorr = corr/3+corr_pref*x**4/(3*np.log(10))
    mprime = a*U+excess+1.44381e-5*z**2.39+1.55468e-12*z**5.35
    radius_cell = (3*z/(4*np.pi*ne))**(1/3)
    # B6/B8 convention uses M' (electron-binding-corrected atomic mass).
    hw = HBARC*np.sqrt(4*np.pi*z*E2*ne/mprime)
    gamma_pl = z*z*E2/(radius_cell*hw)
    s = 1/(1+0.01*np.log(z)**1.5+0.097/z**2)
    b1 = 1-1.1866/z**0.267+0.27/z
    b2 = 1+2.25/z**(1/3)*(1+0.684*z**5+0.222*z**6)/(1+0.222*z**6)
    b3 = 41.5/(1+np.log(z))
    b4 = 0.395*np.log(z)+0.347/z**1.5
    finf = 54/175*(12/np.pi)**(1/3)*ALPHA*z**(2/3)*b1*np.sqrt(1+b2/x**2)
    screening_a = (b3+17.9*x*x)/(1+b4*x*x)
    bb = (0.205/gamma_pl)**s
    base = (4*np.pi/3)**(1/3)*E2*z**(2/3)*ne**(4/3)
    screen = -base*finf*(1+screening_a*bb)
    df = -finf*b2/(3*(x*x+b2))
    da = 2*x*x*(17.9-b4*b3)/(3*(1+b4*x*x)**2)
    ds = df*(1+screening_a*bb)+finf*bb*(da+screening_a*s/6)
    pscreen = screen/3-base*ds
    lattice = -BCC*base
    finite_size = -lattice*5*radius**2/(9*radius_cell**2)
    zero_point = ne/z*1.5*U1*hw
    ionic = lattice+finite_size+zero_point
    pionic = lattice/3+finite_size+zero_point/2
    electron = ee+ex+corr+screen
    electron_pressure = pe+px+pcorr+pscreen
    epsilon = ne/z*(mprime-z*ME)+electron+ionic
    pressure = electron_pressure+pionic
    mu_e = (electron+electron_pressure)/ne
    return epsilon, pressure, mu_e

def isotope_state(ne: float, z: int, a: int, excess: float, radius: float) -> dict:
    if not all(math.isfinite(v) for v in (ne,z,a,excess,radius)):
        raise ValueError("finite isotope inputs required")
    if ne <= 0 or z <= 0 or a < z or radius < 0 or z != int(z) or a != int(a):
        raise ValueError("invalid isotope state")
    x = HBARC*(3*np.pi**2*ne)**(1/3)/ME
    if not 0.1 <= x <= 100:
        raise ValueError("outer-crust computational support requires 0.1 <= x_e <= 100")
    e,p,mu = _energy_pressure(ne,z,a,excess,radius)
    nb = ne*a/z
    return {"n_b":nb,"n_e":ne,"pressure":float(p),"epsilon":float(e),
            "g":float((e+p)/nb),"mu_e":float(mu),"A":int(a),"Z":int(z),
            "Y_n":0.0,"Y_p":0.0,"Y_h":1.0,"Y_e":z/a,"Y_mu":0.0,
            "cluster_number_fraction":1/a}

class OuterCrust:
    """Global finite-catalogue minimization at fixed pressure, not fixed density."""
    def __init__(self, catalogue: MassCatalogue | None = None):
        self.catalogue = masses() if catalogue is None else catalogue

    def candidates(self, pressure: float):
        if not math.isfinite(pressure) or not 1e-13 <= pressure <= 1e-2:
            raise ValueError("supported pressure range is [1e-13,1e-2] MeV fm^-3")
        c = self.catalogue
        # Bracket all isotope roots on the positive-compressibility branch.
        # Vector bisection gives an auditable bound independent of initial guesses.
        lo = np.full(len(c.z), 1e-11)
        hi = np.full(len(c.z), 2e-3)
        plo = _energy_pressure(lo,c.z,c.a,c.excess,c.radius)[1]
        phi = _energy_pressure(hi,c.z,c.a,c.excess,c.radius)[1]
        if np.any(plo >= pressure) or np.any(phi <= pressure):
            raise RuntimeError("pressure bracket does not cover every nuclear candidate")
        for _ in range(44):
            mid = np.sqrt(lo*hi)
            pm = _energy_pressure(mid,c.z,c.a,c.excess,c.radius)[1]
            lo = np.where(pm < pressure, mid, lo)
            hi = np.where(pm >= pressure, mid, hi)
        ne = np.sqrt(lo*hi)
        ee, pp, mu = _energy_pressure(ne,c.z,c.a,c.excess,c.radius)
        nb = ne*c.a/c.z
        return ne, ee, pp, mu, (ee+pp)/nb

    def state(self, pressure: float, *, allow_above_drip: bool = False) -> dict:
        ne, ee, pp, mu, gg = self.candidates(pressure)
        i = int(np.argmin(gg))
        c = self.catalogue
        if c.radius[i] == 0:
            raise RuntimeError("selected nucleus lacks a source charge radius")
        if gg[i] > MN+1e-9 and not allow_above_drip:
            raise ValueError("requested pressure is above neutron drip")
        row = isotope_state(float(ne[i]),int(c.z[i]),int(c.a[i]),float(c.excess[i]),float(c.radius[i]))
        row["pressure"] = pressure
        row["g"] = (row["epsilon"]+pressure)/row["n_b"]
        row["mass_source"] = c.source[i]
        row["mass_index"] = i
        row["second_best_gap_mev"] = float(np.partition(gg,1)[1]-gg[i]) if len(gg)>1 else math.inf
        row["pressure_relative_residual"] = float(pp[i]/pressure-1)
        return row

    def neutron_drip(self) -> dict:
        root = brentq(lambda logp: self.state(math.exp(logp),allow_above_drip=True)["g"]-MN,
                      math.log(1e-5),math.log(2e-3),xtol=2e-12)
        return self.state(math.exp(root),allow_above_drip=True)

    def sequence(self, points: int = 257, p_min: float | None = None):
        """Resolve transitions bracketed by a logarithmic pressure census.

        The complete mass catalogue is checked at every sample and refined
        transition. A finite census can miss a thin layer; compare resolutions.
        At each crossing both one-sided density limits are exported.
        """
        if points < 16:
            raise ValueError("at least16 pressure samples required")
        drip = self.neutron_drip()
        if p_min is None:
            c=self.catalogue
            ii=np.flatnonzero((c.z==26)&(c.a==56))
            if len(ii)!=1: raise ValueError("default lower boundary requires Fe56")
            k=int(ii[0])
            p_min=isotope_state(1e-9*26/56,26,56,float(c.excess[k]),float(c.radius[k]))["pressure"]
        if not 1e-13 <= p_min < drip["pressure"]:
            raise ValueError("invalid lower pressure")
        logs = np.linspace(math.log(p_min),math.log(drip["pressure"]),points)
        rows = [self.state(math.exp(p),allow_above_drip=True) for p in logs]
        transitions=[]
        c=self.catalogue
        self.intermediate_layer_refinements=0
        def refine(left,right,depth=0):
            i,j=left["mass_index"],right["mass_index"]
            if i==j:
                return
            if depth>32:
                raise RuntimeError("nuclear transition refinement depth exceeded")
            def delta(lp):
                gg=self.candidates(math.exp(lp))[4]
                return float(gg[i]-gg[j])
            lp=brentq(delta,math.log(left["pressure"]),math.log(right["pressure"]),xtol=2e-12)
            p=math.exp(lp)
            ne,ee,pp,mu,gg=self.candidates(p)
            winner=int(np.argmin(gg))
            if gg[winner] < min(gg[i],gg[j])-1e-8:
                # A third nuclide wins at the apparent crossing: it is a
                # real intermediate layer, not a failed physical model.
                self.intermediate_layer_refinements+=1
                middle=self.state(p,allow_above_drip=True)
                refine(left,middle,depth+1)
                refine(middle,right,depth+1)
                return
            pair=[]
            for k in (i,j):
                r=isotope_state(float(ne[k]),int(c.z[k]),int(c.a[k]),float(c.excess[k]),float(c.radius[k]))
                r["pressure"]=p
                r["g"]=(r["epsilon"]+p)/r["n_b"]
                r.update(mass_source=c.source[k],mass_index=k,second_best_gap_mev=0.0,
                         pressure_relative_residual=float(pp[k]/p-1))
                pair.append(r)
            transitions.append({"pressure":p,"before":pair[0],"after":pair[1],
                                "g_difference_mev":float(gg[i]-gg[j])})
            rows.extend(pair)
        sampled=rows.copy()
        for left,right in zip(sampled[:-1],sampled[1:]):
            refine(left,right)
        rows.sort(key=lambda r:r["n_b"])
        # Preserve phase-side limits, but remove repeated sampled drip endpoints.
        clean=[]
        for r in rows:
            if not clean or r["n_b"] > clean[-1]["n_b"]*(1+1e-13):
                clean.append(r)
        return clean,transitions,drip

def eos_data(rows: list[dict]) -> EOSData:
    n=np.array([r["n_b"] for r in rows])
    value=lambda key:np.array([r[key] for r in rows])
    zero=np.zeros(len(rows))
    # Each discrete nucleus is a separate thermodynamic branch. Native stencils
    # must not cross density discontinuities at first-order layer transitions.
    labels={}; phases=[]
    for r in rows:
        key=(r["A"],r["Z"])
        if key not in labels: labels[key]=10+len(labels)
        phases.append(labels[key])
    return EOSData(nb=n,mn=MN,mp=MP,Q1=value("pressure")/n,Q2=zero.copy(),
        Q3=value("g")/MN-1,Q4=-value("mu_e")/MN,Q5=zero.copy(),
        Q6=value("epsilon")/(n*MN)-1,Q7=value("epsilon")/(n*MN)-1,q_extra=zero.copy(),
        Iphase=np.array(phases),Yn=zero.copy(),Yp=zero.copy(),Ye=value("Y_e"),
        Ymu=zero.copy(),Aav=value("A"),Zav=value("Z"),Yh=np.ones(len(rows)),
        composition_file_present=True,composition_available=True,beta_equilibrated=True)
