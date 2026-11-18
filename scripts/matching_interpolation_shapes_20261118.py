"""Exploratory cold-EOS matching closures. Never reinterpret a sharp Maxwell doublet.

P(mu) = P_L + integral n(mu)dmu, epsilon(mu) = mu*n(mu)-P(mu).
Only the low-density phenomenological bridge is replaced in memory.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

_GX, _GW = leggauss(64)
_GX, _GW = (_GX+1)/2, _GW/2


def _arr(x):
    x = np.asarray(x, dtype=float)
    if np.any(~np.isfinite(x)) or np.any(x<0) or np.any(x>1):
        raise ValueError("unit-coordinate support required")
    return x


class CausalShapeInterval:
    """Exact endpoint-determined monotone y=n/mu bridge (power or rational)."""
    def __init__(self, n_left, p_left, epsilon_left, n_right, p_right, epsilon_right,
                 *, family):
        v=np.asarray([n_left,p_left,epsilon_left,n_right,p_right,epsilon_right],dtype=float)
        if v.shape!=(6,) or np.any(~np.isfinite(v)) or np.any(v<=0):
            raise ValueError("positive finite endpoint states required")
        if family not in ("power","rational"):
            raise ValueError("explicit known shape required")
        (self.n_left,self.p_left,self.epsilon_left,
         self.n_right,self.p_right,self.epsilon_right)=map(float,v)
        self.mu_left=(self.epsilon_left+self.p_left)/self.n_left
        self.mu_right=(self.epsilon_right+self.p_right)/self.n_right
        self.width=self.mu_right-self.mu_left
        self.y_left=self.n_left/self.mu_left
        self.y_right=self.n_right/self.mu_right
        self.delta_y=self.y_right-self.y_left
        self.delta_p=self.p_right-self.p_left
        if min(self.width,self.delta_y,self.delta_p)<=0:
            raise ValueError("positive ordered endpoints and increasing n/mu required")
        self.causal_lower=self.y_left*self.width*(self.mu_left+self.width/2)
        self.causal_upper=self.y_right*self.width*(self.mu_left+self.width/2)
        if not self.causal_lower<self.delta_p<self.causal_upper:
            raise ValueError("strict general causal endpoint area bounds not met")
        self.family=family
        target=(self.delta_p-self.causal_lower)/(self.width*self.delta_y)
        if family=="power":
            def area(a):return self.mu_left/(a+1)+self.width/(a+2)
            self.shape_parameter=float(brentq(lambda a:area(a)-target,1e-9,1e8,xtol=2e-13))
        else:
            def area(logq):
                q=np.exp(logq)
                x=_GX
                phi=x/(q+(1-q)*x)
                return float(np.dot(_GW,(self.mu_left+self.width*x)*phi))
            self.shape_parameter=float(np.exp(brentq(lambda t:area(t)-target,-25,25,xtol=2e-13)))
        self.pressure_integral_residual=self._raw_pressure(1.)-self.p_right

    def _phi_and_prime(self,x):
        x=_arr(x)
        if self.family=="power":
            a=self.shape_parameter
            phi=x**a
            deriv=np.where(x==0,0 if a>1 else (1 if a==1 else np.inf),a*np.maximum(x,1e-300)**(a-1))
            return phi,deriv
        q=self.shape_parameter
        den=q+(1-q)*x
        return x/den,q/(den*den)

    def _integrals(self,x):
        x=_arr(x)
        if self.family=="power":
            a=self.shape_parameter
            return x**(a+1)/(a+1),x**(a+2)/(a+2)
        q=self.shape_parameter
        m=1-q
        if abs(m)<0.04:
            # Avoid catastrophic cancellation of rational primitives around q=1.
            nodes=np.expand_dims(x,-1)*_GX
            den=q+m*nodes
            phi=nodes/den
            w=np.expand_dims(x,-1)*_GW
            return np.sum(w*phi,axis=-1),np.sum(w*nodes*phi,axis=-1)
        u=np.log1p(m*x/q)
        f0=x/m-q*u/m**2
        f1=x*x/(2*m)-q*x/m**2+q*q*u/m**3
        return f0,f1

    def _raw_pressure(self,x):
        x=_arr(x)
        f0,f1=self._integrals(x)
        return self.p_left+self.width*(self.y_left*(self.mu_left*x+self.width*x*x/2)
                   +self.delta_y*(self.mu_left*f0+self.width*f1))

    def pressure(self,t):
        x=_arr(t)
        p=self._raw_pressure(x)
        return np.where(x==0,self.p_left,np.where(x==1,self.p_right,p)).item() if x.ndim==0 else np.where(x==0,self.p_left,np.where(x==1,self.p_right,p))

    def number(self,t):
        x=_arr(t)
        mu=self.mu_left+self.width*x
        phi,_=self._phi_and_prime(x)
        n=mu*(self.y_left+self.delta_y*phi)
        n=np.where(x==0,self.n_left,np.where(x==1,self.n_right,n))
        return n.item() if x.ndim==0 else n

    def susceptibility(self,t):
        x=_arr(t)
        mu=self.mu_left+self.width*x
        phi,dp=self._phi_and_prime(x)
        value=(self.y_left+self.delta_y*phi)+mu*self.delta_y*dp/self.width
        return value.item() if x.ndim==0 else value

    def evaluate_mu(self,chemical_potential):
        mu=np.asarray(chemical_potential,dtype=float)
        if np.any(~np.isfinite(mu)) or np.any(mu<self.mu_left) or np.any(mu>self.mu_right):
            raise ValueError("chemical potential outside bridge support")
        x=(mu-self.mu_left)/self.width
        n=self.number(x);p=self.pressure(x)
        e=np.where(x==0,self.epsilon_left,np.where(x==1,self.epsilon_right,mu*n-p))
        c2=n/(mu*self.susceptibility(x))
        return {key: (np.asarray(val).item() if np.asarray(val).ndim==0 else val) for key,val in
                dict(n=n,p=p,epsilon=e,mu=mu,cs2_over_c2=c2).items()}

    def validate(self):
        x=np.linspace(0,1,4001)
        st=self.evaluate_mu(self.mu_left+self.width*x)
        c=np.asarray(st["cs2_over_c2"])
        full=np.all(np.isfinite(c)) and np.min(c)>0 and np.max(c)<=1+5e-12
        err=abs(self.pressure_integral_residual)
        bound=1e-9*max(self.delta_p,self.causal_upper)
        if not (np.min(np.diff(st["n"]))>0 and err<bound):
            full=False
        return dict(status="PASS" if full else "FAIL",family=self.family,
                    parameter=self.shape_parameter,min_cs2_over_c2=float(np.min(c)),
                    max_cs2_over_c2=float(np.max(c)),endpoint_pressure_residual=float(err),
                    endpoint_pressure_bound=float(bound),
                    first_law="ANALYTIC_BY_CONSTRUCTION; FINITE_DIFFERENCE_TESTED",
                    causal_margin="mu^2*delta_y*Phi_prime / width >=0")

    def contract(self):
        return dict(kind="exploratory_causal_shape",family=self.family,
                    parameter=self.shape_parameter,
                    endpoint_states=[dict(n=self.n_left,p=self.p_left,epsilon=self.epsilon_left),
                                     dict(n=self.n_right,p=self.p_right,epsilon=self.epsilon_right)],
                    physical_scope="lower bridge only; not microscopic TF source",
                    validation=self.validate())


def naive_loglog_polytrope_endpoint(n_left,p_left,epsilon_left,n_right,p_right,epsilon_right):
    """Negative control: P(n) log-log line, first-law epsilon reconstructed from L."""
    gamma=np.log(p_right/p_left)/np.log(n_right/n_left)
    ratio=n_right/n_left
    if abs(gamma-1)<1e-10:
        integral=p_left/n_left*np.log(ratio)
    else:
        integral=p_left/n_left*np.expm1((gamma-1)*np.log(ratio))/(gamma-1)
    derived_epsilon=n_right*(epsilon_left/n_left+integral)
    return dict(gamma=float(gamma),pressure_endpoint_residual=float(p_left*ratio**gamma-p_right),
                energy_endpoint_residual_MeV_fm3=float(derived_epsilon-epsilon_right),
                right_mu_shift_MeV=float((derived_epsilon-epsilon_right)/n_right),
                thermodynamically_valid_for_fixed_endpoints=abs(derived_epsilon-epsilon_right)<1e-10)


class SlopeMatchedBernsteinInterval(CausalShapeInterval):
    """Nonnegative Bernstein derivative matches endpoint c_s^2 and pressure area.

    The independent data are exact endpoint (n,P,epsilon) and explicitly
    labeled reference-side sound speeds from the adjacent source intervals.
    An LP certificate with nonnegative Bernstein coefficients enforces
    Phi'(x)>=0 for all x. No microscopic cells are reconstructed.
    """
    def __init__(self,n_left,p_left,epsilon_left,n_right,p_right,epsilon_right,
                 *,left_cs2,right_cs2):
        from scipy.optimize import linprog, minimize, Bounds, LinearConstraint
        v=np.asarray([n_left,p_left,epsilon_left,n_right,p_right,epsilon_right],float)
        if np.any(~np.isfinite(v)) or np.any(v<=0):
            raise ValueError("positive finite thermodynamic endpoints required")
        (self.n_left,self.p_left,self.epsilon_left,
         self.n_right,self.p_right,self.epsilon_right)=map(float,v)
        self.mu_left=(self.epsilon_left+self.p_left)/self.n_left
        self.mu_right=(self.epsilon_right+self.p_right)/self.n_right
        self.width=self.mu_right-self.mu_left
        self.y_left=self.n_left/self.mu_left
        self.y_right=self.n_right/self.mu_right
        self.delta_y=self.y_right-self.y_left
        self.delta_p=self.p_right-self.p_left
        if min(self.width,self.delta_y,self.delta_p)<=0:
            raise ValueError("ordered nondegenerate endpoints required")
        self.causal_lower=self.y_left*self.width*(self.mu_left+self.width/2)
        self.causal_upper=self.y_right*self.width*(self.mu_left+self.width/2)
        if not self.causal_lower<self.delta_p<self.causal_upper:
            raise ValueError("general causal endpoint area bounds fail")
        if (not np.isfinite(left_cs2+right_cs2) or
                not 0<left_cs2<=1 or not 0<right_cs2<=1):
            raise ValueError("finite strictly positive causal side speeds required")
        self.left_cs2=float(left_cs2)
        self.right_cs2=float(right_cs2)
        self.g0=self.width*self.y_left/(self.mu_left*self.delta_y)*(1/self.left_cs2-1)
        self.g1=self.width*self.y_right/(self.mu_right*self.delta_y)*(1/self.right_cs2-1)
        if min(self.g0,self.g1)<0:
            raise ValueError("negative derivative at endpoint prohibited")
        target=(self.delta_p-self.causal_lower)/(self.width*self.delta_y)
        accepted=None
        attempts=[]
        for D in (4,6,8,12,16,24,32):
            i=np.arange(D+1,dtype=float)
            mean=(i+1)/(D+2)
            mean2=(i+1)*(i+2)/((D+2)*(D+3))
            weights=(self.mu_left*(1-mean)+self.width*(1-mean2)/2)/(D+1)
            A=np.vstack([np.ones(D+1),weights])
            rhs=np.array([D+1,target])
            bounds=[(self.g0,self.g0)]+[(0,None)]*(D-1)+[(self.g1,self.g1)]
            lp=linprog(np.zeros(D+1),A_eq=A,b_eq=rhs,bounds=bounds,method="highs")
            attempts.append(dict(degree=D,feasible=bool(lp.success)))
            if lp.success:
                accepted=(D,A,rhs,bounds,lp.x)
                break
        if accepted is None:
            raise ValueError("Bernstein finite-degree nonnegative derivative LP infeasible")
        D,A,rhs,bounds,b=accepted
        def curvature(z):
            dif=np.diff(z)
            return float(np.dot(dif,dif))
        opt=minimize(curvature,b,method="SLSQP",
            bounds=Bounds([q[0] for q in bounds],
                          [np.inf if q[1] is None else q[1] for q in bounds]),
            constraints=[LinearConstraint(A,rhs,rhs)],
            options=dict(maxiter=500,ftol=1e-12))
        if (opt.success and np.min(opt.x)>=-1e-10 and
                np.max(abs(A@opt.x-rhs))<=1e-8):
            b=np.maximum(opt.x,0)
            self.selection="minimum discrete derivative roughness, SLSQP after LP feasibility"
        else:
            self.selection="nonnegative LP feasible certificate"
        self.degree=D
        self.shape_parameter=D
        self.family="bernstein_slope_matched"
        self.bernstein_coefficients=np.asarray(b,float)
        self.lp_attempts=attempts
        self.area_constraint_max_abs_residual=float(np.max(abs(A@b-rhs)))
        self.pressure_integral_residual=self._raw_pressure(1.)-self.p_right

    def _phi_and_prime(self,x):
        from scipy.special import betainc, comb
        xx=_arr(x)
        D=self.degree
        i=np.arange(D+1)
        t=np.expand_dims(xx,-1)
        c=self.bernstein_coefficients
        phi=np.sum(c/(D+1)*betainc(i+1,D-i+1,t),axis=-1)
        B=comb(D,i)*(t**i)*((1-t)**(D-i))
        derivative=np.sum(c*B,axis=-1)
        return phi,derivative

    def _integrals(self,x):
        from scipy.special import betainc
        xx=_arr(x)
        D=self.degree
        i=np.arange(D+1)
        t=np.expand_dims(xx,-1)
        c=self.bernstein_coefficients/(D+1)
        phi=np.sum(c*betainc(i+1,D-i+1,t),axis=-1)
        M1=np.sum(c*(i+1)/(D+2)*betainc(i+2,D-i+1,t),axis=-1)
        M2=np.sum(c*(i+1)*(i+2)/((D+2)*(D+3))*betainc(i+3,D-i+1,t),axis=-1)
        return xx*phi-M1,xx*xx*phi/2-M2/2

    def validate(self):
        check=super().validate()
        x=self.evaluate_mu(self.mu_left)
        y=self.evaluate_mu(self.mu_right)
        check.update(degree=self.degree,bernstein_coefficients=self.bernstein_coefficients.tolist(),
            endpoint_cs2_expected=[self.left_cs2,self.right_cs2],
            endpoint_cs2_found=[x["cs2_over_c2"],y["cs2_over_c2"]],
            endpoint_derivative_g=[self.g0,self.g1],
            area_constraint_max_abs_residual=self.area_constraint_max_abs_residual,
            selected_basis=self.selection,lp_attempts=self.lp_attempts)
        if (not np.isclose(x["cs2_over_c2"],self.left_cs2,rtol=1e-8,atol=1e-12)
          or not np.isclose(y["cs2_over_c2"],self.right_cs2,rtol=1e-8,atol=1e-12)
          or self.area_constraint_max_abs_residual>1e-8
          or np.min(self.bernstein_coefficients)<-1e-10):
            check["status"]="FAIL"
        return check
