"""Independent pressure-coordinate TOV + unscaled W,V Cowling comparison."""
from pathlib import Path
import sys,os,json,time,hashlib,traceback
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python_frontend'))
sys.path.insert(0,str(ROOT/'python_frontend/tests'))
import numpy as np
from scipy.integrate import solve_ivp
from test_stellar_sharp import contract
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,matching_residual
from compose_eos.stellar_cowling import CowlingSettings,C_KM_S
from compose_eos.stellar_background import MEV_FM3_TO_GEOM_KM2
OUT=Path(__file__).parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run():
 bar=WholeCrustBarotrope.from_contract(contract())
 eos=SharpContractEnthalpyEOS(bar)
 bg=solve_sharp_background(eos,5.6,r0_km=.001,rtol=2e-11,atol=1e-13,max_step_h=.001)
 cfg=CowlingSettings(rtol=2e-9,atol=1e-12,surface_h=1e-5,max_step_h=.0005,match_h_fraction=.15)
 ps=eos.state_at_enthalpy(cfg.surface_h)['p']
 pt=eos.transitions[0]['p']
 pc0=eos.state_at_enthalpy(bg.h_start)['p']
 ec=eos.state_at_enthalpy(bg.h_c)['epsilon']*MEV_FM3_TO_GEOM_KM2
 geom0=[.001,4*np.pi*ec*.001**3/3]
 geom_segments=[]
 for a,b,side in [(pc0,pt,'high'),(pt,ps,'low')]:
  def rhs(p,y):
   s=bar.state_at_pressure(max(p,ps),transition_side=side)
   r,m=y;pp=p*MEV_FM3_TO_GEOM_KM2;ee=s['epsilon']*MEV_FM3_TO_GEOM_KM2
   dr=-r*(r-2*m)/((pp+ee)*(m+4*np.pi*r**3*pp))*MEV_FM3_TO_GEOM_KM2
   return [dr,4*np.pi*r*r*ee*dr]
  sol=solve_ivp(rhs,(a,b),geom0,method='DOP853',rtol=2e-11,atol=1e-13,max_step=.005,dense_output=True)
  if not sol.success:raise RuntimeError(sol.message)
  geom_segments.append(sol);geom0=sol.y[:,-1]
 # Exterior Phi constant from full enthalpy star; pressure-coordinate background
 # truncated at the common mode surface. M/R are compared at exactly that surface.
 rs,ms=geom0
 etgt=bg.state_at_enthalpy(cfg.surface_h)
 metric_constant=.5*np.log1p(-2*bg.mass_geom_km/bg.radius_km)
 def state(p,side):
  s=bar.state_at_pressure(float(p),transition_side=side)
  gs=geom_segments[0 if p>pt or (p==pt and side=='high') else 1]
  r,m=gs.sol(p); pp=p*MEV_FM3_TO_GEOM_KM2
  phi=metric_constant-s['enthalpy']
  gp=(m+4*np.pi*r**3*pp)/(r*(r-2*m))
  return s,float(r),float(m),float(phi),float(gp),float((1-2*m/r)**-.5)
 pm=eos.state_at_enthalpy(cfg.match_h_fraction*bg.h_c)['p']
 rows=[]
 for conv in ['slow','rapid']:
  for f in [100.,300.,900.]:
   omega=2*np.pi*f/C_KM_S
   def path(a,b,y):
    direction=np.sign(b-a)
    bounds=[a,*([pt] if min(a,b)<pt<max(a,b) else []),b]
    for pa,pb in zip(bounds[:-1],bounds[1:]):
     side='high' if (pa+pb)/2>pt else 'low'
     def rhs(p,y):
      s,r,m,phi,gp,el=state(max(p,ps),side)
      W,V=y
      dW=(omega**2*r*r*el*np.exp(-2*phi)*V+gp*W)/s['cs2_over_c2']-6*el*V
      dV=2*gp*V-el*W/(r*r)
      dpdr=-(s['epsilon']+s['p'])*gp
      return np.array([dW,dV])/dpdr
     sol=solve_ivp(rhs,(pa,pb),y,method='DOP853',rtol=2e-9,atol=1e-12,max_step=.002)
     if not sol.success:raise RuntimeError(sol.message)
     y=sol.y[:,-1]
     if pb==pt and pb!=b:
      fromside,toside=('high','low') if direction<0 else ('low','high')
      ss,r,m,phi,gp,el=state(pt,fromside)
      ts,*_=state(pt,toside)
      # Independent 2x2 physical constraints in unscaled W,V:
      # traction = (eps+P)*(Phi' W+omega^2*r^2*e^(Lambda-2Phi)V)
      acoef=omega**2*r*r*el*np.exp(-2*phi)
      traction=(ss['epsilon']+pt)*(gp*y[0]+acoef*y[1])
      matrix=np.array([[(ts['epsilon']+pt)*gp,(ts['epsilon']+pt)*acoef],
                       [1.,0.] if conv=='slow' else [0.,1.]])
      y=np.linalg.solve(matrix,[traction,y[0] if conv=='slow' else y[1]])
    return y
   ui=path(pc0,pm,[1.,-1/(2*.001)])
   s,r,m,phi,gp,el=state(ps,'low')
   uo=path(ps,pm,[1.,-gp/(omega**2*el*np.exp(-2*phi)*r*r)])
   _,r,*_=state(pm,'low' if pm<pt else 'high')
   ui=ui/np.array([r**3,r**2]);uo=uo/np.array([r**3,r**2])
   independent=float(np.linalg.det(np.column_stack((ui,uo)))/(np.linalg.norm(ui)*np.linalg.norm(uo)))
   production=matching_residual(bg,f,cfg,conversion=conv)
   rows.append(dict(conversion=conv,frequency_hz=f,pressure_unscaled_residual=independent,
                    enthalpy_scaled_residual=production,absolute_difference=abs(independent-production)))
 result={'schema':'G.independent.sharp-control.v1','scope':'Synthetic stable finite density jump; pressure-coordinate background and unscaled W,V with direct linear junction solve',
  'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'python_frontend/compose_eos/stellar_sharp.py',Path(__file__)]},
  'background_at_common_surface':{'r_pressure_km':rs,'r_enthalpy_km':etgt['r_km'],'mass_pressure_geom_km':ms,'mass_enthalpy_geom_km':etgt['m_geom_km'],
      'relative_radius_difference':abs(rs-etgt['r_km'])/etgt['r_km'],'relative_mass_difference':abs(ms-etgt['m_geom_km'])/etgt['m_geom_km']},
  'controls':rows,'max_absolute_mode_residual_difference':max(x['absolute_difference'] for x in rows)}
 assert result['background_at_common_surface']['relative_radius_difference']<2e-7
 assert result['background_at_common_surface']['relative_mass_difference']<2e-7
 assert result['max_absolute_mode_residual_difference']<2e-6
 result['status']='PASS'
 return result
start=time.monotonic()
try:
 result=run();result.update(actual_exit_code=0,wall_seconds=time.monotonic()-start)
except Exception as e:
 result={'status':'FAIL','actual_exit_code':1,'error':str(e),'traceback':traceback.format_exc(),'wall_seconds':time.monotonic()-start}
(OUT/'INDEPENDENT_SHARP_CONTROL.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
sys.exit(result['actual_exit_code'])
