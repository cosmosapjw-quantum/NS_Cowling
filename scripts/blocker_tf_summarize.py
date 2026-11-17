#!/usr/bin/env python3
"""Summarize actual native samples and derivative diagnostics without EOS physics."""
import csv,json,hashlib,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_blockers_20261116/tf'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open()))
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def derivative_summary(mesh):
 base=BASE/f'branch_derivative_N{mesh}'
 if not (base/'SUMMARY.json').exists():return {'status':'NOT_RUN'}
 summary=json.loads((base/'SUMMARY.json').read_text());details=[];states=collections.Counter()
 for rp in sorted(base.glob('*/RECEIPT.json')):
  meta=json.loads(rp.read_text());lines=(rp.parent/'stdout.txt').read_text().splitlines()
  d=[x.split() for x in lines if x.startswith('DERIVATIVE ')];s=[x.split() for x in lines if x.startswith(('CENTER ','LOCAL '))]
  if not d:
   details.append({'id':meta['id'],'status':meta['status']});continue
  cols=('h_fm3','mu_from_energy_h','mu_from_energy_half_h','P_from_energy_h','P_from_energy_half_h','native_center_P','P_residual_h','P_residual_half_h','P_difference_halving_h')
  entry={'id':meta['id'],'status':meta['status'],'density_fm3':meta['density_fm3'],'model':meta['model'],'mesh':mesh,**dict(zip(cols,map(float,d[0][1:]))),'states':[]}
  for row in s:
   states[row[1]]+=1
   columns=('n_b','offset','epsilon','pressure','mu_b','KKT','radius_stationarity','radius_fm','weighted_shape_distance','relative_radius_change','scalar_residual','vector_residual','electron_residual')
   entry['states'].append({'kind':row[0],'numerical_status':row[1],'geometry':row[2],**dict(zip(columns,map(float,row[3:16]))),'iterations':int(row[16]),'evaluations':int(row[17]),'optimizer_mode':int(row[18])})
  entry['unchanged_geometry']=len(set(x['geometry'] for x in entry['states']))==1
  delta=next(x.split() for x in lines if x.startswith('REFINEMENT_DELTA '))
  entry['center_refinement_delta_epsilon']=float(delta[1]);entry['center_refinement_delta_pressure']=float(delta[2]);entry['center_refinement_delta_mu']=float(delta[3])
  entry['receipt_sha256']=sha(rp);entry['stdout_sha256']=sha(rp.parent/'stdout.txt');details.append(entry)
 valid=[x for x in details if 'P_residual_half_h' in x]
 return {'status':summary['status'],'jobs':len(details),'native_stationary_states':dict(states),'max_abs_P_residual_half_h':max((abs(x['P_residual_half_h']) for x in valid),default=None),'max_abs_P_change_halving_h':max((abs(x['P_difference_halving_h']) for x in valid),default=None),'max_weighted_shape_distance':max((s['weighted_shape_distance'] for x in valid for s in x['states']),default=None),'max_abs_relative_radius_change':max((abs(s['relative_radius_change']) for x in valid for s in x['states']),default=None),'geometry_changed_jobs':[x['id'] for x in valid if not x['unchanged_geometry']],'results':details}
def main():
 plan=json.loads((BASE/'PLAN.json').read_text())
 report={'baseline':plan['baseline'],'scope':'50 original N16 negative intervals; exact endpoint union; finite native branch exploration, not continuous/global approval','producer':'Fortran','units':{'n_b':'fm^-3','epsilon':'MeV fm^-3','pressure':'MeV fm^-3','mu_b':'MeV'},'models':{},'derivatives':{str(n):derivative_summary(n) for n in (64,128)}}
 samples=[]
 for m,v in plan['models'].items():
  rows=[]
  for src in (ROOT/v['reuse_selected_path'],BASE/m/'output/generated/selected.csv'):
   rows += [{**r,'source':str(src.relative_to(ROOT))} for r in read(src)]
  sel={mesh:{float(r['n_b']):r for r in rows if int(r['mesh'])==mesh} for mesh in (32,64)}
  extra=BASE/'ambiguous_N128'/m/'output/generated/selected.csv';n128={float(r['n_b']):r for r in read(extra)} if extra.exists() else {}
  intervals=[]
  for k in v['remaining_intervals']:
   left,right=map(float,(k['left'],k['right']));interval={**k,'meshes':{}}
   for mesh in (32,64):
    a,b=sel[mesh][left],sel[mesh][right]
    interval['meshes'][str(mesh)]={'delta_pressure':float(b['pressure'])-float(a['pressure']),'geometry':[a['geometry'],b['geometry']],'epsilon':[float(a['epsilon']),float(b['epsilon'])],'pressure':[float(a['pressure']),float(b['pressure'])]}
   if left in n128 and right in n128:
    keys=sorted(n for n in n128 if left<=n<=right);rs=[n128[n] for n in keys]
    interval['meshes']['128']={'density':keys,'delta_pressure':[float(b['pressure'])-float(a['pressure']) for a,b in zip(rs,rs[1:])],'geometry':[r['geometry'] for r in rs],'epsilon':[float(r['epsilon']) for r in rs],'pressure':[float(r['pressure']) for r in rs]}
   d64=interval['meshes']['64']['delta_pressure']
   interval['classification']='N16_DROP_REMOVED_AT_N64_SAMPLED_ENDPOINTS' if d64>=0 else 'HOLD_N64_PRESSURE_DROP'
   if '128' in interval['meshes']:
    interval['classification']='N64_DROP_REMOVED_AT_N128_THREE_SAMPLES' if min(interval['meshes']['128']['delta_pressure'])>=0 else 'HOLD_N128_PRESSURE_DROP'
   intervals.append(interval)
  attempts=read(BASE/m/'output/generated/attempts.csv')
  statuses=dict(collections.Counter(r['status'] for r in attempts))
  diff=[];changes=[]
  for n in map(float,v['endpoint_union']):
   a,b=sel[32][n],sel[64][n];relative=abs(float(a['pressure'])-float(b['pressure']))/max(abs(float(b['pressure'])),1e-300);diff.append(relative)
   if a['geometry']!=b['geometry']:changes.append({'n_b':n,'N32':a['geometry'],'N64':b['geometry']})
   for mesh,r in ((32,a),(64,b)):
    samples.append({'model':m,'n_b':r['n_b'],'mesh':mesh,'geometry':r['geometry'],'epsilon':r['epsilon'],'pressure':r['pressure'],'mu_b':r['mu_b'],'KKT':r['KKT'],'radius_stationarity':r['radius_stationarity'],'source_selected':r['source'],'source_node':r['node'],'source_attempt':r['attempt']})
  candidate_file=BASE/'COMMON_TANGENT_SUMMARY.json'
  candidates=json.loads(candidate_file.read_text())['models'] if candidate_file.exists() else {}
  for interval in intervals:
   interval['raw_source_support_status']='CROSS_PHASE_UNSUPPORTED' if interval['classification']=='HOLD_N128_PRESSURE_DROP' else 'SAMPLED_ORDERING_ONLY'
   interval['repair_outcome']='SAMPLED_PRESSURE_ORDERING_RECOVERED' if not interval['classification'].startswith('HOLD') else 'UNRESOLVED_PRESSURE_DROP'
   if interval['classification']=='HOLD_N128_PRESSURE_DROP' and m in candidates:
    interval['repair_outcome']='FINITE_TWO_BRANCH_COEXISTENCE_CANDIDATE'
    interval['candidate_numerical_target_met']=candidates[m]['numerical_target_met']
    interval['candidate_evidence']='COMMON_TANGENT_SUMMARY.json#'+m
    interval['production_cross_phase_supported']=False
  report['models'][m]={'remaining_intervals':len(intervals),'endpoint_union_count':len(v['endpoint_union']),'new_mesh_density_nodes':v['new_mesh_density_nodes'],'new_attempt_status_counts':statuses,'N64_negative_intervals':sum(x['meshes']['64']['delta_pressure']<0 for x in intervals),'max_N32_N64_relative_pressure_difference':max(diff),'N32_N64_geometry_changes':changes,'intervals':intervals}
 report['unique_new_N32_N64_mesh_density_nodes']=sum(v['new_mesh_density_nodes'] for v in plan['models'].values())
 report['reused_endpoint_mesh_nodes']=sum(2*len(v['reused_endpoint_values']) for v in plan['models'].values())
 report['requested_finite_campaign_status']='COMPLETED'
 report['repair_outcomes']=dict(collections.Counter(x['repair_outcome'] for v in report['models'].values() for x in v['intervals']))
 report['interval_classifications']=dict(collections.Counter(x['classification'] for v in report['models'].values() for x in v['intervals']))
 write(BASE/'ANALYSIS.json',report)
 with (BASE/'refined_endpoint_samples.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(samples[0]));writer.writeheader();writer.writerows(samples)
 print(json.dumps({'new_mesh_density_nodes':report['unique_new_N32_N64_mesh_density_nodes'],'classifications':report['interval_classifications'],'derivative_P_residual':report['derivatives']['64']['max_abs_P_residual_half_h']}))
if __name__=='__main__':main()
