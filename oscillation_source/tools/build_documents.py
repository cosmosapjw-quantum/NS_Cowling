#!/usr/bin/env python3
"""Build the lecture notes, course plan and assignment cue sheets with XeLaTeX."""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,re,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
CHAPTERS={
1:['ch01','ch01_bootstrap'],2:['ch02_units','ch02_numerics'],
3:['ch03_eos','ch03_tov','ch03_connection'],4:['ch04'],5:['ch05'],
6:['ch06','ch06_reasoning'],7:['ch07','ch07_bridge'],8:['ch08'],
9:['ch09','ch09_interpretation'],10:['ch10','ch10_research'],11:['ch11'],12:['ch12']}
TITLES=['평형 별에서 진동하는 별로','단위와 수치법','열역학 EOS와 TOV 배경',
'해석적 고유값 문제','상대론적 Cowling의 유도와 계산','모드 식별과 수렴',
'EOS 접합과 물질적 계면','핵물질·조성·액적 모형','표의 열역학적 표현',
'통제 비교와 주파수 차이','자신의 코드로 실행하기','누적 음향시간']
HEADER=r'''\usepackage{xeCJK}
\setCJKmainfont{Noto Serif CJK KR}
\setCJKsansfont{Noto Sans CJK KR}
\setCJKmonofont{Noto Sans Mono CJK KR}
\xeCJKsetup{CJKspace=true}
\usepackage{xurl}
\usepackage{fvextra}
\usepackage{fancyhdr}
\usepackage{needspace}
\AtBeginDocument{\DefineVerbatimEnvironment{verbatim}{Verbatim}{breaklines=true,breakanywhere=true,fontsize=\small}}
\setlength{\emergencystretch}{4em}
\setlength{\headheight}{15pt}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small 중성자별 구조와 비방사 진동}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\setcounter{tocdepth}{2}
'''
LUA=r'''function Code(el)
 if FORMAT:match('latex') and #el.text > 18 and not el.text:find('[\128-\255]') and not el.text:find('|',1,true) then
 return pandoc.RawInline('latex','\\path|' .. el.text .. '|') end
end'''
def front(title,subtitle=''):
 return '---\ntitle: '+json.dumps(title,ensure_ascii=False)+'\nsubtitle: '+json.dumps(subtitle,ensure_ascii=False)+'\ndate: ""\nlang: ko-KR\n---\n\n'
def demote(s):return re.sub(r'^(#{1,5}) ',lambda m:m[1]+'# ',s,flags=re.M)
def lecture():
 text=front('중성자별 구조와 비방사 진동','평형 배경에서 Cowling 고유값과 통제 비교까지')
 for i,names in CHAPTERS.items():
  text+=('\n\\clearpage\n\n' if i==1 else '\n\\Needspace{12\\baselineskip}\n\n')+'# '+TITLES[i-1]+'\n\n'
  for j,n in enumerate(names):
   s=(ROOT/'notes/chapters'/f'{n}.md').read_text()
   s=re.sub(r'^(#{1,6}) (?:제\d+장\.?\s*|\d+(?:\.\d+)*\.?\s+)',r'\1 ',s,flags=re.M)
   if j==0:s=re.sub(r'^# [^\n]+\n','',s,count=1)
   else:s=demote(s)
   if i==8 and j==0:
    # Teach uniform matter and its composition before moving to clustered matter.
    marker='## 군집화된 저밀도 물질의 최소 liquid-drop toy'
    if marker not in s:raise ValueError('Missing nuclear/cluster section boundary')
    core,cluster=s.split(marker,1)
    text+=core+'\n\n'+demote((ROOT/'notes/worked/nuclear_calibration.md').read_text())+'\n\n'+marker+cluster+'\n\n'
   else:text+=s+'\n\n'
 text+='\n\\clearpage\n\n'+(ROOT/'notes/references.md').read_text()
 return text.replace('(../figures/','(assets/').replace('(notes/figures/','(assets/')
def jobs():
 a=json.loads((ROOT/'curriculum/assignments.json').read_text())['assignments']
 result={'Lecture_Notes':(lecture(),False),
 'Course_Plan':(front('중성자별 구조와 비방사 진동: 강의계획','TOV 구현 경험자를 위한 16회 과정')+(ROOT/'curriculum/course_plan.md').read_text()+'\n\n\\begingroup\\small\n'+(ROOT/'notes/references.md').read_text()+'\n\\endgroup\n',False),
 'Precision_Study':(front('제공 수치 방법의 사용 안내','함수 호출과 수렴 결과 읽기')+(ROOT/'notes/precision_study.md').read_text(),False),
 'Mentor_Guide':(front('설명에서 독립적인 판단으로','지도와 피드백을 위한 안내')+(ROOT/'curriculum/instructor/mentor_guide.md').read_text(),False),
 'FPM_Note':(front('fpm.toml의 역할')+(ROOT/'docs/FPM_NOTE.md').read_text(),True),
 'Code_Comparison':(front('독립 계산과 해석적 검산')+(ROOT/'notes/comparison.md').read_text(),False)}
 for row in a:
  text=(ROOT/row['cue_source']).read_text()
  result[row['id']]=(front(row['id']+' · '+row['title'],'과제 큐시트')+re.sub(r'^# [^\n]+\n','',text,count=1),True)
 return result

def run(cmd,cwd,log,timeout=150):
 with log.open('w') as f:
  f.write('$ '+' '.join(map(str,cmd))+'\n');f.flush()
  p=subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
 if p.returncode:raise RuntimeError(str(log)+': exit '+str(p.returncode))

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--only');ap.add_argument('--limit',type=int,default=4);ap.add_argument('--work',type=Path,default=ROOT.parent/'document_build');ap.add_argument('--force',action='store_true');args=ap.parse_args()
 w=args.work.resolve();w.mkdir(parents=True,exist_ok=True);logs=w/'logs';logs.mkdir(exist_ok=True)
 out=ROOT/'documents';out.mkdir(exist_ok=True);(out/'cues').mkdir(exist_ok=True)
 (w/'header.tex').write_text(HEADER);(w/'paths.lua').write_text(LUA);shutil.copytree(ROOT/'notes/figures',w/'assets',dirs_exist_ok=True)
 book=jobs();selected={args.only:book[args.only]} if args.only else book
 receipt=w/'results.json';rec=json.loads(receipt.read_text()) if receipt.exists() else {};count=0
 for name,(text,cue) in selected.items():
  header=HEADER if name=='Lecture_Notes' else HEADER.replace('\\usepackage{needspace}\n','')
  sig=hashlib.sha256((text+header+LUA+('cue' if cue else 'stable-toc')).encode()).hexdigest()
  target=out/('cues/'+name+'.pdf' if re.fullmatch('[AP][0-9]{1,2}',name) else name+'.pdf')
  if not args.force and name in rec and rec[name]['input_sha256']==sig and target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==rec[name]['sha256']:continue
  if count>=args.limit:break
  (w/'header.tex').write_text(header)
  (w/(name+'.md')).write_text(text)
  cmd=['pandoc',name+'.md','--from=markdown+tex_math_single_backslash','--to=latex','--standalone','--no-highlight','--lua-filter=paths.lua','--include-in-header=header.tex','-V','documentclass=article','-V','papersize=a4','-V','fontsize=10pt','-V','geometry:margin=22mm','-V','mainfont=Liberation Serif','-V','sansfont=Liberation Sans','-V','monofont=DejaVu Sans Mono','-V','colorlinks=false','-V','linestretch=1.12','-o',name+'.tex']
  if not cue:cmd+=['--toc','--number-sections']
  run(cmd,w,logs/(name+'.pandoc.log'))
  previous=None
  for i in range(1,6):
   run(['xelatex','-interaction=nonstopmode','-halt-on-error',name+'.tex'],w,logs/(name+f'.xelatex{i}.log'))
   current=tuple((w/(name+ext)).read_bytes() if (w/(name+ext)).exists() else b'' for ext in ('.aux','.toc','.out'))
   if i>=(2 if cue else 3) and current==previous:break
   previous=current
  else:raise RuntimeError(name+': cross-references did not stabilize')
  shutil.copy2(w/(name+'.pdf'),target)
  warnings=[l for l in (w/(name+'.log')).read_text(errors='replace').splitlines() if 'Overfull' in l or 'Missing character:' in l]
  rec[name]={'input_sha256':sig,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'size':target.stat().st_size,'warnings':warnings,'pdf':str(target.relative_to(ROOT)),'latex_passes':i}
  receipt.write_text(json.dumps(rec,ensure_ascii=False,indent=2)+'\n');count+=1;print(name,'built',len(warnings),'warnings',flush=True)
 print('built this call',count,'total',len(rec),flush=True)
if __name__=='__main__':main()
