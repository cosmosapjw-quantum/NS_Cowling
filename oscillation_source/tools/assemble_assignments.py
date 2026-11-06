#!/usr/bin/env python3
from pathlib import Path
import json,shutil,sys
import fitz
R=Path(__file__).resolve().parents[1]
import argparse
ap=argparse.ArgumentParser(description='Build the assignment index, combine cue sheets and retain standalone LaTeX.')
ap.add_argument('--work',type=Path,default=R.parent/'document_build')
W=ap.parse_args().work.resolve()
sys.path.insert(0,str(R));from tools.build_documents import front,run
rows=json.loads((R/'curriculum/assignments.json').read_text())['assignments']
text=front('중성자별 구조와 비방사 진동','과제 큐시트 모음')
text+='## 본문 순서로 읽기\n\n과제 번호는 식별자다. EOS 구간 계산 A20은 제9장, 전체 수렴 A21은 제10장의 본 과제이며 아래 순서로 읽는다. 제공 방법을 사용하는 과제와 직접 물리 함수를 작성하는 과제는 각 카드에서 구별한다.\n\n| 과제 | 제목 | 합본 시작 쪽 |\n|:--|:--|--:|\n'
page=2
for row in rows:
 text+=f"| {row['id']} | {row['title']} | {page} |\n"
 page+=len(fitz.open(R/'documents/cues'/(row['id']+'.pdf')))
(W/'Assignments_Index.md').write_text(text)
cmd=['pandoc','Assignments_Index.md','--from=markdown+tex_math_single_backslash','--to=latex','--standalone','--no-highlight','--lua-filter=paths.lua','--include-in-header=header.tex','-V','documentclass=article','-V','papersize=a4','-V','fontsize=10pt','-V','geometry:margin=22mm','-V','mainfont=Liberation Serif','-V','sansfont=Liberation Sans','-V','monofont=DejaVu Sans Mono','-o','Assignments_Index.tex']
run(cmd,W,W/'logs/Assignments_Index.pandoc.log')
for i in [1,2]:run(['xelatex','-interaction=nonstopmode','-halt-on-error','Assignments_Index.tex'],W,W/f'logs/Assignments_Index.xelatex{i}.log')
index=fitz.open(W/'Assignments_Index.pdf');assert len(index)==1,'Index grew; fix page references before combining'
combined=fitz.open();combined.insert_pdf(index);toc=[[1,'과제 목록',1]]
for row in rows:
 page=len(combined)+1;combined.insert_pdf(fitz.open(R/'documents/cues'/(row['id']+'.pdf')));toc.append([1,row['id']+' '+row['title'],page])
combined.set_toc(toc);combined.save(R/'documents/Assignments.pdf',garbage=4,deflate=True)
# Current, standalone LaTeX source only, with its actual figure assets.
source=R/'documents/latex'
if source.exists():shutil.rmtree(source)
source.mkdir()
for pat in ['*.tex','*.md','*.lua']:
 for p in W.glob(pat):shutil.copy2(p,source/p.name)
shutil.copytree(W/'assets',source/'assets')
(source/'BUILD.md').write_text('''# PDF를 다시 만들기\n\n이 디렉터리에서 `xelatex -interaction=nonstopmode -halt-on-error Lecture_Notes.tex`를 목차와 참조가 안정될 때까지 최소 세 번 실행한다. 다른 문서도 같은 방식으로 빌드한다. Pandoc을 이용한 원문 재구성은 패키지 루트의 `python tools/build_documents.py --limit 40`을 사용하고, `python tools/assemble_assignments.py`로 합본을 만든다. 합본 도구에는 requirements-docs.txt의 PyMuPDF가 필요하다. Noto CJK KR, Liberation, DejaVu 계열 글꼴과 XeLaTeX가 필요하다. 글꼴 파일은 포함하지 않는다.\n\n큐시트 합본은 Assignments_Index.pdf 다음에 A00부터 본문 순서의 각 과제 PDF를 결합한다. 순서는 curriculum/assignments.json을 따른다. 개별 큐시트는 독립된 두 쪽 문서다.\n''')
print('Assignments',len(combined),'pages; standalone sources',len(list(source.glob('*.tex'))))
