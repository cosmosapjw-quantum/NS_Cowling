"""Create a separate working directory without overwriting another calculation."""
from __future__ import annotations
import argparse,hashlib,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PREREQUISITES=[
 'src/numerics/ns_bisection.f90','src/numerics/ns_rk4.f90','src/numerics/ns_adaptive_rk.f90',
 'src/eos/ns_polytrope.f90','src/structure/ns_tov.f90','src/structure/ns_target_mass.f90',
 'src/modes/ns_scan.f90','src/modes/ns_shooting_toy.f90']
IGNORE=shutil.ignore_patterns('solutions','build','.git','__pycache__','.pytest_cache','.runtime','*.pyc','*.o','*.mod','.runs','results')
WORKSPACE_README="""# 중성자별 구조와 비방사 진동

이 디렉터리는 직접 작성하는 `learner/`와, 제공된 `calculation/`, `methods/`, `support/`를 포함한다. 실행은 `START_HERE.md`에서 시작한다. 강의 노트와 회차별 계획, 과제별 PDF는 `documents/`에 있다.

제공된 선행 함수는 `WORKSPACE.json`에 명시되어 있다. 직접 작성한 함수와 제공된 함수를 구분한다. 자동으로 다른 풀이를 불러오지 않으며, 아직 작성하지 않은 함수는 `not_implemented`로 응답한다. 이미 함수를 완성한 뒤 복사한 작업공간에서는 작성한 소스가 그대로 유지된다.
"""
WORKSPACE_START="""# 작업 디렉터리에서 시작하기

이 디렉터리에는 실행할 `learner/`와 `support/`가 이미 들어 있다. 작업공간을 다시 만들 필요가 없다. GNU Fortran, fpm, Python을 준비하고 아래 명령으로 시작한다. `fpm`이 PATH에 없다면 `NS_FPM`에 실행파일 경로를 지정한다.

```bash
python -m venv ../ns-env
source ../ns-env/bin/activate
python -m pip install -r requirements.txt
python tools/run_assignment.py A00 --run --workspace learner --output ../A00.json
python tools/run_assignment.py A01 --run --workspace learner --output ../A01.json
python tools/run_assignment.py FIRST_TASK --plan
```

TOV에서 시작한 독자의 첫 직접 구현은 A09다. `learner/src/modes/ns_center_state.f90`, `ns_cowling.f90`, `ns_surface_residual.f90`를 강의 노트 제5장과 함께 작성한다. 기초 수치법부터 시작하는 독자는 A02부터 순서대로 진행한다. 이미 완성한 과제는 기록을 확인하고 다음 과제로 넘어간다.

```bash
python tools/run_assignment.py FIRST_TASK --run --workspace learner --output ../FIRST_TASK.json
```

A10까지 완성한 뒤에는 아래 명령으로 자신의 코드에서 별과 두 모드를 계산한다. 미구현 상태에서는 이 명령이 실패하는 것이 정상이다.

```bash
python tools/example_star.py --workspace learner --output ../my_star
```

강의 노트 `documents/Lecture_Notes.pdf`, 전체 계획 `documents/Course_Plan.pdf`, 큐시트 `documents/Assignments.pdf`, 빌드 파일 설명 `documents/FPM_Note.pdf`를 함께 읽는다. 제9·10장의 A20·A21은 제공된 calculation 경로에서 바로 수행할 수 있다. 직접 작성한 물리 함수를 사용하려면 별도로 명시하여 연결한다. 제공 adapter와 모델 생성기는 `curriculum/precision_adapters.json`에 구분되어 있다. 방법의 유도는 본문 제3·4·5·7·9장, 전체 수렴은 제10장에 있다. `python tools/stellar_methods.py --output ../my_study`로 제공 계산을 실행하고 함수 사용 안내는 `documents/Precision_Study.pdf`에서 확인한다. 모든 출력은 새 경로를 사용한다. 전체 자료의 `solutions/`를 사용하는 설명은 별도로 제공된 풀이 실행에 해당하며 이 작업 디렉터리의 실행 명령이 아니다.
"""
def create(output:Path,start:str='foundations')->dict:
 output=output.resolve()
 if output.exists():raise FileExistsError(f'Destination already exists: {output}')
 if output.is_relative_to(ROOT):raise ValueError('Choose a destination outside the source directory')
 if start not in {'foundations','cowling'}:raise ValueError('Unknown starting point')
 if start=='cowling':
  for rel in PREREQUISITES:
   if not (ROOT/'solutions'/rel).is_file():raise FileNotFoundError(f'Prerequisite unavailable: {rel}')
 shutil.copytree(ROOT,output,ignore=IGNORE)
 supplied=[]
 prior=ROOT/'WORKSPACE.json'
 if prior.is_file() and start=='foundations':
  supplied=json.loads(prior.read_text()).get('supplied_prerequisites',[])
 if start=='cowling':
  for rel in PREREQUISITES:
   src=ROOT/'solutions'/rel;dst=output/'learner'/rel
   shutil.copy2(src,dst)
   supplied.append({'file':rel,'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
 result={'start':start,'supplied_prerequisites':supplied,'exercise_workspace':'learner',
         'unmodified_exercises':['ns_center_state.f90','ns_cowling.f90','ns_surface_residual.f90'],
         'supplied_precision_adapters':json.loads((ROOT/'curriculum/precision_adapters.json').read_text())['files']}
 (output/'README.md').write_text(WORKSPACE_README)
 (output/'START_HERE.md').write_text(WORKSPACE_START.replace('FIRST_TASK', 'A09' if start=='cowling' or supplied else 'A02'))
 (output/'WORKSPACE.json').write_text(json.dumps(result,indent=2)+'\n')
 return result
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--start',choices=['foundations','cowling'],default='cowling')
 a=ap.parse_args()
 try:print(json.dumps(create(a.output,a.start),indent=2));return 0
 except (OSError,ValueError) as e:print(str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
