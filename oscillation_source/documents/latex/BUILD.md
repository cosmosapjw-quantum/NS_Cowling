# PDF를 다시 만들기

이 디렉터리에서 `xelatex -interaction=nonstopmode -halt-on-error Lecture_Notes.tex`를 목차와 참조가 안정될 때까지 최소 세 번 실행한다. 다른 문서도 같은 방식으로 빌드한다. Pandoc을 이용한 원문 재구성은 패키지 루트의 `python tools/build_documents.py --limit 40`을 사용하고, `python tools/assemble_assignments.py`로 합본을 만든다. 합본 도구에는 requirements-docs.txt의 PyMuPDF가 필요하다. Noto CJK KR, Liberation, DejaVu 계열 글꼴과 XeLaTeX가 필요하다. 글꼴 파일은 포함하지 않는다.

큐시트 합본은 Assignments_Index.pdf 다음에 A00부터 본문 순서의 각 과제 PDF를 결합한다. 순서는 curriculum/assignments.json을 따른다. 개별 큐시트는 독립된 두 쪽 문서다.
