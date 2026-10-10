# quantum-computing-upstream

SHA: `03a23cdffe4b8f27e1b848e83653c705cff46b9e`  
판정: **동일 끝점 PR 병합 확인**  
분야: Wavefunction methods (검색용 분류)  
마지막 commit: 2026-06-20T17:54:44+09:00 / Cheol Ho Choi  
제목: Merge branch 'main' into quantum-computing-upstream

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / quantum-computing-upstream](https://github.com/karmachoi/openqp/tree/03a23cdffe4b8f27e1b848e83653c705cff46b9e)
- [gitlab / quantum-computing-upstream](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/03a23cdffe4b8f27e1b848e83653c705cff46b9e)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `bac3655845911f7083e9cf6b4d50d1c85434ba28`
- 앞선 커밋 4 / 뒤처진 커밋 700
- non-merge patch: main과 일치 0, 다름 3
- merge/empty 등 patch 비교 제외: 1
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 13; 그중 현재 main과 동일 9, 다름 4

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #215: feat(quantum): add FCIDUMP export bridge](https://github.com/Open-Quantum-Platform/openqp/pull/215) — state=closed; merged=2026-06-20T12:34:13Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `03a23cdffe4b` | 2026-06-20 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into quantum-computing-upstream |
| `97012c357adc` | 2026-06-20 | 다름 | commit 보존 또는 비교 제외 | fix(quantum): make ERI export safe for MPI and reruns |
| `664fcb0eb0e0` | 2026-06-20 | 다름 | commit 보존 또는 비교 제외 | fix(quantum): keep upstream CFFI declarations current |
| `b6726332be5f` | 2026-06-20 | 다름 | commit 보존 또는 비교 제외 | Add oqp.quantum FCIDUMP export for quantum computing |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| A | `examples/QUANTUM/export_fcidump.py` | 동일 |
| A | `examples/QUANTUM/h2.inp` | 동일 |
| M | `include/oqp.h` | 다름 |
| M | `pyoqp/oqp/library/__init__.py` | 다름 |
| A | `pyoqp/oqp/library/ints_2e.py` | 동일 |
| A | `pyoqp/oqp/quantum/README.md` | 동일 |
| A | `pyoqp/oqp/quantum/__init__.py` | 동일 |
| A | `pyoqp/oqp/quantum/fcidump.py` | 동일 |
| A | `pyoqp/oqp/quantum/hamiltonian.py` | 동일 |
| A | `pyoqp/oqp/quantum/integrals.py` | 동일 |
| A | `source/modules/int2e.F90` | 다름 |
| M | `source/tagarray_driver.F90` | 다름 |
| A | `tests/test_quantum_fcidump.py` | 동일 |
