# draft/thermo-symmetry-number

SHA: `1cfc745df1ff37cd240b886c7e2f25ec6e37c86d`  
판정: **동일 끝점 PR 병합 확인**  
분야: Symmetry (검색용 분류)  
마지막 commit: 2026-08-12T08:29:15+09:00 / Cheol Ho Choi  
제목: Merge main after PR #303 into thermochemistry branch

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / draft/thermo-symmetry-number](https://github.com/karmachoi/openqp/tree/1cfc745df1ff37cd240b886c7e2f25ec6e37c86d)
- [gitlab / draft/thermo-symmetry-number](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/1cfc745df1ff37cd240b886c7e2f25ec6e37c86d)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `ee058e9f36c4c276152c9b90f890dca21d43b182`
- 앞선 커밋 14 / 뒤처진 커밋 479
- non-merge patch: main과 일치 1, 다름 12
- merge/empty 등 patch 비교 제외: 1
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 8; 그중 현재 main과 동일 3, 다름 5

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #320: Thermochemistry: apply the rotational symmetry number, fix linear rotors and the Gibbs sign](https://github.com/Open-Quantum-Platform/openqp/pull/320) — state=closed; merged=2026-08-11T23:57:31Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `1cfc745df1ff` | 2026-08-12 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge main after PR #303 into thermochemistry branch |
| `325890d301a1` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): harden the tolerance guard and the screen against overflow |
| `8103d91af4f3` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): make the tolerance guard and the near-linear test axis-free |
| `27301ed3b69e` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): refuse a tolerance that cannot be matched against |
| `01ae7c6cb8c9` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | fix(thermo): close the near-linear hole in the sigma screen, and make its test real |
| `fbf3cf5fb516` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | perf(thermo): skip the symmetry scan when nothing can be permuted |
| `89b683bc20c4` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | fix(thermo): keep sigma out of the positional argument list |
| `e43a0cc200f2` | 2026-08-06 | 다름 | commit 보존 또는 비교 제외 | fix(thermo): verify sigma is a group order, and share the moment mask with the printer |
| `1fe20f70e40c` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | test(thermo): contain the deliberate divide-by-zero in the nonlinear-branch test |
| `ba6fbfb07816` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | docs(thermo): record the two documented limits of the symmetry number |
| `2b72b6b3c093` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(thermo): select the vanishing linear moment on the inertia, not on rc/rt |
| `e9188a6678f8` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(thermo): G = H - TS, not H + TS |
| `c6532ac552b9` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(thermo): apply the rotational symmetry number, and make linear rotors finite |
| `c35898705363` | 2026-08-05 | 일치 | commit 보존 또는 비교 제외 | docs: [draft] Rotational symmetry number is missing from the thermochemistry |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| A | `docs/planned/thermo-sigma.md` | 다름 |
| M | `pyoqp/oqp/library/frequency.py` | 동일 |
| M | `pyoqp/oqp/library/single_point.py` | 다름 |
| M | `pyoqp/oqp/library/symmetry_detect.py` | 동일 |
| M | `pyoqp/oqp/utils/file_utils.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `tests/test_symmetry_parser_checker.py` | 다름 |
| A | `tests/test_thermochemistry_symmetry_number.py` | 동일 |
