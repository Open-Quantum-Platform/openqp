# codex/casscf-analytic-gradient-review-v2-20260816

SHA: `3bdfb60dbd9c94ffdb0d8c487d8fcdd9a2ebcc78`  
판정: **동일 끝점 PR 병합 확인**  
분야: Wavefunction methods (검색용 분류)  
마지막 commit: 2026-08-17T07:23:10+09:00 / Cheol Ho Choi  
제목: Address CASSCF gradient review findings

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / codex/casscf-analytic-gradient-review-v2-20260816](https://github.com/karmachoi/openqp/tree/3bdfb60dbd9c94ffdb0d8c487d8fcdd9a2ebcc78)
- [gitlab / codex/casscf-analytic-gradient-review-v2-20260816](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/3bdfb60dbd9c94ffdb0d8c487d8fcdd9a2ebcc78)
- `local-15 / codex/casscf-analytic-gradient-review-v2-20260816`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `a22661e77ec58c8d0e73a3496beaf020cddc8fc8`
- 앞선 커밋 3 / 뒤처진 커밋 465
- non-merge patch: main과 일치 0, 다름 2
- merge/empty 등 patch 비교 제외: 1
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 26; 그중 현재 main과 동일 5, 다름 21

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #346: Analytic state-specific CASSCF nuclear gradient](https://github.com/Open-Quantum-Platform/openqp/pull/346) — state=closed; merged=2026-08-16T22:55:30Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `3bdfb60dbd9c` | 2026-08-17 | 다름 | commit 보존 또는 비교 제외 | Address CASSCF gradient review findings |
| `5b73b9516c38` | 2026-08-16 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge upstream main after MP2 gradient integration |
| `6dd537ee5d94` | 2026-08-16 | 다름 | commit 보존 또는 비교 제외 | Complete analytic state-specific CASSCF gradients |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| A | `docs/casscf_analytic_gradient.md` | 다름 |
| A | `examples/WF_methods/H2O_CASSCF_CAS44_grad.inp` | 동일 |
| A | `examples/WF_methods/H2O_CASSCF_CAS44_grad.json` | 다름 |
| A | `examples/WF_methods/H2O_CASSCF_CAS44_grad.oqp` | 다름 |
| A | `examples/WF_methods/H4_CASSCF_CAS22_ROOT1_grad.inp` | 동일 |
| A | `examples/WF_methods/H4_CASSCF_CAS22_ROOT1_grad.json` | 다름 |
| A | `examples/WF_methods/H4_CASSCF_CAS22_ROOT1_grad.oqp` | 다름 |
| M | `examples/WF_methods/LiH_CASSCF_grad.inp` | 동일 |
| M | `examples/WF_methods/LiH_CASSCF_grad.oqp` | 다름 |
| A | `examples/WF_methods/LiH_CASSCF_optimize.inp` | 동일 |
| A | `examples/WF_methods/LiH_CASSCF_optimize.json` | 다름 |
| A | `examples/WF_methods/LiH_CASSCF_optimize.oqp` | 다름 |
| M | `examples/WF_methods/README.md` | 다름 |
| M | `include/oqp.h` | 다름 |
| M | `pyoqp/oqp/library/casscf.py` | 다름 |
| A | `pyoqp/oqp/library/casscf_gradient.py` | 다름 |
| M | `pyoqp/oqp/library/single_point.py` | 다름 |
| M | `pyoqp/oqp/library/wf_numgrad.py` | 다름 |
| M | `pyoqp/oqp/openqp.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `source/modules/casscf_driver.F90` | 다름 |
| A | `source/modules/casscf_gradient.F90` | 동일 |
| A | `tests/test_casscf_gradient.py` | 다름 |
| M | `tests/test_casscf_numgrad.py` | 다름 |
| M | `tests/test_openqp_api.py` | 다름 |
| M | `tests/test_oqp_input.py` | 다름 |

## 이 끝점을 사용 중인 로컬 worktree

- `/Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-pr344-review-maintenance` — exists=True, status 항목 0; 
