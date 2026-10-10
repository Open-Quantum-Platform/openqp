# codex/casscf-numgrad-readme-20260815

SHA: `ad2fe53a0d2c662852287b5cd1d649787072b9f1`  
판정: **동일 끝점 PR 병합 확인**  
분야: Wavefunction methods (검색용 분류)  
마지막 commit: 2026-08-16T06:59:48+09:00 / Cheol Ho Choi  
제목: Address CASSCF gradient review feedback

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / codex/casscf-numgrad-readme-20260815](https://github.com/karmachoi/openqp/tree/ad2fe53a0d2c662852287b5cd1d649787072b9f1)
- [gitlab / codex/casscf-numgrad-readme-20260815](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/ad2fe53a0d2c662852287b5cd1d649787072b9f1)
- `local-15 / codex/casscf-numgrad-readme-20260815`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `71f672f4f84bfdb30fb9d883dc0a1b28c12e4395`
- 앞선 커밋 3 / 뒤처진 커밋 467
- non-merge patch: main과 일치 0, 다름 3
- merge/empty 등 patch 비교 제외: 0
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 19; 그중 현재 main과 동일 1, 다름 18

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #343: Add numerical nuclear gradients for CASSCF and SA-CASSCF](https://github.com/Open-Quantum-Platform/openqp/pull/343) — state=closed; merged=2026-08-16T01:21:27Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `ad2fe53a0d2c` | 2026-08-16 | 다름 | commit 보존 또는 비교 제외 | Address CASSCF gradient review feedback |
| `4579fb5c641f` | 2026-08-15 | 다름 | commit 보존 또는 비교 제외 | Complete CASSCF numerical gradient interface |
| `46183430562e` | 2026-08-15 | 다름 | commit 보존 또는 비교 제외 | Add numerical CASSCF nuclear gradients |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `README.md` | 다름 |
| A | `examples/WF_methods/LiH_CASSCF_grad.inp` | 다름 |
| A | `examples/WF_methods/LiH_CASSCF_grad.json` | 다름 |
| A | `examples/WF_methods/LiH_CASSCF_grad.oqp` | 다름 |
| A | `examples/WF_methods/LiH_SA-CASSCF_grad.inp` | 다름 |
| A | `examples/WF_methods/LiH_SA-CASSCF_grad.json` | 다름 |
| A | `examples/WF_methods/LiH_SA-CASSCF_grad.oqp` | 다름 |
| M | `examples/WF_methods/README.md` | 다름 |
| M | `pyoqp/oqp/library/pt2_numgrad.py` | 동일 |
| M | `pyoqp/oqp/library/single_point.py` | 다름 |
| A | `pyoqp/oqp/library/wf_numgrad.py` | 다름 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/openqp.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `pyoqp/oqp/utils/oqp_input.py` | 다름 |
| A | `tests/test_casscf_numgrad.py` | 다름 |
| M | `tests/test_openqp_api.py` | 다름 |
| M | `tests/test_oqp_input.py` | 다름 |
| M | `tests/test_oqp_input_schema_manifest.py` | 다름 |

## 이 끝점을 사용 중인 로컬 worktree

- `/Users/cheolhochoi/Documents/Codex/2026-08-15/openqp-casscf-numgrad-readme` — exists=True, status 항목 0; 
