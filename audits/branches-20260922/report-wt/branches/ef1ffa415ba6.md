# feat/advanced-guess

SHA: `ef1ffa415ba6ce28b7bf4e73638daa937616ee9a`  
판정: **동일 끝점 PR 병합 확인**  
분야: SCF / response (검색용 분류)  
마지막 commit: 2026-05-25T10:32:35+09:00 / Cheol Ho Choi  
제목: Merge branch 'main' into feat/advanced-guess

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / feat/advanced-guess](https://github.com/karmachoi/openqp/tree/ef1ffa415ba6ce28b7bf4e73638daa937616ee9a)
- [github-private / feat/advanced-guess](https://github.com/karmachoi/openqp-private/tree/ef1ffa415ba6ce28b7bf4e73638daa937616ee9a)
- [gitlab / feat/advanced-guess](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/ef1ffa415ba6ce28b7bf4e73638daa937616ee9a)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `60c42ba7186201b0cc633fc2769ae8783118a67e`
- 앞선 커밋 8 / 뒤처진 커밋 766
- non-merge patch: main과 일치 0, 다름 7
- merge/empty 등 patch 비교 제외: 1
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 21; 그중 현재 main과 동일 3, 다름 18

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #132: feat: add advanced OpenQP guess modes](https://github.com/Open-Quantum-Platform/openqp/pull/132) — state=closed; merged=2026-05-25T02:18:16Z; merge commit in main=True

## 남은 non-merge patch 집합이 같은 다른 끝점

- [feat/advanced-guess](a3ec89074053.md)

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `ef1ffa415ba6` | 2026-05-25 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into feat/advanced-guess |
| `a3ec89074053` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | chore: remove geomeTRIC changes from advanced guess PR |
| `9bae918f65d9` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | fix: use env flag for Docker credentials |
| `693c66a20f03` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | fix: skip Docker Hub push without secrets |
| `3217575ccd54` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | fix: stabilize Docker build workflow |
| `2b70824e843a` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | test: add advanced guess OpenQP examples |
| `0485d9ca3841` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: replace MOKIT guess export with native PySCF exporter |
| `e40bc47d9de4` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: add advanced OpenQP guess modes |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `README.md` | 다름 |
| A | `examples/other/h2o_rhf_3-21g_sad.inp` | 다름 |
| A | `examples/other/h2o_rhf_3-21g_sad.json` | 다름 |
| A | `examples/other/h2o_rhf_3-21g_sap.inp` | 동일 |
| A | `examples/other/h2o_rhf_3-21g_sap.json` | 다름 |
| M | `pyoqp/oqp/library/external.py` | 다름 |
| M | `pyoqp/oqp/library/guess.py` | 다름 |
| M | `pyoqp/oqp/library/libdlfind.py` | 다름 |
| M | `pyoqp/oqp/library/libscipy.py` | 다름 |
| M | `pyoqp/oqp/library/single_point.py` | 다름 |
| M | `pyoqp/oqp/periodic_table/__init__.py` | 동일 |
| M | `pyoqp/oqp/pyoqp.py` | 다름 |
| M | `pyoqp/oqp/utils/file_utils.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `pyoqp/setup.py` | 다름 |
| M | `pyproject.toml` | 다름 |
| A | `tests/test_advanced_guess.py` | 다름 |
| A | `tests/test_advanced_guess_examples.py` | 다름 |
| A | `tests/test_native_pyscf_exporter.py` | 다름 |
| A | `tests/test_periodic_table.py` | 동일 |
| A | `tests/test_single_point_scf_fallback.py` | 다름 |
