# feat/mp2

SHA: `e550dd0be8df8e1a5ba67ea2065bdd08a2808c22`  
판정: **고유 변경·반영 여부 검토**  
분야: Wavefunction methods (검색용 분류)  
마지막 commit: 2026-07-03T20:12:44+09:00 / Cheol Ho Choi  
제목: feat: add MP2 Python theory helper

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [gitlab / feat/mp2](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/e550dd0be8df8e1a5ba67ea2065bdd08a2808c22)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `f2a24c1ce97077e2b286e0d4baf121526cc40161`
- 앞선 커밋 6 / 뒤처진 커밋 659
- non-merge patch: main과 일치 1, 다름 3
- merge/empty 등 patch 비교 제외: 2
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 13; 그중 현재 main과 동일 0, 다름 13

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.


## 이 끝점을 포함하는 같은 분야의 후속 브랜치 후보

- [github-personal/pr252-mp2-pythonic, gitlab/pr252-mp2-pythonic](5d7b4cc01cef.md) — 추가 커밋 5

이 목록은 ancestry 관계이며 기능 대체나 최선의 통합 대상을 보증하지 않는다.

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `e550dd0be8df` | 2026-07-03 | 다름 | commit 보존 또는 비교 제외 | feat: add MP2 Python theory helper |
| `3c4bbc2a6298` | 2026-07-02 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'origin/main' into feat/mp2 |
| `1a2103d82237` | 2026-07-02 | 다름 | commit 보존 또는 비교 제외 | Address MP2 review issues |
| `da12c58bc9bc` | 2026-07-02 | 다름 | commit 보존 또는 비교 제외 | Add standalone MP2 ground-state method (method=mp2) |
| `df7bda425f62` | 2026-05-26 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'Open-Quantum-Platform:main' into main |
| `e6468d103837` | 2026-05-22 | 일치 | commit 보존 또는 비교 제외 | fix: printing electric moments |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.gitignore` | 다름 |
| A | `examples/MP2/README.md` | 다름 |
| A | `examples/MP2/h2o_ump2_6-31g.inp` | 다름 |
| A | `examples/MP2/h2o_ump2_6-31g.ref.txt` | 다름 |
| M | `include/oqp.h` | 다름 |
| M | `pyoqp/oqp/library/single_point.py` | 다름 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/openqp.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| A | `source/modules/mp2_energy.F90` | 다름 |
| A | `source/mp2_lib.F90` | 다름 |
| A | `tests/test_mp2_input_checker.py` | 다름 |
| M | `tests/test_openqp_api.py` | 다름 |
