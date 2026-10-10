# feat/builtin-optimizer

SHA: `4c0cdbcbeaf169b6035a6d60651f98facdfd34e9`  
판정: **동일 끝점 PR 병합 확인**  
분야: Geometry optimization (검색용 분류)  
마지막 commit: 2026-06-07T14:10:11+09:00 / Cheol Ho Choi  
제목: Merge branch 'main' into feat/builtin-optimizer

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / feat/builtin-optimizer](https://github.com/karmachoi/openqp/tree/4c0cdbcbeaf169b6035a6d60651f98facdfd34e9)
- [gitlab / feat/builtin-optimizer](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/4c0cdbcbeaf169b6035a6d60651f98facdfd34e9)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `10c7e86e0d61a43c4a6e0f2aa43a4c311c726efc`
- 앞선 커밋 16 / 뒤처진 커밋 729
- non-merge patch: main과 일치 0, 다름 13
- merge/empty 등 patch 비교 제외: 3
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 17; 그중 현재 main과 동일 4, 다름 13

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #186: OQP-native NumPy/SciPy geometry optimizer (lib=oqp): minima, TS, MECI/MECP, three-state CI, NEB, IRC, MEP; TRIC/DLC coordinates](https://github.com/Open-Quantum-Platform/openqp/pull/186) — state=closed; merged=2026-06-07T05:37:51Z; merge commit in main=True
- [personal PR #8: OQP-native NumPy/SciPy geometry optimizer (lib=oqp): minima, TS, MECI/MECP, three-state CI, NEB, IRC, MEP; TRIC/DLC coords; default backend](https://github.com/karmachoi/openqp/pull/8) — state=closed; merged=아니오; merge commit in main=False

## 남은 non-merge patch 집합이 같은 다른 끝점

- [feat/builtin-optimizer](cf1f3949163c.md)

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `4c0cdbcbeaf1` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into feat/builtin-optimizer |
| `cf1f3949163c` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Address PR review (P1): judge climbing-image NEB convergence by its climbing force |
| `fc0841be6971` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'karmachoi/feat/builtin-optimizer' into feat/builtin-optimizer |
| `534ad7911551` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Address PR review: reject mismatched NEB product atoms; resolve product path via input dir |
| `4c2b62e75a70` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into feat/builtin-optimizer |
| `60983b622fd1` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Log the requested coordsys by name (TRIC/DLC/RIC/CART) with fallback flag |
| `8c14ee5ee741` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Rename the in-house optimizer builtin -> oqp |
| `8f419d02e31e` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Make builtin the default optimizer (geometric still available) |
| `9c364373ddaa` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Add TRIC and DLC coordinate systems (TRIC is the default) |
| `ab47fdf28301` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Rename native optimizer -> builtin; add MEP (lib=builtin, runtype=mep) |
| `c8e927954d80` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Harden native NEB: endpoint pre-optimization + exposed FIRE/climb controls |
| `ab54a2095c51` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Add native IRC (Gonzalez-Schlegel) -- lib=native, runtype=irc |
| `f8a4466702b4` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Add native NEB reaction-path optimizer (lib=native, runtype=neb) |
| `decd28f42f19` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Add three-state conical intersection (TCI) search via native optimizer |
| `0cf91279cd98` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Wire native MECI/MECP (penalty objective) into lib=native |
| `f4f0b3abf6b4` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Add native NumPy/SciPy geometry optimizer (lib=native) |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| A | `examples/OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.inp` | 다름 |
| A | `examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.inp` | 동일 |
| A | `examples/OPT/HCN_RHF-DFT_IRC_OQP.inp` | 동일 |
| A | `examples/OPT/HCN_RHF-DFT_NEB_OQP.inp` | 동일 |
| A | `examples/OPT/HCN_RHF-DFT_NEB_OQP_product.xyz` | 동일 |
| A | `pyoqp/oqp/library/liboqp.py` | 다름 |
| M | `pyoqp/oqp/library/libscipy.py` | 다름 |
| A | `pyoqp/oqp/library/oqp_coords.py` | 다름 |
| A | `pyoqp/oqp/library/oqp_engine.py` | 다름 |
| A | `pyoqp/oqp/library/oqp_irc.py` | 다름 |
| A | `pyoqp/oqp/library/oqp_neb.py` | 다름 |
| M | `pyoqp/oqp/library/runfunc.py` | 다름 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/pyoqp.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `tests/test_geometric_optimizer.py` | 다름 |
| A | `tests/test_oqp_optimizer.py` | 다름 |
