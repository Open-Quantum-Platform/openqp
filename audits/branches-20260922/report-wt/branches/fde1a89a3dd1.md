# feat/giao-nmr

SHA: `fde1a89a3dd1110a60af9557d1fde8606b8e90fd`  
판정: **고유 변경·반영 여부 검토**  
분야: NMR (검색용 분류)  
마지막 commit: 2026-06-03T07:15:46+09:00 / cheolhochoi  
제목: Document validated GIAO two-electron debug checkpoint

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-private / feat/giao-nmr](https://github.com/karmachoi/openqp-private/tree/fde1a89a3dd1110a60af9557d1fde8606b8e90fd)
- [gitlab / feat/giao-nmr](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/fde1a89a3dd1110a60af9557d1fde8606b8e90fd)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `6ab4b86701bf972b56acc7c79bd40338a51b13ac`
- 앞선 커밋 16 / 뒤처진 커밋 748
- non-merge patch: main과 일치 0, 다름 16
- merge/empty 등 patch 비교 제외: 0
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 33; 그중 현재 main과 동일 1, 다름 32

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.


## 이 끝점을 포함하는 같은 분야의 후속 브랜치 후보

- [github-private/claude/cool-ritchie-vaMB6, gitlab/github-private/claude/cool-ritchie-vaMB6, gitlab/github-private/pull/1/head](1f4c0da61efe.md) — 추가 커밋 18

이 목록은 ancestry 관계이며 기능 대체나 최선의 통합 대상을 보증하지 않는다.

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `fde1a89a3dd1` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | Document validated GIAO two-electron debug checkpoint |
| `71a88a878f4a` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | fix: advance native GIAO two-electron contraction |
| `04d6a39ee5bc` | 2026-06-02 | 다름 | commit 보존 또는 비교 제외 | test: harden GIAO two-electron oracle |
| `4c3bce280630` | 2026-06-02 | 다름 | commit 보존 또는 비교 제외 | chore(nmr): ignore Python cache artifacts |
| `62a9d19e94ad` | 2026-06-02 | 다름 | commit 보존 또는 비교 제외 | feat(nmr): checkpoint native GIAO h10 scaffolding |
| `00383a72474e` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | bench: use pip-installed OpenQP for NMR matrix |
| `95458b4b0e05` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | bench: run OpenQP NMR CGO baseline with memory |
| `d8597b038fce` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | bench: add GIAO NMR validation scaffold |
| `de0d4fceac6f` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | feat: gate NMR GIAO interface |
| `9690deccf1f1` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | docs: feat/dft-nmr branch handoff (RHF/DFT CGO NMR base; moved to private) |
| `d5de87f3f48e` | 2026-05-30 | 다름 | commit 보존 또는 비교 제외 | docs: Phase-0 validation summary + magnetic-symmetry derivation note |
| `61fbda1c6cec` | 2026-05-30 | 다름 | commit 보존 또는 비교 제외 | Phase 0: coupled HF/hybrid ground-state NMR magnetic response (validation) |
| `d1458b100c4f` | 2026-05-30 | 다름 | commit 보존 또는 비교 제외 | docs: add MRSF-TDDFT NMR shielding design draft |
| `8e24156d5bd4` | 2026-05-30 | 다름 | commit 보존 또는 비교 제외 | Fix NMR PSO integral antisymmetry for off-center nuclei |
| `24ff43963f8a` | 2026-05-30 | 다름 | commit 보존 또는 비교 제외 | chore: ignore build/ and .venv/ directories |
| `77f03a1a22c9` | 2026-05-30 | 다름 | commit 보존 또는 비교 제외 | Add native RHF/DFT NMR shielding prototype |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.gitignore` | 다름 |
| A | `HANDOFF_DFT_NMR.md` | 다름 |
| A | `examples/NMR/H2O_RHF-NMR.inp` | 동일 |
| M | `include/oqp.h` | 다름 |
| M | `pyoqp/oqp/library/runfunc.py` | 다름 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| A | `scripts/nmr_giao_benchmark_matrix.py` | 다름 |
| M | `source/integrals/int1.F90` | 다름 |
| M | `source/integrals/int_rys.F90` | 다름 |
| M | `source/integrals/mod_1e_primitives.F90` | 다름 |
| A | `source/modules/MRSF_NMR_DESIGN.md` | 다름 |
| A | `source/modules/NMR_BUILD_NOTES.md` | 다름 |
| A | `source/modules/NMR_MAGNETIC_SYMMETRY_NOTE.md` | 다름 |
| A | `source/modules/NMR_SHIELDING_STATUS.md` | 다름 |
| A | `source/modules/nmr_giao_debug.F90` | 다름 |
| A | `source/modules/nmr_shielding.F90` | 다름 |
| A | `tests/fixtures/nmr/benchmark_results/nmr_giao_benchmark_results.csv` | 다름 |
| A | `tests/fixtures/nmr/benchmark_results/nmr_giao_benchmark_results.json` | 다름 |
| A | `tests/fixtures/nmr/benchmark_results/nmr_giao_benchmark_results.md` | 다름 |
| A | `tests/fixtures/nmr/generate_pyscf_cgo_reference.py` | 다름 |
| A | `tests/fixtures/nmr/generate_pyscf_giao_reference.py` | 다름 |
| A | `tests/fixtures/nmr/giao_benchmark_matrix.json` | 다름 |
| A | `tests/fixtures/nmr/pyscf_cgo_reference.json` | 다름 |
| A | `tests/fixtures/nmr/pyscf_giao_reference.json` | 다름 |
| A | `tests/test_nmr_coupled.py` | 다름 |
| A | `tests/test_nmr_gauge_interface.py` | 다름 |
| A | `tests/test_nmr_giao_benchmark_path.py` | 다름 |
| A | `tests/test_nmr_giao_h10_live.py` | 다름 |
| A | `tests/test_nmr_giao_h10_twoe_live.py` | 다름 |
| A | `tests/test_nmr_giao_native_scaffold.py` | 다름 |
| A | `tests/test_nmr_shielding.py` | 다름 |
| A | `tests/test_rohf_status_and_interface.py` | 다름 |
