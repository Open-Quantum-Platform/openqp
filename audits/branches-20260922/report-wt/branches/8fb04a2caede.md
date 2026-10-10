# feat/scf-cleanup

SHA: `8fb04a2caede38ffdda6911c9239ef74393d9211`  
판정: **고유 변경·반영 여부 검토**  
분야: SCF / response (검색용 분류)  
마지막 commit: 2026-06-08T17:27:47+09:00 / Cheol Ho Choi  
제목: ci: retry the build to survive transient 504 dependency downloads

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-private / feat/scf-cleanup](https://github.com/karmachoi/openqp-private/tree/8fb04a2caede38ffdda6911c9239ef74393d9211)
- [gitlab / github-private/feat/scf-cleanup](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/8fb04a2caede38ffdda6911c9239ef74393d9211)
- `local-15 / pr-188-check`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `da646bf76164f5161c3e9b17d5913ded4daca48c`
- 앞선 커밋 31 / 뒤처진 커밋 721
- non-merge patch: main과 일치 0, 다름 26
- merge/empty 등 patch 비교 제외: 5
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 23; 그중 현재 main과 동일 4, 다름 19

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.


## 같은 이름의 PR — 끝점이 다르므로 전체 반영 증거가 아님

- [upstream #188: SCF: converger cleanup, robust escalation ladder, and a documented TRAH solver](https://github.com/Open-Quantum-Platform/openqp/pull/188) — closed; PR head `d1ec2b8fdc67`; merged=2026-06-08T10:08:24Z

## 이 끝점을 포함하는 같은 분야의 후속 브랜치 후보

- [github-personal/feat/scf-cleanup, gitlab/feat/scf-cleanup, local-15/feat/scf-cleanup](3ba40812f00f.md) — 추가 커밋 3

이 목록은 ancestry 관계이며 기능 대체나 최선의 통합 대상을 보증하지 않는다.

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `8fb04a2caede` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | ci: retry the build to survive transient 504 dependency downloads |
| `a388704bcc5e` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | scf: address review of the escalation ladder + fix OTR arch flags on ARM |
| `35782682808a` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | feat(scf): selector sees SCF reference (is_rohf/is_uhf); refresh shipped model |
| `7b81cfd45483` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | feat(scf): ship distilled ML converger selector (converger_type=ml goes live) |
| `11b09d6f0e41` | 2026-06-08 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'origin/main' into feat/scf-cleanup |
| `f0ad77c2286b` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | docs: README — SOC and ddPCM are merged capabilities, not upcoming |
| `77fe78a66c7c` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | WIP: SCF manager (converger_type=auto\|ml) -- feature-based converger selection |
| `1665cf9034f3` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | build: fix OpenTrustRegion AVX-512 SIGILL on AVX2 nodes (portable arch flags) |
| `cd7357992558` | 2026-06-08 | 다름 | commit 보존 또는 비교 제외 | SCF: incremental-Fock refresh to remove the DIIS noise-floor (match PySCF) |
| `9f71506d965e` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | SCF: in-loop stagnation detection -> early escalation handoff (orchestration) |
| `d181b3bf861f` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'karmachoi/feat/scf-cleanup' into feat/scf-cleanup |
| `d71eabef21df` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | build: ENABLE_OPENTRAH CMake selector + README compile-section trim |
| `f2fea1ffbae7` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into feat/scf-cleanup |
| `094a8f3b32a9` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'origin/main' into feat/scf-cleanup |
| `4fce8d2a7ee3` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | TRAH: drop the extra Fock build per macro-iteration in Steihaug-CG (#15) |
| `72c941f65055` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | TRAH: make the final canonicalization MOM-aware (#14) |
| `748c61b98dc5` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | TRAH: route micro-solver subspace algebra through BLAS (dgemm/dgemv) |
| `74561c0d85be` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | SOSCF: keep the no-Hessian-reset default after dropping soscf_reset_mod |
| `54cc19cd0940` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Clean up SCF converger options + de-'native' TRAH naming + cite Helmich-Paris |
| `2d9aac9b4d32` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Hybrid converger: SOSCF escalation step before TRAH in the robust SCF ladder (#20) |
| `e8eab992dc07` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | build: propagate CMAKE_PREFIX_PATH to the OpenTrustRegion external subbuild |
| `6a9c7e054836` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'origin/main' into feat/scf-cleanup |
| `75d6d37a496a` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Fix garbage nuclear repulsion: compute nenergy in Fortran every SCF; zero ecp_zn for non-ECP |
| `91f8b8926633` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Make native TRAH the default for ground-state SCF energy (trh_impl=auto) |
| `fb7c6915135e` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Native TRAH: tighten orbitals past FP-energy convergence + fresh final Fock |
| `9841863ffa23` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Native TRAH: stability check (lowest-eig of H) escapes unstable solutions |
| `abdd6b3bc759` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Native TRAH: stability-aware convergence + Fock canonicalization (native-gated) |
| `81b7572b8af0` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Native TRAH: MOM-aware, aug-Hessian+RTV default; keep OTR as the SCF default |
| `cbe6385764dd` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Native TRAH: random trial vectors in aug-Hessian Davidson -> tier-3 parity |
| `2274344b7e3d` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | WIP: experimental native-TRAH micro-solvers (opt-in) |
| `7accbb3b7163` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | SCF: conditional options logging + native Fortran TRAH solver |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.github/workflows/CI.yml` | 다름 |
| M | `CMakeLists.txt` | 다름 |
| M | `README.md` | 다름 |
| M | `examples/SCF/h2o_rhf_6-31g_pbe_vdiis.inp` | 동일 |
| M | `examples/SCF/h2o_rohf_6-31g_pbe_vshift.inp` | 동일 |
| M | `examples/SCF/h2o_uhf-s_6-31g_pbe_vdiis.inp` | 동일 |
| M | `examples/SCF/h2o_uhf-t_6-31g_pbe_vdiis.inp` | 동일 |
| M | `external/CMakeLists.txt` | 다름 |
| M | `include/oqp.h` | 다름 |
| A | `pyoqp/oqp/library/scf_selector_model.py` | 다름 |
| M | `pyoqp/oqp/library/single_point.py` | 다름 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/utils/file_utils.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `pyoqp/setup.py` | 다름 |
| M | `pyproject.toml` | 다름 |
| M | `source/CMakeLists.txt` | 다름 |
| M | `source/basis_api.F90` | 다름 |
| M | `source/scf.F90` | 다름 |
| M | `source/scf_converger.F90` | 다름 |
| A | `source/trah_converger.F90` | 다름 |
| M | `source/types.F90` | 다름 |
| M | `tests/test_single_point_scf_fallback.py` | 다름 |
