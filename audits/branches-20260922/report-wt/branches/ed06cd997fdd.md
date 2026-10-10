# perf/xc-numerical-kernel

SHA: `ed06cd997fdd3be9f383611d4e49652ecae4742e`  
판정: **동일 끝점 PR 병합 확인**  
분야: DFT / grid (검색용 분류)  
마지막 commit: 2026-06-07T07:51:17+09:00 / Cheol Ho Choi  
제목: Merge branch 'main' into perf/xc-numerical-kernel

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / perf/xc-numerical-kernel](https://github.com/karmachoi/openqp/tree/ed06cd997fdd3be9f383611d4e49652ecae4742e)
- [gitlab / perf/xc-numerical-kernel](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/ed06cd997fdd3be9f383611d4e49652ecae4742e)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `0f396ec8250f937a83465f2de24bf6c0cae486c3`
- 앞선 커밋 9 / 뒤처진 커밋 734
- non-merge patch: main과 일치 0, 다름 7
- merge/empty 등 patch 비교 제외: 2
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 10; 그중 현재 main과 동일 2, 다름 8

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #182: Speed up grid XC integration ~2x at full-core OpenMP](https://github.com/Open-Quantum-Platform/openqp/pull/182) — state=closed; merged=2026-06-06T23:18:12Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `ed06cd997fdd` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into perf/xc-numerical-kernel |
| `81b9512ef09f` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | dftlib: clamp pruning-layer index at the last radial bucket |
| `57cb55e1e6ba` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into perf/xc-numerical-kernel |
| `14f9e70a1494` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | Fix Docker (glibc 2.31) build: define _GNU_SOURCE for RTLD_DEFAULT |
| `bfe1662bee55` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | dftlib: split outer grid layers into finer angular patches |
| `a14e127489b7` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | dftlib: pair-bound AO threshold and live-only significance scan |
| `07aa95aed319` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | dftlib: prescreen shells per grid slice before AO evaluation |
| `72edb962c8d8` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | dftlib: evaluate per-point density reductions with BLAS ddot/dgemv |
| `6d9e8b626700` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | dftlib: speed up grid XC integration ~1.7x at full-core OpenMP |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `source/CMakeLists.txt` | 다름 |
| M | `source/basis_tools.F90` | 다름 |
| M | `source/dftlib/dft.F90` | 다름 |
| M | `source/dftlib/dft_gridint.F90` | 다름 |
| M | `source/dftlib/dft_gridint_energy.F90` | 다름 |
| M | `source/dftlib/dft_gridint_fxc.F90` | 다름 |
| M | `source/dftlib/dft_molgrid.F90` | 다름 |
| M | `source/mathlib/CMakeLists.txt` | 동일 |
| A | `source/mathlib/blas_thread.F90` | 동일 |
| A | `source/mathlib/blas_thread_ctl.c` | 다름 |
