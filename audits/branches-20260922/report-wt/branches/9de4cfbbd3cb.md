# feat/windows-intel-build

SHA: `9de4cfbbd3cb17f75acb36d77651cba2ea5575bc`  
판정: **동일 끝점 PR 병합 확인**  
분야: Build / CI / release (검색용 분류)  
마지막 commit: 2026-08-21T07:51:56+09:00 / Cheol Ho Choi  
제목: README: drop the standalone-archive section

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / feat/windows-intel-build](https://github.com/karmachoi/openqp/tree/9de4cfbbd3cb17f75acb36d77651cba2ea5575bc)
- [gitlab / feat/windows-intel-build](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/9de4cfbbd3cb17f75acb36d77651cba2ea5575bc)
- `local-11 / feat/windows-intel-build`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `cbc4eff856cec812058c4f68866ef95fa634d2de`
- 앞선 커밋 71 / 뒤처진 커밋 453
- non-merge patch: main과 일치 0, 다름 70
- merge/empty 등 patch 비교 제외: 1
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 19; 그중 현재 main과 동일 12, 다름 7

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #359: Windows support: build and run OpenQP with Intel oneAPI (ifx + icx + MKL ILP64)](https://github.com/Open-Quantum-Platform/openqp/pull/359) — state=closed; merged=2026-08-20T23:17:50Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `9de4cfbbd3cb` | 2026-08-21 | 다름 | commit 보존 또는 비교 제외 | README: drop the standalone-archive section |
| `9bc1f3fe429d` | 2026-08-21 | 다름 | commit 보존 또는 비교 제외 | README: point platform detail at the installation guide |
| `85099c6711b1` | 2026-08-21 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge main after #360 |
| `b2735803492c` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: mandate MKL ILP64 on Windows, reject ddX there, test the toolchain |
| `071545c0171e` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Hand the standalone bundles to OQP Studio |
| `84b7d22b55b7` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Pass the module array itself, not a pointer to it |
| `5927c1ab91aa` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | macOS bundle: advertise the floor it actually has, and verify it |
| `eb9fb4a9c9cc` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Declare the Win32 signatures for the module enumeration |
| `7a2078ff3f59` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Fix the guard: keep the metadata helpers at module level |
| `1bb64838c336` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: run the DFT-D4 checks on Windows too |
| `a6a86ae7fcfe` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: pin the MKL runtime, and load our own DFT-D4 DLLs first |
| `531d1cd4a557` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: detect the auto-selected BLAS, and reject static Windows outright |
| `e7f179831719` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: RPATH for ddX, D4-free cdefs, and reject multi-config Windows |
| `9061d051c224` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: tag checks for bundles, Windows RTLD, and BLAS mangling |
| `ed181ac2624c` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Promote Windows to a required wheel platform |
| `54658581f6ed` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: hold the DLL handles, fix the policy test, accept Windows wheels |
| `382d5af60459` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: resolve paths before judging them, and skip only the real ID |
| `6feecf950d5f` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Review: keep ELF strict, and verify the macOS precondition |
| `bfb8caa12c56` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Address the Codex review: four packaging defects |
| `bab21239cd2b` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Accept delocate 0.13's RPATH-free macOS wheel layout |
| `8f747db891f2` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Move the delocate 0.13 macOS fix out to its own PR |
| `48732b1ee074` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Tolerate only the two optional D4 entry points, not every missing symbol |
| `d2b5168a1240` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Drop the temporary branch triggers from the bundle workflows |
| `329ed6ffcb70` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Correct the stale BLAS rationale in the Windows bundle header |
| `bace6dce216b` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Do not read a dylib's own install ID as a dangling @rpath edge |
| `c9c13a3ebfd8` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Freeze the Linux bundle with a shared-libpython interpreter |
| `6ca9d68a0557` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Accept delocate 0.13's RPATH-free macOS wheel layout |
| `361bef1a2988` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | bundles: run the two new workflows on this branch until they land |
| `98ef51a71882` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Add a standalone Linux bundle |
| `d9194e32d4e2` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Add a standalone macOS bundle |
| `5a6a61b46bea` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows bundle: stamp the commit into the archive |
| `97bd66353302` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: fix the wheel smoke test, settle the triggers, document the download |
| `c34f883abe95` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: build the Windows wheel against MKL and let pip supply it |
| `22cd6895c444` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: drop DYNAMIC_ARCH from the Windows OpenBLAS and check its symbols properly |
| `6bc4f9ab2fdd` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: report which dgemm spelling the Windows OpenBLAS exports |
| `ff77772781f6` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: link a static OpenBLAS on Windows |
| `953e0ae21c1e` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: save the Windows OpenBLAS cache as soon as it is built |
| `f72a0966bc38` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: match Fortran symbol naming to the BLAS being linked |
| `7176c85f1ca8` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: cap the Windows OpenBLAS dispatch list at AVX2 |
| `19083d3a7c64` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: let OpenBLAS's x86-64 kernels compile with icx on Windows |
| `de8b156a775a` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: build the Windows OpenBLAS with its own thread pool |
| `c31c41a03e5e` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: build an unrenamed ILP64 OpenBLAS for the Windows leg |
| `37a8eb273061` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: hand the DFT-D4 subprojects the OpenBLAS import library on Windows |
| `5d8a56478e63` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: keep the Windows leg experimental until DFT-D4 links there |
| `3c5799b637d8` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: give the Windows build its externals root and a pkg-config |
| `8f7410e7f9cd` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: pass the Intel compilers through a toolchain file on Windows |
| `e235d459592e` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: point CMake at the absolute Intel compiler paths on Windows |
| `574976fe37eb` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: select ifx/icx through CIBW_ENVIRONMENT on Windows |
| `eab1ab21ecb8` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | wheels: build a Windows wheel too |
| `e5154f8da6ed` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | tests: update the DFT-D4 packaging contract for the Windows guards |
| `7f4241a81257` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | windows-bundle: copy package metadata into the freeze |
| `a9da68fce59f` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | windows-bundle: trim the MKL payload and make the DLL assertions version-agnostic |
| `8e082b512b80` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | windows-bundle: own entry script, and discover the runtime DLLs instead of assuming names |
| `e88bb4b6da79` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | windows-bundle: take the last stdout line when locating the oqp package |
| `2012b7296e7e` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: standalone zip distribution (unzip and run openqp.exe) |
| `fb94dde2b84b` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: use the shared Intel Fortran runtime across liboqp and the DFT-D4 DLLs |
| `faf4ec94344e` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: take the DFT-D4 DLLs from the build tree |
| `aa8133fb7fce` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: let CMake export the DFT-D4 stack's symbols |
| `85a793a10721` | 2026-08-20 | 다름 | commit 보존 또는 비교 제외 | Windows: build the DFT-D4 stack as DLLs instead of dropping it |
| `21b87cc61d95` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | windows-intel CI: run the example suite (informational) |
| `90c0e6ac44a1` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Windows: keep the pyoqp data installs, skip only the cffi API-mode extension |
| `c0492b74c4f6` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Tolerate optional backends missing from liboqp at import time |
| `5f4a2ab32bf5` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Windows: map the routec dynamic-loader seam onto kernel32; give the bundled externals a build type |
| `097910229cfc` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | ENABLE_DFTD4: include the patch SHA256 / cache-key computation in the guard |
| `9eb994f28100` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | ENABLE_DFTD4: also guard the corresponding-source installs and the 'patch' prerequisite |
| `635b67e28d86` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Make the native DFT-D4 backend optional (ENABLE_DFTD4); default OFF on Windows |
| `bc3f9efaa792` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Windows: keep build paths under the 250-character object-path cap |
| `39d387ac3a43` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | windows-intel CI: save the oneAPI install to cache right after installing |
| `848fc7d5f509` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Windows: link the DFT-D4 stack statically into liboqp.dll (OQP_DFTD4_SHARED) |
| `cb9193238c6c` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | windows-intel CI: install ifx+icx+MKL from the unified oneAPI installer (setup-fortran 2026.1 ships no icx) |
| `ad2ed6d2950a` | 2026-08-19 | 다름 | commit 보존 또는 비교 제외 | Windows/Intel experimental build: liboqp.dll naming, symbol export, DLL search path, CI probe |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| A | `.github/actions/setup-oneapi-windows/action.yml` | 동일 |
| M | `.github/scripts/wheel_smoke_test.py` | 동일 |
| M | `.github/workflows/build_wheels.yml` | 동일 |
| A | `.github/workflows/windows-intel.yml` | 동일 |
| M | `CMakeLists.txt` | 다름 |
| M | `README.md` | 다름 |
| M | `cmake/oqp_functions.cmake` | 동일 |
| M | `cmake/sanitize_macos_package_rpaths.cmake.in` | 동일 |
| M | `external/CMakeLists.txt` | 다름 |
| M | `pyoqp/CMakeLists.txt` | 동일 |
| M | `pyoqp/oqp/__init__.py` | 다름 |
| M | `pyoqp/oqp/runtime.py` | 동일 |
| M | `pyoqp/oqp_cffi_build.py` | 동일 |
| M | `pyproject.toml` | 다름 |
| M | `source/CMakeLists.txt` | 다름 |
| M | `source/modules/routec_bridge.F90` | 동일 |
| M | `source/modules/routec_sig.F90` | 동일 |
| M | `tests/test_dftd4_shared_interface.py` | 다름 |
| M | `tests/test_runtime_root_resolution.py` | 동일 |
