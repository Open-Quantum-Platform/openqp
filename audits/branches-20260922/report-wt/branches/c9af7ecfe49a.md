# draft/grd2-petite-optin

SHA: `c9af7ecfe49a270c9d737d0c1fcc3331fdd3645a`  
판정: **동일 끝점 PR 병합 확인**  
분야: Symmetry (검색용 분류)  
마지막 commit: 2026-08-12T08:57:59+09:00 / Cheol Ho Choi  
제목: Merge main after PR #320 into grd2 opt-in branch

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / draft/grd2-petite-optin](https://github.com/karmachoi/openqp/tree/c9af7ecfe49a270c9d737d0c1fcc3331fdd3645a)
- [gitlab / draft/grd2-petite-optin](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/c9af7ecfe49a270c9d737d0c1fcc3331fdd3645a)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `5653b78b1b76906989815cd0e1918f8e658f0628`
- 앞선 커밋 30 / 뒤처진 커밋 478
- non-merge patch: main과 일치 1, 다름 25
- merge/empty 등 patch 비교 제외: 4
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 14; 그중 현재 main과 동일 3, 다름 11

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #323: Symmetry reduction wiring: grd2 needs a per-caller opt-in, and the XC reduction is dead code](https://github.com/Open-Quantum-Platform/openqp/pull/323) — state=closed; merged=2026-08-12T00:25:57Z; merge commit in main=True

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `c9af7ecfe49a` | 2026-08-12 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge main after PR #320 into grd2 opt-in branch |
| `78959c0657ce` | 2026-08-11 | 다름 | commit 보존 또는 비교 제외 | fix(examples): give the XC-reduction deck a standard-frame geometry |
| `db21503528ab` | 2026-08-11 | 다름 | commit 보존 또는 비교 제외 | fix(tests): repair the CI break from the merge, and put a number behind the XC reduction |
| `f45b606b9184` | 2026-08-10 | 다름 | commit 보존 또는 비교 제외 | fix(mrsf): keep IXCORE runs out of the symmetry coverage substitution |
| `0c149542fdd1` | 2026-08-10 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge origin/main into draft/grd2-petite-optin |
| `587b55c0ecc1` | 2026-08-06 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'fork/fix/mrsf-davidson-symmetry-block-coverage' into _p323 |
| `f8d2d65ea294` | 2026-08-06 | 다름 | commit 보존 또는 비교 제외 | docs(symmetry): correct the cost of declining the full tier for DFT |
| `c0567ed24a98` | 2026-08-06 | 다름 | commit 보존 또는 비교 제외 | docs(symmetry): the labelling guard tests T^T S T, and say why that is observable |
| `0055a886cb27` | 2026-08-06 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): decline the full tier when a functional is active |
| `bd346ffe6768` | 2026-08-06 | 다름 | commit 보존 또는 비교 제외 | fix(dft): report a no-op XC reduction as inactive, and make the projector test real |
| `b30b861b318d` | 2026-08-06 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): validate the MO labels with the operator they are built from |
| `5d7022895844` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | refactor(symmetry): take the shell purity flag from the library, not a dimension test |
| `b2590affee37` | 2026-08-05 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge remote-tracking branch 'fork/draft/integral-symmetry' into fix/mrsf-davidson-symmetry-block-coverage |
| `87efacb35890` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(dft): project the reduced XC skeleton with the abelian operations |
| `96fe9c69a0da` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | test(symmetry): skip the petite engagement test without a native runtime |
| `43b9aec0eca5` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | feat(tdhf): warn when the Davidson start vectors leave a symmetry block unseeded |
| `1e00d1e609b7` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(dft): the XC symmetry reduction was dead code; revive it on the live path |
| `55e259e5c296` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): make the petite-list reduction engage on spherical bases |
| `1f831b681419` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(grd2): require a per-caller opt-in for the petite reduction |
| `8eae4edfb64e` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): validate the MO-labelling AO maps against the overlap matrix |
| `e96b42082a98` | 2026-08-05 | 일치 | commit 보존 또는 비교 제외 | docs: [draft] Give grd2 the same per-caller petite opt-in that int2 has |
| `db078b59f4e9` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | docs: record the petite-list integral reduction finding (no code yet) |
| `cfe50c9cc2f2` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(symmetry): stop tagging s and p shells as spherical when labelling MOs |
| `a00227abef28` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | mrsf: warn loudly when a symmetry block cannot be given a trial vector |
| `4a0fd9da9b92` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | Revert "mrsf: refuse to leave a symmetry block unseeded, and say so" |
| `af77131f4228` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | mrsf: refuse to leave a symmetry block unseeded, and say so |
| `5f8c40370a2c` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | mrsf: never displace the last trial vector of an irrep |
| `908496c007e1` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(mrsf): only substitute trial vectors when the guess has slack |
| `9f517067b2ef` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(mrsf): cover every symmetry block in the Davidson guess; symmetry on by default |
| `c4ab6dd2be5e` | 2026-08-05 | 다름 | commit 보존 또는 비교 제외 | fix(mrsf): seed the Davidson guess widely enough to reach every symmetry block |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.gitignore` | 다름 |
| A | `docs/planned/grd2-optin.md` | 다름 |
| A | `examples/other/h2o_rhf_6-31g_bhhlyp_intsym.inp` | 동일 |
| A | `examples/other/h2o_rhf_6-31g_bhhlyp_intsym.json` | 다름 |
| M | `source/dftlib/dft_gridint_energy.F90` | 동일 |
| M | `source/integrals/grd2.F90` | 다름 |
| M | `source/modules/hf_gradient.F90` | 다름 |
| M | `source/modules/tdhf_gradient.F90` | 다름 |
| M | `source/modules/tdhf_mrsf_gradient.F90` | 다름 |
| M | `source/modules/tdhf_sf_gradient.F90` | 다름 |
| M | `source/scf_addons.F90` | 다름 |
| M | `source/tdhf_mrsf_lib.F90` | 다름 |
| A | `tests/test_grd2_petite_optin.py` | 다름 |
| A | `tests/test_xc_symmetry_reduction_live.py` | 동일 |
