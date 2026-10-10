# feat/solvent-backend-spike

SHA: `85c4caed3ea450bc10660259cb94801538e25a4a`  
판정: **동일 끝점 PR 병합 확인**  
분야: Solvation / PCM (검색용 분류)  
마지막 commit: 2026-06-07T15:46:51+09:00 / Cheol Ho Choi  
제목: Merge origin/main into feat/solvent-backend-spike

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / feat/solvent-backend-spike](https://github.com/karmachoi/openqp/tree/85c4caed3ea450bc10660259cb94801538e25a4a)
- [gitlab / feat/solvent-backend-spike](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/85c4caed3ea450bc10660259cb94801538e25a4a)
- `local-15 / pr176`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `628f2b0ca5ddfd5ccc09bf532e7c73d518855ad0`
- 앞선 커밋 140 / 뒤처진 커밋 727
- non-merge patch: main과 일치 0, 다름 137
- merge/empty 등 patch 비교 제외: 3
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 37; 그중 현재 main과 동일 9, 다름 28

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #176: feat(pcm): energy-only ddPCM (ddX) solvent — runtime path, energy fix, independent validation](https://github.com/Open-Quantum-Platform/openqp/pull/176) — state=closed; merged=2026-06-07T07:12:17Z; merge commit in main=True

## 같은 이름의 PR — 끝점이 다르므로 전체 반영 증거가 아님

- [private #2: feat(pcm): energy-only ddPCM (ddX) solvent — runtime path, energy fix, independent validation](https://github.com/karmachoi/openqp-private/pull/2) — open; PR head `49b717fe59a1`; merged=아니오

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `85c4caed3ea4` | 2026-06-07 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge origin/main into feat/solvent-backend-spike |
| `8b8b1cc07060` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | build(ddx): build ddX as ILP64 from source via CMake (DDX_ROOT optional) |
| `204c27a5ceee` | 2026-06-07 | 다름 | commit 보존 또는 비교 제외 | docs->tests: drop PR-added PCM design docs, pin invariants in a source test |
| `ae6232939b5c` | 2026-06-04 | 다름 | commit 보존 또는 비교 제외 | fix: gate PCM SCF scope on enabled flag |
| `5122745f9720` | 2026-06-04 | 다름 | commit 보존 또는 비교 제외 | test: drop removed thymine NEB example expectation |
| `084648af7db9` | 2026-06-04 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge upstream main into ddPCM solvent PR branch |
| `e5c02f9763a0` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | refactor(pcm): de-name PySCF in PCM provenance strings + dependent tests |
| `624ac3962b16` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): genericize PySCF wording in PR docs (keep benchmarks/provenance) |
| `63d3bef269a8` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): state l<=2 (quadrupole) scope; record higher-l as to-be-done |
| `c1056e6b6ad0` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): report e_pcm as -0.5<phi_exact,q_cav> — independently validated |
| `39e46fb37d0a` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | validation(pcm): root-cause refinement — multipole source is a dead end |
| `80ae954772b2` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | validation(pcm): independent ddPCM cross-check (PySCF) — H2O not validated |
| `b3ac2150bf6a` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | test(pcm): defer DFT (BHHLYP/PBE) benchmark rows out of the verified scope |
| `300424de28d8` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | test(pcm): trim benchmark matrix to the defensible subset |
| `3b72823bcf8d` | 2026-06-03 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): review note for the QM source-term fix |
| `7e93ca83af3e` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | Populate closed-shell PCM benchmark matrix |
| `4e14577e810e` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | Stabilize HF DFT PCM diagnostics |
| `662d14dc6eb3` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | Derive ddX PCM Fock sign scale |
| `434c04641257` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | Add QM multipole PCM source validation |
| `0f6659396b6f` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | Fix ddX PCM QM source consistency |
| `99bc32d4839f` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): handoff spec for the QM source-term fix |
| `c28e92bd190e` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | test(pcm): ddX-enabled build + run results for the QM PCM benchmark gate |
| `b5f9180ad6f5` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | test(pcm): first Fortran-driven QM PCM benchmark gate (H2O) + convention diagnostics |
| `9c2d604b66d3` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | test(pcm): populate Tier-1 ddX trusted-reference regression gate |
| `ad9b64d75833` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | test(pcm): add ddX literature/reference validation gate |
| `973699b2d062` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | chore: ignore Python bytecode caches (__pycache__/, *.pyc) |
| `47226fec18a9` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | refactor(pcm): collapse duplicate runtime PCM path |
| `903f0eb4cd8d` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): add session handoff note |
| `38f2186aad2e` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record dual-path reconciliation task |
| `cdddfbaffb31` | 2026-05-31 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | merge: integrate ddX-gated PCM energy path onto private solvent line |
| `9d3026ecdaf1` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record provisional conventions and ddX validation gate |
| `9f0777d9216d` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): wire ddX-gated PCM energy path |
| `1f5b173eece0` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): preserve ddx payload roundtrip boundary |
| `377dc62bfe74` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs: preserve PCM call-site shape audit boundary |
| `42ce6e857edf` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): audit ddx call-site shape metadata |
| `109316a9df8f` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): audit final call-site bridge shape |
| `9fb3c5c35dee` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record final calc-fock bridge in ddx seam |
| `d31a98849a74` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): reject malformed molecule handoff payloads |
| `e333a474ccb1` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs: record calc_jk_xc PCM prototype boundary |
| `b1996ed0fed0` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): preserve molecule payload mapping guard |
| `f5e45fafd51b` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): guard molecule runtime payload type |
| `a86a92146fb0` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record old-buffer provenance in validation matrix |
| `074ca175fec8` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record native calc_fock shape guard |
| `df84aec8a2ef` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): require literal ddx provisional opt-in |
| `2bf7a6b5498a` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record old-buffer bridge diagnostics |
| `65fad4575759` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record ddx payload consumer guard |
| `b918a0755668` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record non-mapping payload guard |
| `f3217fdee797` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): reject non-mapping runtime payloads |
| `bbfe8c532f98` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record molecule payload roundtrip contract |
| `b85961d70a39` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record required runtime payload metadata |
| `54aea24accf7` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): reject boolean runtime payload numerics |
| `fe851ad403a5` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): require nbf in runtime payload consumer |
| `93ae414dbbd1` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): preserve runtime payload shape metadata |
| `a7381aee5e9c` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record backend validation payload gate |
| `f8a59f5402aa` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): require backend validation provenance in runtime payload |
| `376cd0a879b9` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): require runtime payload shape metadata |
| `580470bac41f` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record runtime payload shape contract |
| `a3b7501dc422` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): validate runtime payload shape metadata |
| `3f8cfa290d5c` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): preserve runtime payload shape metadata |
| `0157809e95d8` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): preserve reviewed payload shape metadata |
| `57106a210ff9` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): expose calc-fock shape metadata |
| `c3ad77846ea6` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): document calc-fock call-site shape guard |
| `ab1672347863` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): add guarded calc_fock call-site bridge |
| `6eb85f74a851` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): expose packed AO length in calc_fock handoff |
| `231899c7138e` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): guard calc_fock reaction potential shape |
| `e02eebb33f23` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record incremental fock validation guard |
| `98a5a3c7ddc8` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label disabled molecule handoff validation status |
| `3398128a0f86` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): describe guarded calc_fock handoff |
| `4161f4e2b2e3` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): report native old-buffer trigger fields |
| `e2d824d85579` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): report reference fock old-buffer triggers |
| `523d9ca9007d` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): preserve fock update scope metadata |
| `74f584b32527` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): report scf old-buffer provenance |
| `32c4f5d8a13d` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): detail calc_fock old-buffer guard diagnostics |
| `ba712343bcc1` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): clarify incremental fock old-buffer guard |
| `3161044064c0` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): factor SCF incremental audit metadata |
| `88ef24eda1a2` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): expose SCF old-buffer request metadata |
| `50718a8075d1` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): require integer runtime payload nbf |
| `05f215a310cb` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record SCF-state calc_fock request guard |
| `b7594ce34135` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): guard calc fock request from scf state |
| `96cc54204997` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | test(pcm): guard calc_fock request mode |
| `772a708f8735` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): guard old Fock state in calc_fock handoff |
| `541147740bea` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | docs: record non-incremental PCM calc_fock request guard |
| `25e44ea11256` | 2026-05-29 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): gate reference calc_fock requests |
| `9e79319439b8` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): block unvalidated incremental fock handoff |
| `454358beb57b` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label empty calc-fock handoff scope |
| `0e8b97958fab` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): reject density leakage in runtime payload |
| `9ff963ca96ce` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | docs(pcm): record guarded calc_fock handoff chain |
| `8772d9937c5c` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): gate calc_fock handoff from molecule payload |
| `6f66d33bfc40` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): validate runtime payload nbf |
| `d4b06a905188` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): carry nbf through calc_fock handoff |
| `022863024a24` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): package calc_fock reaction handoff |
| `95e63446fac5` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): require reference-density runtime payloads |
| `e5998efb0314` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): gate restored reaction potential payload |
| `ad2265347123` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): round-trip runtime payload metadata |
| `a000bcb8f99c` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): record opt-in reference reaction energy |
| `10a474533fb3` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): add reference scf runtime payload |
| `94cb678b844c` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): add SCF PCM energy bookkeeping slot |
| `d306ecd8a2f2` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): package reference energy handoff |
| `2443aad95921` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label provisional ddX first-scope handoff |
| `db07f49a9302` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label reference energy bookkeeping scope |
| `133627edb77a` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label provisional ddx charge mapping scope |
| `3585042cf248` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label reference fock update scope |
| `dfdf3ab4ef9d` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): label disabled runtime coupling contract |
| `9785dc82e4b4` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): guard reference fock block count |
| `f62ed4970348` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): thread reference reaction field through calc_fock |
| `41347c6261e4` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): add opt-in SCF reaction handoff |
| `4e255177d943` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): label reference coupling scope |
| `5ffaf2264d01` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): add reference reaction field fock helper |
| `0c636f502d6c` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): guard reference reaction fock updates |
| `71dab6f59e2e` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): guard reference SCF energy terms |
| `b67a6cfc6d40` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(solvent): package reference PCM coupling inputs |
| `7db39a444a2b` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(solvent): reject nonfinite PCM seam values |
| `9c6a6f6bf497` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): expose packed density for phi cav handoff |
| `24bdf95515df` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): guard reference phi_cav handoff |
| `1d3685753265` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): guard reference reaction-field contract |
| `f34b1db10f7b` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | test(pcm): reject empty reference density blocks |
| `b70650c50b12` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): guard reference density block count |
| `4432a1e59cb8` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): guard reference scf total density handoff |
| `8e8551a07006` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): split provisional ddx cavity charges |
| `ce7ed7ea0de2` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | fix(pcm): reject empty ddX cavity handoff |
| `d27204e2cd58` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): validate provisional ddX reaction-field inputs |
| `a10775b783c3` | 2026-05-28 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): guard provisional ddX reaction charges |
| `a9aff7af06e7` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | docs: add PCM validation matrix |
| `be9f8883f894` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | fix: preserve malformed PCM epsilon for diagnostics |
| `0601b23ba04a` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | fix: guard pcm dielectric validation |
| `15b8f1fe920c` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | feat: tighten PCM first-scope guardrails |
| `f66f87280299` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | feat(pcm): validate backend model pairings |
| `65e224b92958` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | feat: record ddX cavity charge derivative seam |
| `d7c608a717e7` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: expose ddX q cavity handoff |
| `cc30d1582d20` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | docs: label ddX xi as projected q |
| `f2400997ed0c` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: add explicit ddX PCM adapter smoke |
| `2f928d89c37e` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: expose unweighted electrostatic potential |
| `3d5ffb2b3251` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | docs: map ddX SCF coupling quantities |
| `923deb973caa` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: expose external-charge potential seam |
| `b06394888eb4` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: add OpenQP ddX adapter API |
| `52892c79e644` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | test: add ddX adapter lifecycle smoke |
| `331b191a462d` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | build: add optional ddX link smoke test |
| `448b0b7ff416` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | docs: record ddX source-build probe |
| `849080d62481` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | docs: add ddX backend API probe |
| `761b4d5d0bf2` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | feat: scaffold PCM solvent input |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.gitignore` | 다름 |
| M | `CMakeLists.txt` | 다름 |
| A | `cmake/FindDDX.cmake` | 다름 |
| A | `cmake/patches/ddx-v0.8.0-ilp64.patch` | 동일 |
| M | `external/CMakeLists.txt` | 다름 |
| M | `include/oqp.h` | 다름 |
| A | `pyoqp/oqp/library/solvent.py` | 동일 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| M | `pyoqp/oqp/utils/input_parser.py` | 다름 |
| A | `scripts/pcm_benchmark_matrix.py` | 다름 |
| A | `scripts/pcm_independent_ddpcm_validation.py` | 다름 |
| A | `scripts/pcm_pyscf_reference.py` | 다름 |
| A | `scripts/pcm_trusted_reference_diagnostics.py` | 다름 |
| A | `scripts/pcm_vacuum_pyscf_validation.py` | 다름 |
| A | `scripts/pcm_validate_all.py` | 다름 |
| M | `source/CMakeLists.txt` | 다름 |
| M | `source/integrals/int1.F90` | 다름 |
| M | `source/scf_addons.F90` | 다름 |
| A | `source/solvent_ddx_adapter.c` | 다름 |
| A | `source/solvent_ddx_adapter.h` | 다름 |
| A | `source/solvent_pcm.F90` | 다름 |
| M | `source/types.F90` | 다름 |
| A | `spikes/001-ddx-api-probe/README.md` | 다름 |
| A | `spikes/001-ddx-api-probe/probe_pyddx_point_charges.py` | 다름 |
| A | `tests/data/pcm_literature_benchmarks.json` | 다름 |
| A | `tests/data/pcm_trusted_reference_diagnostics.json` | 동일 |
| A | `tests/data/pcm_vacuum_pyscf_validation.json` | 동일 |
| A | `tests/ddx_adapter_smoke.c` | 동일 |
| A | `tests/ddx_link_smoke.c` | 동일 |
| A | `tests/test_ddx_cmake_scaffold.py` | 다름 |
| A | `tests/test_ddx_scf_integration_seam.py` | 동일 |
| M | `tests/test_geometric_optimizer.py` | 다름 |
| A | `tests/test_pcm_canonical_runtime_path.py` | 다름 |
| A | `tests/test_pcm_energy_path.py` | 동일 |
| A | `tests/test_pcm_literature_benchmarks.py` | 다름 |
| A | `tests/test_pcm_scaffold.py` | 동일 |
