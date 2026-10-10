# OpenQP 개인·private 개발 브랜치 전수 조사

기준일: 2026-09-22. 조사 대상은 GitHub `karmachoi/openqp`,
`karmachoi/openqp-private`, 현재 GitLab `open-quantum-platform/internal/openqp`,
그리고 Ultra에서 확인한 해당 저장소의 로컬 개발 브랜치다.
GitLab의 과거 `cheol/openqp`는 현재 내부 그룹 저장소이며,
별도 `cheol/openqp-private` 프로젝트는 이번 GitLab 목록에서 발견되지 않았다.

**두 개인 저장소는 단순 중복이 아니다. 개발 이력 대부분은 GitLab에 보존되어
있지만, 아직 GitLab에 없는 최신 REKS 변경과 로컬 변경이 있다. 브랜치 이름이나
마지막 커밋 제목만 보고 삭제·archive하면 안 된다.**


## 전체 집계

스냅샷 UTC: `2026-09-21T22:55:39.472889+00:00`  
비교 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`

| 조사 대상 | 브랜치/ref 수 |
| --- | ---: |
| GitHub karmachoi/openqp | 244 |
| GitHub karmachoi/openqp-private | 152 |
| GitLab internal/openqp | 390 (개발·일반 379, PR 보존 11) |
| Ultra 로컬 16개 common directory | 326 |
| 합계 / 고유 SHA | 1112 / 533 |

| 고유 끝점 판정 | 개수 |
| --- | ---: |
| 고유 변경·반영 여부 검토 | 314 |
| 동일 끝점 PR 병합 확인 | 152 |
| main 이력에 포함 | 58 |
| non-merge patch 일치·merge 검토 | 5 |
| 별도 이력 | 4 |

## 상세 자료

- [검색·필터 가능한 전체 조사표](index.html)
- [원본별 전체 브랜치 CSV](all-branches.csv)
- [GitLab 미보존 로컬 브랜치](local-only.md)
- [215개 등록 worktree 상태](worktrees.md)
- [원시 스냅샷·patch·PR 증거](../evidence/)

## 조사 방법과 판정 범위

모든 유효한 브랜치를 원본별로 열거하고 같은 SHA를 묶었다. 각 끝점마다 현재
GitLab `main`과의 공통 조상, 앞선/뒤처진 커밋 수, 전체 추가 커밋 목록,
공통 조상 이후 변경된 파일, 해당 파일의 현재 main과의 일치 여부를 조사했다.
5,037개 코드 이력 커밋을 대상으로 non-merge patch-id를 계산했으며,
GitHub upstream·개인·private PR과 GitLab MR 기록을 대조했다.

판정은 서로 다르다.

- **main 이력에 포함:** 브랜치 끝점 자체가 현재 main의 조상이다.
- **동일 끝점의 PR 병합 확인:** PR head SHA가 브랜치 끝점과 일치하고,
  PR의 merge commit이 현재 main 이력에 있다. squash merge도 이 기준으로 찾았다.
- **non-merge patch 일치:** nonempty patch는 main에 있으나 merge commit의
  충돌 해결까지 증명한 것은 아니다. 정리 전에 최종 파일 비교가 필요하다.
- **검토 필요:** 위 조건으로 전체 반영을 증명하지 못했다. 브랜치 전체가
  미반영이라는 뜻도, 모든 변경이 새로운 기능이라는 뜻도 아니다.
- **별도 이력:** `gh-pages`, README 전용 등 main과 공통 조상이 없는 이력이다.

패치가 달라도 같은 기능을 다른 방식으로 구현했거나 여러 커밋을 합쳐 반영했을 수
있다. 따라서 `GitLab에 없는 patch` 수는 코드 비교의 우선순위를 정하는 증거이며,
그 수만큼 독립 기능이 빠졌다는 뜻이 아니다. branch date는 마지막 커밋 날짜이며
최근에 사람이 작업한 시각이나 branch 생성일을 뜻하지 않는다.

전수 조사는 Git ref·commit·patch·file·PR 수준에서 수행했다. 중요한 GPU,
NMR, PCM, SSC 등은 실제 소스와 handoff의 구현 경계도 확인했다. 모든 소스 행에
대한 수학적 증명, 전체 조합의 빌드, 분자 계산과 성능 재측정은 이번 조사 범위에
포함하지 않았다. M5 Air·M4 mini·Windows·HPC의 미push 로컬 브랜치는 조사하지
않았으므로, 그 장치까지 포함한 전체 데이터 보존을 보증하지 않는다.

## 우선 보존할 변경

### REKS: GitHub 최신 변경이 GitLab보다 앞서 있다

GitHub 두 저장소의 `feature/reks22-scf`는 모두 `53d1d7e8f1`을 가리킨다.
GitLab 일반 이름의 끝점 `3dfacba0`보다 16개, 보존용
`github-private/feature/reks22-scf`의 `06318be0`보다 14개 커밋 앞선다.
두 비교 모두 GitLab 쪽에만 있는 커밋은 0개다.

누락된 변경에는 trust-region gradient scale, spin-adapted configuration의
reference averaging, reference block과 core-to-virtual singles의 coupling,
REKS microstate ensemble의 오비탈 최적화가 포함된다. 해당 commit들은
2026-09-20~21에 작성됐다. 보존 우선순위가 가장 높지만, 이 조사에서는 push나
과학적 타당성 판단을 하지 않았다.

### 로컬에만 남은 브랜치

GitLab의 어느 브랜치 이력에도 끝점이 없는 로컬 ref는 63개이며 SHA로 묶으면
62개다. 57개 끝점에는 GitLab 전체 코드 이력과 patch-id가 일치하지 않는 변경이
있다. 4개는 non-merge 변경이 다른 SHA로 GitLab에 존재하며, `resolve-241`은
별도 merge commit 검토가 필요하다. 자세한 경로·변경 목록은 local-only 표에 있다.

특히 다음 작업은 보존·대조를 먼저 해야 한다.

- `feat/gpu-workspace-manager`: workspace 관리·검증 관련 3개 patch.
- `feat/gpu-metc-regression-test`: CPU/GPU regression 및 kernel-only timing 2개 patch.
- `feat/routec-dfjk-bridge`, `mac-rot-seam`: Route-C 연결부와 macOS ERI 연결의 로컬 변형.
- MRSF Hessian·response·screening 진단, S1 frequency, ddX lifetime 수정.
- QMRSF H4/Thiel/manuscript 조사, EKT 진단, CASSCF·CASPT2 개발 변형.

GPU workspace 관련 파일 자체는 다른 GPU 브랜치에도 있다. 위 patch 개수는
해당 로컬 변경과 동일한 patch를 찾지 못했다는 뜻이며, 기능 전체가 GitLab에
없다는 뜻은 아니다.

## 분야별 세부 판단

### GPU / rotation ERI

`gpu-mrsf-seam`의 끝점 `b4505281`은 upstream PR #265로 병합됐고 그 merge
commit은 현재 main에 있다. 현재 GPU 라이브러리 연결을 다시 통째로 가져올
필요는 없다. 반면 `feat/routec-dfjk-bridge`는 SCF J/K, GauXC Vxc,
2e-gradient, MRSF low-rank J/K를 연결하는 별도 개발선이다.

METC는 기본 구현, profile breakdown, kernel efficiency, accumulation layout,
warp atomics, two-phase accumulation으로 이어지는 이력이 있다. 단순 복사본으로
취급하지 말고 공통 조상과 옵션별 kernel을 비교해야 한다.
`feat/gpu-metc-two-phase-accum`에는 workspace/cache/Fortran bridge 코드도 있어,
앞서 찾은 로컬 GPU 관리 코드를 옮기기 전에 이 브랜치와 대조해야 한다.

`feat/gpu-xc-response`의 CUDA 파일은 packed slot마다
`response[idx] += density[idx] * kernel[idx]`를 실행하며, 소스 주석도 full XC
quadrature/cache 연결 이전의 scaffold라고 명시한다. 완성된 XC response로
표시하면 안 된다. private 쪽의 같은 이름은 다른 SHA의 더 긴 이력이다.

`feat/gpu-porting-unified`는 `GPU_PORTING_UNIFIED_HANDOFF.md`만 추가한다.
브랜치 이름과 달리 그 끝점 자체는 GPU 코드 통합 결과가 아니다.
`rot-seam-v12`는 rotation ERI dispatch와 ERI store/RAM budget 작업이고,
`perf/int`에는 rotation backend를 제거한 커밋도 있으므로 성능 브랜치를
무조건 최신 통합본으로 선택하면 안 된다.

권장 목적지는 GPU application code의 경우 그룹 `openqp-gpu`, ERI kernel은
독립 `libintRot`, OpenQP 입력·Fortran 연결은 `internal/openqp`다.

### NMR

CGO NMR, GIAO ground-state NMR, singlet TDA NMR, MRSF amplitude/CSF response가
서로 다른 개발선이다. `feat/mrsf-nmr-gate5a`는 amplitude response의 deferred
조건을 포함한다. `codex/mrsf-nmr-csf-response-20260907`에는 interacting
response와 integral-response covariance 진단이 있다.

`feature/mrsf-nmr-shielding`의 handoff는 SOMO rotation 관련 수정 이후에도
shielding 크기와 induced current에 대한 미해결 사항을 기록한다. 또한 그
브랜치의 당시 테스트 환경에서 `OPENQP_ROOT`가 없으면 실질 검사를 하지 않고
통과할 수 있다고 경고한다. 이는 해당 handoff의 역사적 기록이며 현재 모든
NMR 테스트에 같은 문제가 있다는 판정은 아니다. ground-state NMR의 PR 병합과
MRSF NMR 완성 여부를 구분해야 한다.

### HF/DFT/TDDFT/MRSF Hessian

`feat/analytic-hessian`, `feat/hf-dft-analytic-hessian`,
`feat/response-analytic-hessian-private`, `feat/mrsf-analytical-hessian-private`,
그리고 9월의 relaxed-density·XC-performance·screening·residual 진단이 공존한다.
일부 기본 Hessian은 upstream PR로 반영됐지만, private 브랜치는 항별 구현과
진단 자료를 포함하므로 단일 'Hessian 완료' 상태로 묶으면 안 된다.
가장 최근 날짜보다 현재 main과 비교한 실제 추가 항, regression, 계산 기록을
기준으로 선택해야 한다.

### QMRSF / doublet-quartet / REKS

`review/qmrsf-dk-original`, `feat/qmrsf-dk-covariant-seam`,
`feat/dk-fullspace-singles`, `feat/qmrsf-dk-dft`, `feat/qmrsf-dual-pathways`,
`alireza/*`, `feat/doublet-quartet-mrsf`는 서로 다른 이론·구현·검증 가정을 갖는
개발선이다. `review/qmrsf-dk-original_2`와 `feat/qmrsf-dk-covariant-seam`은
현재 같은 SHA다. `icPT2.STATUS.md`에는 active-only correlation과 큰 basis에서의
scalability 한계가 기록돼 있다. 이번에는 이론의 옳고 그름을 재판정하지 않았다.

연구 code와 benchmark 자료를 보존하고, 구현 비교에 들어갈 때 MRSF theory skill의
검증 절차를 적용해야 한다. '버전이 뒤이므로 대체 가능'이라는 가정은 부적절하다.

### SOC / SSC / X2C

`feat/x2c-scalar`, `fix/soc2e-reference-density`, `feat/spin-spin-coupling`,
`ssc_tensor`, `feat/mrsf-soc-analytic-gradient`를 별도로 보존해야 한다.
private PR #3/#4/#5는 아직 open이며, PR의 존재는 병합 증거가 아니다.
`ssc_tensor` handoff는 aromatic-triplet ZFS의 reference-path 문제를 기록하고,
이전 method-limit 해석을 철회한다. SSC 결과를 검증 완료로 일반화하면 안 된다.
SOC gradient의 과거 benchmark commit 역시 이번에 재실행한 결과가 아니다.

### Solvation / PCM

ground-state ddX 연결, runtime/메모리 수명 수정, MRSF-PCM 진단을 구분해야 한다.
`feat/mrsf-pcm-spike`의 상태 문서에는 `energy_coupled=false`, state-specific /
transition / relaxed density가 deferred라고 명시돼 있으며, 비진단 요청을
거부하는 경계가 있다. 이 브랜치를 완성된 MRSF-PCM 계산 기능으로 소개하면 안 된다.
최근 로컬 `codex/ddx-lifetime-deferred-20260920` 변경도 별도 보존 후보다.

### EKT / PBC / spectroscopy

`feat/ea-fock-rebuild-k`와 `claude/ekt-mrsf-pbc-*`에는 k-point, PBC J/K,
EA Fock rebuild, Dyson pole-strength 작업이 함께 들어 있다. 이름에 EKT만
있다고 작은 property patch로 판단하면 안 된다. `feature/mrsf-pecd`는 외부
`openqp-spec` 소비와 export/API 코드가 있으며, GitLab MR !1은 spectroscopy와
Gelius 모델 통합 작업으로 아직 open이다. 별도 그룹 companion repository의
현재 코드와 비교한 뒤 OpenQP 연결부만 선택적으로 가져와야 한다.

### Wavefunction / DFTB / optimization / SCF / build

FCI/CASSCF/SA-CASSCF/CASPT2/NEVPT2, MP2, CCSD(T), quantum Hamiltonian export가
공존하며 `openqp-fci-option`은 수백 커밋의 장기 개발선이다. 작은 기능 이름과
실제 변경 범위가 다르다. DFTB/xTB, QM/MM, geomeTRIC/native optimizer,
SCF/TRAH/Davidson, symmetry/grid, CI/build는 많은 PR이 이미 병합돼 있다.
branch별 PR head 일치 여부를 상세표에 기록했으므로, 반영된 기능을 다시 합치기
전에 그 기록을 확인하면 된다.

특히 private PR #6 `Allow macOS LP64 BLAS builds`는 현재의 OpenQP ILP64-only
개발 규칙과 충돌한다. 역사적 자료로 보존하되 현재 main으로 그대로 합치면 안 된다.

## 중복과 로컬 작업공간

GitHub 두 저장소에서 같은 branch name이 서로 다른 SHA인 경우가 17개다.
GPU XC/persistent buffers, NMR, Hessian, PCM, NAMD 등이 포함된다.
GitLab의 `github-private/*`, `github-karmachoi/*`는 이런 충돌을 보존한 것이므로,
일괄 중복 삭제 대상이 아니다.

Ultra의 16개 Git common directory에 등록된 worktree 215개 중 202개가 존재했고,
50개에서 tracked 변경 또는 untracked 항목을 확인했다. 이 개수에는 같은 브랜치의
중복 작업공간, 연구 기록, 빌드 산출물 등이 함께 포함된다. 실제 파일 diff와
소유 작업을 확인하기 전에는 정리하면 안 된다.

`gradient-fix/openqp/.git/refs/heads/feat/gpu-porting-unified 2`라는 유효하지 않은
중복 ref 파일도 있다. 그 안의 commit `67df0c10`은 여러 GitLab METC 브랜치에서
도달 가능하므로 해당 commit 자체는 유실 상태가 아니다. 파일은 수정하지 않았다.

## 권장 순서

1. GitHub REKS 최신 14개 커밋과 로컬 미보존 변경을 먼저 개별 브랜치로 보존한다.
2. GPU METC·workspace·XC scaffold·Route-C 연결의 실제 차이를 비교하여
   GPU/ERI/OpenQP host 위치별로 선택적으로 통합한다.
3. NMR, SSC, MRSF Hessian, QMRSF, PCM은 구현 경계와 검증 조건을 가진 별도
   개발 과제로 유지한다. 미해결 연구 가정은 코드 정리로 해소되지 않는다.
4. main 포함 또는 정확한 PR 병합이 확인된 브랜치는 로컬 변경과 추가 커밋이
   없는지 확인한 다음 보관·정리 후보로 삼는다.

이번 조사에서는 원본 저장소의 branch·commit·remote·worktree를 변경하거나,
push·merge·archive·삭제하지 않았다. 아래 상세 자료는 조사 시점의 스냅샷이다.


## 분야별 브랜치 찾아보기

### GPU (21)

- [claude/elegant-hypatia-wZpoB · a1c4f06d](branches/a1c4f06de592.md) — 고유 변경·반영 여부 검토
- [feat/gpu-metc · 5e4b23eb](branches/5e4b23eba181.md) — 고유 변경·반영 여부 검토
- [feat/gpu-metc-accumulation-layout · a3b8a4a6](branches/a3b8a4a67914.md) — 고유 변경·반영 여부 검토
- [feat/gpu-metc-atomics-v2 · c03935a1](branches/c03935a17b6e.md) — 고유 변경·반영 여부 검토
- [feat/gpu-metc-kernel-efficiency · 70220d95](branches/70220d95deeb.md) — 고유 변경·반영 여부 검토
- [feat/gpu-metc-regression-test · fc92929b](branches/fc92929bbdbb.md) — 고유 변경·반영 여부 검토
- [feat/gpu-metc-two-phase-accum · beb51a26](branches/beb51a2617ad.md) — 고유 변경·반영 여부 검토
- [feat/gpu-porting-unified · ffabe60c](branches/ffabe60cb2f4.md) — 고유 변경·반영 여부 검토
- [feat/gpu-profile-breakdown · 1cacad57](branches/1cacad57621f.md) — 고유 변경·반영 여부 검토
- [feat/gpu-workspace-manager · 6413215e](branches/6413215eeb41.md) — 고유 변경·반영 여부 검토
- [feat/gpu-xc-response · c03b17e8](branches/c03b17e83566.md) — 고유 변경·반영 여부 검토
- [feat/gpu-xc-response · 76366540](branches/763665403a4b.md) — 고유 변경·반영 여부 검토
- [feat/gpu-xc-response-rebase · 7f9eda08](branches/7f9eda08b208.md) — 고유 변경·반영 여부 검토
- [feat/routec-dfjk-bridge · 4d8afd6e](branches/4d8afd6e122f.md) — 고유 변경·반영 여부 검토
- [feat/routec-dfjk-bridge · d9a752d4](branches/d9a752d4002a.md) — 고유 변경·반영 여부 검토
- [fix/oqp-timer-parser-blockers · 901c3e2b](branches/901c3e2bc723.md) — 고유 변경·반영 여부 검토
- [gpu-mrsf-seam · b4505281](branches/b45052816768.md) — 동일 끝점 PR 병합 확인
- [mac-rot-seam · 3d4984f8](branches/3d4984f8c2d5.md) — 고유 변경·반영 여부 검토
- [perf/gpu-metc-persistent-buffers · fb8309e1](branches/fb8309e16725.md) — 고유 변경·반영 여부 검토
- [perf/gpu-metc-persistent-buffers · dd4f6351](branches/dd4f63515d10.md) — 고유 변경·반영 여부 검토
- [rot-seam-v12 · 9636a98b](branches/9636a98b85dc.md) — 고유 변경·반영 여부 검토
### ERI / integrals (17)

- [bench/rys-cart · 91626817](branches/91626817b923.md) — 고유 변경·반영 여부 검토
- [chore/spherical-digest-gate · 437d8ba0](branches/437d8ba0981e.md) — main 이력에 포함
- [codex/ispher-rot-rys · 24f83cbf](branches/24f83cbf8149.md) — 동일 끝점 PR 병합 확인
- [codex/ispher-rot-rys · 8a3b2610](branches/8a3b261020ec.md) — 고유 변경·반영 여부 검토
- [codex/pr199-ispher-fixes · fd7e8fb4](branches/fd7e8fb412d7.md) — 고유 변경·반영 여부 검토
- [draft/integral-symmetry · 96fe9c69](branches/96fe9c69a0da.md) — 고유 변경·반영 여부 검토
- [feat/ispher-basis-option · 6f51e89c](branches/6f51e89cc189.md) — 고유 변경·반영 여부 검토
- [feat/rotaxis-direct-pure · 20f504d4](branches/20f504d410c2.md) — 고유 변경·반영 여부 검토
- [fix/int2-blas-pin-serialises-openmp · fe68ed23](branches/fe68ed23f126.md) — 동일 끝점 PR 병합 확인
- [fix/int2-openmp-runtime-defaults · a9519e17](branches/a9519e1738e7.md) — 고유 변경·반영 여부 검토
- [fix/int2-openmp-workshare-barrier · 035ba274](branches/035ba274b4ef.md) — 동일 끝점 PR 병합 확인
- [fix/int2-openmp-workshare-sync · e6a46d19](branches/e6a46d19b9bc.md) — main 이력에 포함
- [fix/macos-int2-flat-workshare · d2f6a27e](branches/d2f6a27e67dd.md) — 동일 끝점 PR 병합 확인
- [perf/int · 0d8cda74](branches/0d8cda74da43.md) — 고유 변경·반영 여부 검토
- [perf/int-upstream · f479b7bb](branches/f479b7bbb9b4.md) — 동일 끝점 PR 병합 확인
- [perf/int2-pair-decode-opt · 392e95a5](branches/392e95a54c20.md) — non-merge patch 일치·merge 검토
- [perf/int2-pair-decode-opt · 7e179be1](branches/7e179be1cafe.md) — 동일 끝점 PR 병합 확인
### NMR (13)

- [claude/cool-ritchie-vaMB6 · 6c7f0fa8](branches/6c7f0fa8e750.md) — 동일 끝점 PR 병합 확인
- [claude/cool-ritchie-vaMB6 · 1f4c0da6](branches/1f4c0da61efe.md) — 고유 변경·반영 여부 검토
- [claude/tda-nmr-cgo-20260905 · 194a0504](branches/194a0504befc.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-nmr-csf-response-20260907 · fff81acc](branches/fff81acc29a9.md) — 고유 변경·반영 여부 검토
- [codex/nmr-integral-response-validation-20260908 · edd73016](branches/edd7301654eb.md) — 고유 변경·반영 여부 검토
- [codex/nmr-metric-validation-20260908 · a7c6d07b](branches/a7c6d07b66e9.md) — 고유 변경·반영 여부 검토
- [codex/tda-nmr-origin-regression-20260906 · 98930fd7](branches/98930fd73cf1.md) — 고유 변경·반영 여부 검토
- [feat/dft-nmr · 9690decc](branches/9690deccf1f1.md) — 고유 변경·반영 여부 검토
- [feat/giao-nmr · fde1a89a](branches/fde1a89a3dd1.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-nmr-gate5a · 71e064cc](branches/71e064cc4e01.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-nmr-gate5a · 4faaa117](branches/4faaa117420e.md) — 고유 변경·반영 여부 검토
- [feature/mrsf-nmr-shielding · 2ba421ee](branches/2ba421ee2734.md) — 고유 변경·반영 여부 검토
- [fix/nmr-debug-print-gating · 468a6840](branches/468a68404f35.md) — 동일 끝점 PR 병합 확인
### REKS (3)

- [feature/reks22-scf · 53d1d7e8](branches/53d1d7e8f1ab.md) — 고유 변경·반영 여부 검토
- [feature/reks22-scf · 3dfacba0](branches/3dfacba07f18.md) — 고유 변경·반영 여부 검토
- [github-private/feature/reks22-scf · 06318be0](branches/06318be09e97.md) — 고유 변경·반영 여부 검토
### QMRSF / higher spin (31)

- [agent/qmrsf-fth-dc · 24a90fd9](branches/24a90fd95601.md) — 고유 변경·반영 여부 검토
- [alireza/qmrsf-dk-clean · 7d8ea0e4](branches/7d8ea0e4a7b1.md) — 고유 변경·반영 여부 검토
- [alireza/qmrsf-dk-paper-format · 87b4779f](branches/87b4779f2a91.md) — 고유 변경·반영 여부 검토
- [alireza/xqmrsf-on-dk · 25f4c147](branches/25f4c147fa38.md) — 고유 변경·반영 여부 검토
- [claude/all-open-issues-20260821 · 848f8cb8](branches/848f8cb8ced9.md) — 동일 끝점 PR 병합 확인
- [claude/mrsf-adc-code-p2w5ah · 5e66a0a9](branches/5e66a0a956fb.md) — 고유 변경·반영 여부 검토
- [claude/ptc-tenno-umu7ko · 9041f1c7](branches/9041f1c7f14d.md) — 고유 변경·반영 여부 검토
- [codex/fix-s3r-seam-continuity-20260814-019ff5e3 · 0bacd2d6](branches/0bacd2d64794.md) — 고유 변경·반영 여부 검토
- [codex/h4-original2-mac-20260812 · 79b47081](branches/79b470814a3f.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-state-analysis · 827f5b81](branches/827f5b810b90.md) — 동일 끝점 PR 병합 확인
- [codex/qmrsf-be-c-operator-audit-20260901 · f7572f76](branches/f7572f7667dd.md) — 고유 변경·반영 여부 검토
- [codex/qmrsf-dk-original2-ilp64 · 63870142](branches/638701427c50.md) — main 이력에 포함
- [codex/qmrsf-dk-original2-validation · 3e5c61b3](branches/3e5c61b3be96.md) — 고유 변경·반영 여부 검토
- [codex/qmrsf-dk-singlet-vectors · 8025f916](branches/8025f9167c31.md) — 고유 변경·반영 여부 검토
- [codex/qmrsf-h4-scan-20260918-h4agent · 31768813](branches/31768813cde4.md) — 고유 변경·반영 여부 검토
- [codex/qmrsf-paper-audit-20260918-01a0aa9f · 3c7d475e](branches/3c7d475e8536.md) — 고유 변경·반영 여부 검토
- [codex/qmrsf-sr-canonical-roks-20260814-019fff32 · daf0f444](branches/daf0f4448772.md) — 고유 변경·반영 여부 검토
- [codex/qmrsf-thiel-20260916-01a0aa9f · ebf92ef2](branches/ebf92ef2b380.md) — 고유 변경·반영 여부 검토
- [codex/review-qmrsf-dual-pathways · 2c26faab](branches/2c26faabec16.md) — 고유 변경·반영 여부 검토
- [feat/dk-fullspace-singles · 98190b03](branches/98190b03a483.md) — 고유 변경·반영 여부 검토
- [feat/dk-fullspace-singles · bc75f319](branches/bc75f319d8af.md) — 고유 변경·반영 여부 검토
- [feat/doublet-quartet-mrsf · 53b56a20](branches/53b56a207e9e.md) — 고유 변경·반영 여부 검토
- [feat/doublet-quartet-mrsf · 3ee72f82](branches/3ee72f82fc52.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-ensemble-reference · 2804a0ec](branches/2804a0ec3e5d.md) — 고유 변경·반영 여부 검토
- [feat/qmrsf-dk-covariant-seam · c4dd34f8](branches/c4dd34f895a2.md) — 고유 변경·반영 여부 검토
- [feat/qmrsf-dk-dft · befa4023](branches/befa4023bf15.md) — 고유 변경·반영 여부 검토
- [feat/qmrsf-dk-live · 70c9825a](branches/70c9825af713.md) — 고유 변경·반영 여부 검토
- [feat/qmrsf-dual-pathways · 83789078](branches/8378907823cd.md) — 고유 변경·반영 여부 검토
- [fix/dk-drop-delta · e699f80b](branches/e699f80bee2a.md) — 고유 변경·반영 여부 검토
- [github-private/pull/7/head · e2db16da](branches/e2db16da55ec.md) — 고유 변경·반영 여부 검토
- [review/qmrsf-dk-original · 6c78c59e](branches/6c78c59ea57b.md) — 고유 변경·반영 여부 검토
### Hessian / frequency (38)

- [claude/mrsf-hessian-seam-residual-20260904 · 5f20bba4](branches/5f20bba43fb7.md) — 고유 변경·반영 여부 검토
- [claude/mrsf-hessian-xc-perf-20260904 · fdc60421](branches/fdc60421be08.md) — 고유 변경·반영 여부 검토
- [claude/practical-hypatia-4tcBa · 6d967ba4](branches/6d967ba41b54.md) — 고유 변경·반영 여부 검토
- [claude/practical-hypatia-4tcBa · fd529e12](branches/fd529e126182.md) — 고유 변경·반영 여부 검토
- [claude/rks-hessian-fock-deriv-speed-20260904 · a9a3755f](branches/a9a3755fa2c2.md) — 고유 변경·반영 여부 검토
- [claude/rks-hessian-fock-deriv-speed-20260904 · 3896e891](branches/3896e891be24.md) — 고유 변경·반영 여부 검토
- [claude/rys-hessian-skeleton · 89113742](branches/89113742c84d.md) — 고유 변경·반영 여부 검토
- [claude/tddft-hessian-diagnostics-20260904 · f603db28](branches/f603db28bf99.md) — 고유 변경·반영 여부 검토
- [claude/tddft-hessian-relaxed-density-fix-20260904 · 9795db63](branches/9795db633e0e.md) — 고유 변경·반영 여부 검토
- [codex/claude-mrsf-hessian-handoff-20260902 · 74ff3cd7](branches/74ff3cd78f73.md) — 고유 변경·반영 여부 검토
- [codex/import-tddft-hessian-20260831 · 11fe5698](branches/11fe56982622.md) — 동일 끝점 PR 병합 확인
- [codex/import-tddft-hessian-20260831 · 576be1a7](branches/576be1a7a3d0.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-curvature-origin-20260920-01a07f50 · dbd2214d](branches/dbd2214df861.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-curvature-postfix-20260920-01a07f50 · a0b9d20a](branches/a0b9d20ab724.md) — main 이력에 포함
- [codex/mrsf-formaldehyde-debug-20260902 · 46c9b929](branches/46c9b92950f9.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-s0-diagnostics-20260901 · 718bcc33](branches/718bcc3329ca.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-s0-screening-fix-20260901 · 5fcabbd3](branches/5fcabbd3fea6.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-s1-freq-20260907-01a07a3a · 30d8ed9c](branches/30d8ed9c7cc1.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-tddft-hessian-20260901 · 879b6b60](branches/879b6b60953b.md) — 고유 변경·반영 여부 검토
- [codex/private-pr391-closed-shell-hessian-opt-20260902 · c081e7b8](branches/c081e7b8404f.md) — 고유 변경·반영 여부 검토
- [draft/hess-unique-displacements · 55d1a3c3](branches/55d1a3c3da43.md) — 고유 변경·반영 여부 검토
- [draft/hessian-child-reorient · 4ea6283b](branches/4ea6283be570.md) — 동일 끝점 PR 병합 확인
- [feat/analytic-hessian · 947781d7](branches/947781d78bea.md) — 고유 변경·반영 여부 검토
- [feat/analytic-hessian · 977de410](branches/977de410161b.md) — 고유 변경·반영 여부 검토
- [feat/analytic-hessian-mrsf-private · 77ae1d26](branches/77ae1d26a09a.md) — 고유 변경·반영 여부 검토
- [feat/analytic-hessian-public-clean · f9e8bf35](branches/f9e8bf35e857.md) — 고유 변경·반영 여부 검토
- [feat/hf-dft-analytic-hessian · a2e67e99](branches/a2e67e99fdce.md) — 고유 변경·반영 여부 검토
- [feat/hf-dft-analytic-hessian · 4476613f](branches/4476613fd1fb.md) — 고유 변경·반영 여부 검토
- [feat/hf-dft-hessian · 4704c2ed](branches/4704c2edf05f.md) — 동일 끝점 PR 병합 확인
- [feat/hf-dft-hessian · 3c840740](branches/3c8407403161.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-analytical-hessian-private · dc431b55](branches/dc431b559b92.md) — 고유 변경·반영 여부 검토
- [feat/response-analytic-hessian-private · 20ee42c8](branches/20ee42c856e9.md) — 고유 변경·반영 여부 검토
- [feat/uhf-rohf-hessian · 9a576f77](branches/9a576f77b8a8.md) — 동일 끝점 PR 병합 확인
- [feat/uhf-rohf-hessian · a204de53](branches/a204de539551.md) — 고유 변경·반영 여부 검토
- [fix-h2o-analytic-hessian-reference · 801c6562](branches/801c65624ab3.md) — 동일 끝점 PR 병합 확인
- [fix/hessian-followups · 6b483c68](branches/6b483c68b7c5.md) — 동일 끝점 PR 병합 확인
- [mpi · f8f07e23](branches/f8f07e23ad37.md) — main 이력에 포함
- [rys-hessian-skeleton · 733620aa](branches/733620aaaf5a.md) — 동일 끝점 PR 병합 확인
### NAC / NAMD (47)

- [agent/analytic-directional-fssh-20260830 · fa7b709c](branches/fa7b709c5776.md) — 고유 변경·반영 여부 검토
- [agent/nac-exact-perf-visible-20260830 · 492887b8](branches/492887b801be.md) — 고유 변경·반영 여부 검토
- [agent/nac-zpredict-namd-fgr-20260830 · 02d3dc09](branches/02d3dc09d88b.md) — 고유 변경·반영 여부 검토
- [agent/odp-umbrella-native · a04b04b9](branches/a04b04b98dc1.md) — main 이력에 포함
- [agent/odp-umbrella-native · 717ecfbb](branches/717ecfbb6312.md) — main 이력에 포함
- [agent/rigorous-namd-gates · b1412f15](branches/b1412f159aac.md) — main 이력에 포함
- [agent/soc-namd-gates · 8cbb9319](branches/8cbb931976ea.md) — main 이력에 포함
- [agent/uracil-namd-campaign-20260830 · 57fe2c53](branches/57fe2c539f2c.md) — 고유 변경·반영 여부 검토
- [agent/zhu-nakamura-hopping · 33c37859](branches/33c37859f4c7.md) — 고유 변경·반영 여부 검토
- [bugfix-num-NAC · 56e7a16c](branches/56e7a16cc604.md) — main 이력에 포함
- [bugfix-num-NAC · 3289ffd9](branches/3289ffd9a066.md) — 고유 변경·반영 여부 검토
- [bugfix-num-NAC · 04e58a47](branches/04e58a47ac4c.md) — 고유 변경·반영 여부 검토
- [claude/affectionate-ramanujan-ovj7F · 2a33f192](branches/2a33f19298cd.md) — 고유 변경·반영 여부 검토
- [claude/bethe-surface-scattering-20260917 · 1f047d56](branches/1f047d562c6b.md) — main 이력에 포함
- [claude/busy-babbage-mKs8g · cc1daaac](branches/cc1daaacaa33.md) — 고유 변경·반영 여부 검토
- [claude/namd-analytic-nac-handoff-20260907 · 406f5a0a](branches/406f5a0a6576.md) — 고유 변경·반영 여부 검토
- [claude/namd-mo-reuse-20260905 · a4db6b97](branches/a4db6b972430.md) — 고유 변경·반영 여부 검토
- [claude/namd-selected-nac-20260908 · 6523a0f1](branches/6523a0f1e11b.md) — 고유 변경·반영 여부 검토
- [codex/analnac-merged-regression-20260919 · 7b3ba519](branches/7b3ba5190f8b.md) — main 이력에 포함
- [codex/mrsf-response-diagnostic-20260914-01a07f50 · fcdf2091](branches/fcdf2091bc71.md) — 고유 변경·반영 여부 검토
- [codex/nac-analytic-meci · 85c025de](branches/85c025de13c7.md) — 고유 변경·반영 여부 검토
- [codex/nac-performance · 8eb9a8d9](branches/8eb9a8d9730a.md) — 고유 변경·반영 여부 검토
- [codex/nac-zres-refinement · e4e41f41](branches/e4e41f410e2c.md) — 고유 변경·반영 여부 검토
- [codex/nac-zres-tolerance · b6ae311c](branches/b6ae311c6d9c.md) — 고유 변경·반영 여부 검토
- [codex/namd-local-continuation-20260914-01a07f50 · e3cdf526](branches/e3cdf5268a77.md) — 고유 변경·반영 여부 검토
- [codex/namd-no-fort6-20260827 · 7001c65c](branches/7001c65c6d10.md) — 고유 변경·반영 여부 검토
- [codex/pr205-soc-namd-options · c4bc21fc](branches/c4bc21fc5026.md) — 고유 변경·반영 여부 검토
- [codex/soc-namd-options · 39bfc439](branches/39bfc439c8d4.md) — 고유 변경·반영 여부 검토
- [codex/soc-namd-options · 4e014247](branches/4e014247bffb.md) — 고유 변경·반영 여부 검토
- [codex/soc-namd-options · 7b1bb759](branches/7b1bb7597dea.md) — 고유 변경·반영 여부 검토
- [codex/static-analytic-nac-input-20260910-01a07f50 · 1f508e4f](branches/1f508e4f9f2e.md) — 고유 변경·반영 여부 검토
- [codex/student-handoff-20260831 · f78eb044](branches/f78eb044268a.md) — 고유 변경·반영 여부 검토
- [codex/thymine-five-20260919-01a07f50 · eb8877de](branches/eb8877debe50.md) — 고유 변경·반영 여부 검토
- [codex/thymine-namd-20260917-01a07f50 · 3bc8946e](branches/3bc8946efb86.md) — main 이력에 포함
- [codex/thymine-namd-20260917-01a07f50 · f7570dc9](branches/f7570dc92d4c.md) — 고유 변경·반영 여부 검토
- [dev · 7198ce0e](branches/7198ce0e896e.md) — 고유 변경·반영 여부 검토
- [docs/readme-soc-namd-qmmm · d1859852](branches/d185985275da.md) — 동일 끝점 PR 병합 확인
- [feat/mrsf-analytic-nac · ecda4b84](branches/ecda4b84db96.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-analytic-nac-private · 48401534](branches/484015344ce9.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-analytic-nac-private · 5d823d51](branches/5d823d51f6e6.md) — 고유 변경·반영 여부 검토
- [feat/namd-qmmm-compile-gating · 76feb04a](branches/76feb04a729b.md) — 고유 변경·반영 여부 검토
- [feat/rt-mrsf-state-propagation · 925c18e3](branches/925c18e351f7.md) — 고유 변경·반영 여부 검토
- [main · aafc6e13](branches/aafc6e136fb6.md) — 고유 변경·반영 여부 검토
- [nac · b71b864e](branches/b71b864ea404.md) — 고유 변경·반영 여부 검토
- [nac-lagrangian · 81d41e62](branches/81d41e6200fb.md) — 고유 변경·반영 여부 검토
- [namd-qmmm · 8750d7f7](branches/8750d7f7ed7c.md) — 고유 변경·반영 여부 검토
- [pr405 · 3a9157eb](branches/3a9157eb70ef.md) — main 이력에 포함
### SOC / SSC / X2C (22)

- [agent/fix-umrsf-mixed-exchange · e6c52024](branches/e6c520241940.md) — 동일 끝점 PR 병합 확인
- [agent/umrsf-dipole-moment-calculations-fix · fd9ea411](branches/fd9ea411c274.md) — 동일 끝점 PR 병합 확인
- [chore/trim-example-reference-jsons · ed86fa80](branches/ed86fa809f0d.md) — 동일 끝점 PR 병합 확인
- [codex/fix-ssc-reference-contraction · 0b070671](branches/0b070671a4a4.md) — 고유 변경·반영 여부 검토
- [codex/pi-mkl-pack-unpack-workaround · d832e766](branches/d832e7667422.md) — 고유 변경·반영 여부 검토
- [codex/umrsf-gradient-zvector-20260701 · 79e1e77b](branches/79e1e77bd48a.md) — 고유 변경·반영 여부 검토
- [codex/x2c-soc-science-20260826 · e51cf8a0](branches/e51cf8a016bf.md) — 고유 변경·반영 여부 검토
- [codex/x2c-ssc-pes-20260826 · 40d151a7](branches/40d151a7853b.md) — 고유 변경·반영 여부 검토
- [docs/socgrad-channel-audit · 723598af](branches/723598af2f18.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-soc-analytic-gradient · 4d107ec0](branches/4d107ec011df.md) — 고유 변경·반영 여부 검토
- [feat/spin-spin-coupling · 942990f0](branches/942990f01c64.md) — 고유 변경·반영 여부 검토
- [feat/spin-spin-coupling · 9ae677d6](branches/9ae677d66e56.md) — 고유 변경·반영 여부 검토
- [feat/x2c-scalar · b5371300](branches/b5371300dc83.md) — 고유 변경·반영 여부 검토
- [fix/soc2e-reference-density · c0c8a510](branches/c0c8a510f128.md) — 고유 변경·반영 여부 검토
- [fix/umrsf-energy-issues · a3e7878a](branches/a3e7878ad113.md) — 동일 끝점 PR 병합 확인
- [fix/umrsf-energy-issues · c264ffe8](branches/c264ffe8a049.md) — 고유 변경·반영 여부 검토
- [github-private/pull/3/merge · 3229c731](branches/3229c73151e9.md) — 고유 변경·반영 여부 검토
- [github-private/pull/4/merge · cbfba3a5](branches/cbfba3a5938d.md) — 고유 변경·반영 여부 검토
- [github-private/pull/5/merge · 2ef71c25](branches/2ef71c2501d9.md) — 고유 변경·반영 여부 검토
- [soc-mrsf-pr166 · f6003304](branches/f60033049b89.md) — 고유 변경·반영 여부 검토
- [ssc_tensor · 42747b20](branches/42747b200a85.md) — 고유 변경·반영 여부 검토
- [umrsf_energy · 9228f7fa](branches/9228f7fa6918.md) — 고유 변경·반영 여부 검토
### EKT / PBC / spectroscopy (17)

- [agent/spec-gelius-merge-20260921 · b305c460](branches/b305c4608b27.md) — 고유 변경·반영 여부 검토
- [claude/analytic-ht-transition-dipole-20260903 · 060af26a](branches/060af26aa1ba.md) — 고유 변경·반영 여부 검토
- [claude/ekt-mrsf-pbc-bands-fixes-20260825 · 2b8bdf50](branches/2b8bdf502140.md) — 고유 변경·반영 여부 검토
- [claude/ekt-mrsf-pbc-fixes-on-live-20260826 · 0aafbee6](branches/0aafbee6eb29.md) — 고유 변경·반영 여부 검토
- [claude/mrsf-2pa-sos-20260917 · 41b3c7d1](branches/41b3c7d1d4c5.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-circular-dichroism · 1172ce48](branches/1172ce480493.md) — 고유 변경·반영 여부 검토
- [codex/vibrational-circular-dichroism · 2ccd8fad](branches/2ccd8fad90a9.md) — 고유 변경·반영 여부 검토
- [docs/xas-readme-clean · 313acecc](branches/313acecc8c44.md) — 고유 변경·반영 여부 검토
- [docs/xas-readme-energy-clarification · e18e138d](branches/e18e138dbed1.md) — 고유 변경·반영 여부 검토
- [dyson · ee9d988c](branches/ee9d988c4773.md) — 고유 변경·반영 여부 검토
- [feat/ea-fock-rebuild-k · 24649e2f](branches/24649e2f1ec5.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-ekt-ip-ea · a33b0a7e](branches/a33b0a7eae57.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-ekt-ip-ea · a3536076](branches/a3536076c606.md) — 고유 변경·반영 여부 검토
- [feature/mrsf-pecd · e1817074](branches/e18170749ee6.md) — 고유 변경·반영 여부 검토
- [pr393 · 6bb42083](branches/6bb4208333c1.md) — 고유 변경·반영 여부 검토
- [research/ekt-mrsf-w-diagnostics · 0e7038a5](branches/0e7038a5a5c3.md) — 고유 변경·반영 여부 검토
- [research/ekt-mrsf-wtilde-derivation · 9e2b94f7](branches/9e2b94f76120.md) — 고유 변경·반영 여부 검토
### Wavefunction methods (33)

- [claude/caspt2-analytic-gradient-20260817 · 11965fab](branches/11965fab7fa9.md) — 동일 끝점 PR 병합 확인
- [claude/caspt2-analytic-gradient-20260817 · 1b687296](branches/1b6872969208.md) — 고유 변경·반영 여부 검토
- [claude/caspt2-ci-phase-20260822 · fd37c71f](branches/fd37c71fac88.md) — 고유 변경·반영 여부 검토
- [claude/casscf-analytic-gradient-20260815 · 20d23572](branches/20d235721e66.md) — 고유 변경·반영 여부 검토
- [claude/ccsd-t-parallel-lnoacw · a23ec633](branches/a23ec6336d7a.md) — 동일 끝점 PR 병합 확인
- [claude/ccsd-t-parallel-lnoacw · 65e7ca95](branches/65e7ca95ffc4.md) — 고유 변경·반영 여부 검토
- [claude/ci-irrep-selection-20260821 · 96e0cc07](branches/96e0cc0707d1.md) — 동일 끝점 PR 병합 확인
- [claude/ci-irrep-selection-20260821 · a532e2a6](branches/a532e2a6c4cf.md) — 고유 변경·반영 여부 검토
- [claude/ci-spin-cluster-20260821 · 83c48ad5](branches/83c48ad56b54.md) — 동일 끝점 PR 병합 확인
- [claude/ci-spin-cluster-20260821 · 5e3d27d5](branches/5e3d27d502bb.md) — 고유 변경·반영 여부 검토
- [claude/determined-galileo-523gm1 · 78271f90](branches/78271f9011e5.md) — 고유 변경·반영 여부 검토
- [claude/nevpt2-analytic-gradient-20260817 · 1116c056](branches/1116c0568fdc.md) — 동일 끝점 PR 병합 확인
- [claude/nevpt2-analytic-gradient-20260817 · 45aac8f4](branches/45aac8f4611f.md) — 고유 변경·반영 여부 검토
- [claude/sa-casscf-analytic-gradient-20260817 · de565b60](branches/de565b60a805.md) — 고유 변경·반영 여부 검토
- [claude/sa-casscf-zvector-20260816 · 95274875](branches/95274875095c.md) — 고유 변경·반영 여부 검토
- [codex/casscf-analytic-gradient-review-v2-20260816 · 3bdfb60d](branches/3bdfb60dbd9c.md) — 동일 끝점 PR 병합 확인
- [codex/casscf-numgrad-readme-20260815 · ad2fe53a](branches/ad2fe53a0d2c.md) — 동일 끝점 PR 병합 확인
- [codex/mp2-analytic-gradient-20260816 · e7082ac0](branches/e7082ac056bb.md) — 고유 변경·반영 여부 검토
- [codex/mp2-analytic-gradient-specialist-20260816 · 8b7ad53e](branches/8b7ad53e5d7e.md) — 동일 끝점 PR 병합 확인
- [codex/sa-casscf-zvector-pr-20260817 · 076db6cf](branches/076db6cfe604.md) — 동일 끝점 PR 병합 확인
- [codex/xms-caspt2-molcas-benchmark-20260812 · 9bd4d64f](branches/9bd4d64f5880.md) — 고유 변경·반영 여부 검토
- [codex/xms-caspt2-molcas-upstream-20260812 · 28a6e2f4](branches/28a6e2f41eee.md) — main 이력에 포함
- [feat/mp2 · e550dd0b](branches/e550dd0be8df.md) — 고유 변경·반영 여부 검토
- [feature/ncsf · 40790ca6](branches/40790ca63077.md) — 고유 변경·반영 여부 검토
- [feature/ncsf_test · e615a960](branches/e615a9609d89.md) — 고유 변경·반영 여부 검토
- [fix/readme-sa-casscf-analytic · c1f378d9](branches/c1f378d9786d.md) — 동일 끝점 PR 병합 확인
- [fix/rohf-mp2-singles · 17539632](branches/1753963221e4.md) — 동일 끝점 PR 병합 확인
- [fix/scnevpt2-route-contract · 5f7dc646](branches/5f7dc6468bf7.md) — 동일 끝점 PR 병합 확인
- [openqp-fci-option · 51a75c8c](branches/51a75c8c60d5.md) — 고유 변경·반영 여부 검토
- [pr252-mp2-pythonic · 5d7b4cc0](branches/5d7b4cc01cef.md) — 동일 끝점 PR 병합 확인
- [quantum-computing · 04ec3064](branches/04ec3064cef1.md) — 고유 변경·반영 여부 검토
- [quantum-computing-upstream · 03a23cdf](branches/03a23cdffe4b.md) — 동일 끝점 PR 병합 확인
- [wf-native-stack · 6580cba4](branches/6580cba440b8.md) — 동일 끝점 PR 병합 확인
### Solvation / PCM (14)

- [backup/ddx-prerebase · d086594f](branches/d086594f6ec1.md) — 고유 변경·반영 여부 검토
- [backup/local-pcm-before-private-merge · 9d3026ec](branches/9d3026ecdaf1.md) — 고유 변경·반영 여부 검토
- [backup/private-solvent-backend-spike · 1f5b173e](branches/1f5b173eece0.md) — 고유 변경·반영 여부 검토
- [claude/adoring-cerf-zHrke · 99bc32d4](branches/99bc32d4839f.md) — 고유 변경·반영 여부 검토
- [codex/ddx-lifetime-deferred-20260920 · 04846864](branches/04846864d4fe.md) — 고유 변경·반영 여부 검토
- [codex/ddx-linux-libdir · 12ac0f3a](branches/12ac0f3ae8e0.md) — 고유 변경·반영 여부 검토
- [codex/ddx-windows-runtime · 41bbc46b](branches/41bbc46be4a9.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-pcm-spike · 89b5e312](branches/89b5e312ee08.md) — 고유 변경·반영 여부 검토
- [feat/solvent-backend-spike · 85c4caed](branches/85c4caed3ea4.md) — 동일 끝점 PR 병합 확인
- [feat/solvent-backend-spike · 49b717fe](branches/49b717fe59a1.md) — 고유 변경·반영 여부 검토
- [feat/solvent-backend-spike · 7e93ca83](branches/7e93ca83af3e.md) — 고유 변경·반영 여부 검토
- [feat/solvent-backend-spike-private-integrated · 903f0eb4](branches/903f0eb4cd8d.md) — 고유 변경·반영 여부 검토
- [fix/ddx-ci-enable · 5ffa8c06](branches/5ffa8c06aba8.md) — 고유 변경·반영 여부 검토
- [fix/ddx-linux-link-and-ci-enable · ed6c0905](branches/ed6c0905e551.md) — 동일 끝점 PR 병합 확인
### DFTB / xTB / QM/MM (38)

- [agent/dftb-input-logging · 34608427](branches/346084274e50.md) — main 이력에 포함
- [agent/finite-water-droplet-boundary · 515f4f78](branches/515f4f781941.md) — main 이력에 포함
- [agent/finite-water-droplet-boundary · f5c755db](branches/f5c755db59a2.md) — main 이력에 포함
- [agent/omit-duplicate-qmmm-pdb · 241c0a59](branches/241c0a5994c6.md) — 동일 끝점 PR 병합 확인
- [agent/openqp-dftb-abi-fallback · 6820c4b3](branches/6820c4b314e3.md) — 고유 변경·반영 여부 검토
- [archive/dtcam-tb-interface-bridge · 4586740f](branches/4586740f956f.md) — 고유 변경·반영 여부 검토
- [claude/dtcam-tb-interface · 9e5d7240](branches/9e5d72405f93.md) — 동일 끝점 PR 병합 확인
- [claude/dtcam-tb-interface · 55734248](branches/55734248e0ca.md) — 동일 끝점 PR 병합 확인
- [claude/espf-boundary-swscale-20260821 · 26a2e0b1](branches/26a2e0b185f1.md) — 동일 끝점 PR 병합 확인
- [claude/espf-boundary-swscale-20260821 · 8b8a123c](branches/8b8a123c18e4.md) — 고유 변경·반영 여부 검토
- [codex/dftb-s0-continuity · 57872b5f](branches/57872b5f6248.md) — 동일 끝점 PR 병합 확인
- [codex/openqp-dftb-cmake-hook · 2b086076](branches/2b08607666ca.md) — 고유 변경·반영 여부 검토
- [codex/openqp-dftb-native · ad1532ed](branches/ad1532ed27bf.md) — 동일 끝점 PR 병합 확인
- [codex/openqp-dftb-native · 929ea903](branches/929ea903b7a9.md) — 고유 변경·반영 여부 검토
- [codex/openqp-dftb-upstream · 8c30c85d](branches/8c30c85d68a1.md) — main 이력에 포함
- [docs/reduce-dftb-readme · 29c57fc3](branches/29c57fc39252.md) — 동일 끝점 PR 병합 확인
- [feat/dftb-conventional-default · db391506](branches/db391506c820.md) — main 이력에 포함
- [feat/dftb-default-conventional-method · 7ea8fe04](branches/7ea8fe04ebf1.md) — 고유 변경·반영 여부 검토
- [feat/dftb-log-overhaul · 458812c2](branches/458812c239c6.md) — main 이력에 포함
- [feat/dftb-mo-printout · dec2eb79](branches/dec2eb791902.md) — main 이력에 포함
- [feat/dftb-reference-keyword · 44ec1d97](branches/44ec1d97f62d.md) — 동일 끝점 PR 병합 확인
- [feat/dftb-state-configurations · 65d9461e](branches/65d9461e8735.md) — main 이력에 포함
- [feat/dftbplus-excited-state-features · 4cb5c748](branches/4cb5c7488d24.md) — 고유 변경·반영 여부 검토
- [feat/dftbplus-external-backend · 19bbbd5c](branches/19bbbd5cbacf.md) — 고유 변경·반영 여부 검토
- [feat/dftbplus-mrsf-tddftb-private · c702f586](branches/c702f5866e84.md) — 고유 변경·반영 여부 검토
- [feat/mrsf-tddft-spawning-private · 74052b90](branches/74052b90b381.md) — 고유 변경·반영 여부 검토
- [feat/openqp-xtb-adapter · 022016fd](branches/022016fdecf8.md) — 고유 변경·반영 여부 검토
- [feat/openqp-xtb-adapter · aa3865d9](branches/aa3865d9d1e3.md) — 동일 끝점 PR 병합 확인
- [feat/qmmm-frontier-charge · ac84aa50](branches/ac84aa50c37e.md) — 동일 끝점 PR 병합 확인
- [feat/xtb-qmmm-fullespf · b217f46d](branches/b217f46d5785.md) — 고유 변경·반영 여부 검토
- [feat/xtb-qmmm-fullespf · a7fac8df](branches/a7fac8dfd547.md) — 고유 변경·반영 여부 검토
- [fix/qmmm-linkatom-energy-conservation · dd233309](branches/dd2333091a62.md) — 고유 변경·반영 여부 검토
- [fix/qmmm-mechanical-embedding · 575f3d46](branches/575f3d4676fe.md) — 고유 변경·반영 여부 검토
- [fix/qmmm-qm-atoms-order · 561aad6e](branches/561aad6ebf03.md) — 동일 끝점 PR 병합 확인
- [fix/skip-missing-dftb-tests · 3ac04a53](branches/3ac04a537b14.md) — 동일 끝점 PR 병합 확인
- [qmmm · 28b10935](branches/28b109351730.md) — 고유 변경·반영 여부 검토
- [wt/pr354 · 2ea62051](branches/2ea620513530.md) — 고유 변경·반영 여부 검토
- [wt/pr368-solo · d0d499fb](branches/d0d499fbb371.md) — 고유 변경·반영 여부 검토
### Geometry optimization (21)

- [codex/concise-geometry-keywords · e085213c](branches/e085213c6172.md) — 동일 끝점 PR 병합 확인
- [codex/dlc-native-default · 9b8a5400](branches/9b8a5400a5de.md) — 동일 끝점 PR 병합 확인
- [codex/dlc-native-default · b88b7fb0](branches/b88b7fb0e8e9.md) — 고유 변경·반영 여부 검토
- [codex/native-opt-recovery-oqp-format · 26425f58](branches/26425f580a38.md) — main 이력에 포함
- [codex/optimizer-nonfinite-step-recovery · f093f8a4](branches/f093f8a4fa18.md) — main 이력에 포함
- [feat/builtin-optimizer · 4c0cdbcb](branches/4c0cdbcbeaf1.md) — 동일 끝점 PR 병합 확인
- [feat/builtin-optimizer · cf1f3949](branches/cf1f3949163c.md) — 고유 변경·반영 여부 검토
- [feat/geometric-default · 92c055e5](branches/92c055e5ef85.md) — 고유 변경·반영 여부 검토
- [feat/geometric-integration · 60c42ba7](branches/60c42ba71862.md) — main 이력에 포함
- [feat/geometric-neb · b42f9ab5](branches/b42f9ab5bc15.md) — 동일 끝점 PR 병합 확인
- [feat/geometric-optimizer · 03e73fad](branches/03e73fadad82.md) — 동일 끝점 PR 병합 확인
- [feat/grad-cutoff-default-on · c704bfbc](branches/c704bfbcd775.md) — 고유 변경·반영 여부 검토
- [feat/meci-sqp · aef07199](branches/aef07199a6e0.md) — 고유 변경·반영 여부 검토
- [feat/mecp-sqp · a89c6a70](branches/a89c6a70af39.md) — 동일 끝점 PR 병합 확인
- [fix/mecp-converging-objectives · b019d4ca](branches/b019d4caa921.md) — 동일 끝점 PR 병합 확인
- [fix/optimizer-coord-degeneracy · 1ce80760](branches/1ce80760eb26.md) — 고유 변경·반영 여부 검토
- [fix/optimizer-internal-coords · 0d87fa26](branches/0d87fa2606f9.md) — 동일 끝점 PR 병합 확인
- [fix/optimizer-internal-coords · c6134ab8](branches/c6134ab861a1.md) — 고유 변경·반영 여부 검토
- [fix/post-v1.1-test-failures · a3575f09](branches/a3575f095363.md) — 고유 변경·반영 여부 검토
- [fix/post-v1.1-test-failures · 8a1edfbe](branches/8a1edfbe4cd5.md) — 고유 변경·반영 여부 검토
- [pr-142 · 5bc362ee](branches/5bc362eefbe9.md) — 고유 변경·반영 여부 검토
### SCF / response (66)

- [agent/replace-nlopt-simplex-qp · 2d357a50](branches/2d357a509281.md) — 동일 끝점 PR 병합 확인
- [claude/Davidson-ZVector · 582172d9](branches/582172d9920e.md) — 동일 끝점 PR 병합 확인
- [claude/Davidson-ZVector · d8c989ce](branches/d8c989cecf3a.md) — 고유 변경·반영 여부 검토
- [claude/cool-ritchie-vaMB6-pyscf-tools · c4cc4a1c](branches/c4cc4a1c5a37.md) — 고유 변경·반영 여부 검토
- [claude/fix-mrsf-rohf-davidson-diagonal-20260817 · fa52fd31](branches/fa52fd31229c.md) — 동일 끝점 PR 병합 확인
- [claude/mrsf-davidson-window-20260822 · 4091a6bb](branches/4091a6bb742a.md) — 고유 변경·반영 여부 검토
- [claude/mrsf-davidson-window-rebased-20260822 · 972f9990](branches/972f999040f4.md) — 동일 끝점 PR 병합 확인
- [claude/mrsf-davidson-window-rebased-20260822 · a31ae4cd](branches/a31ae4cd4c57.md) — 고유 변경·반영 여부 검토
- [claude/mrsf-davidson-window-v2-20260822 · 268c5727](branches/268c5727052e.md) — 고유 변경·반영 여부 검토
- [claude/mrsf-response-phase · 636ef11a](branches/636ef11ae099.md) — 고유 변경·반영 여부 검토
- [claude/mrsf-response-phase-v2 · 2026aa80](branches/2026aa808bf4.md) — 동일 끝점 PR 병합 확인
- [claude/response-phase-convention · 010571d3](branches/010571d340cd.md) — 고유 변경·반영 여부 검토
- [claude/response-phase-convention2 · 563184fd](branches/563184fd4446.md) — 고유 변경·반영 여부 검토
- [claude/response-phase-convention2 · f6b8fc99](branches/f6b8fc99e2c7.md) — 고유 변경·반영 여부 검토
- [claude/sharp-bell-r7piU · 80936c61](branches/80936c61090f.md) — 고유 변경·반영 여부 검토
- [claude/spin-label-variance · 6e731654](branches/6e7316542a0c.md) — 동일 끝점 PR 병합 확인
- [codex/issue-171-rstctmo-trah · a698dda8](branches/a698dda8d020.md) — 동일 끝점 PR 병합 확인
- [codex/issue-268-json-layout · ffc50acc](branches/ffc50acc6488.md) — 동일 끝점 PR 병합 확인
- [codex/move-scf-converger-ml · 17cad16b](branches/17cad16bbe34.md) — 동일 끝점 PR 병합 확인
- [codex/move-scf-converger-ml · 6b9a86c5](branches/6b9a86c58ed9.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-fock-basis-fix-20260901 · 6f1e569c](branches/6f1e569c7a81.md) — 고유 변경·반영 여부 검토
- [codex/native-trah-default · 6bb207cc](branches/6bb207ccee74.md) — 동일 끝점 PR 병합 확인
- [codex/native-trah-default · 5f87db07](branches/5f87db071b41.md) — non-merge patch 일치·merge 검토
- [codex/native-trah-first · 18f84d67](branches/18f84d67bc32.md) — 고유 변경·반영 여부 검토
- [codex/thymine-ic0103-scf-20260920 · b6f3ebf5](branches/b6f3ebf5cdfa.md) — main 이력에 포함
- [codex/trah-strict-dt-20260914-01a07f50 · 4ec64a5d](branches/4ec64a5d3ffd.md) — 고유 변경·반영 여부 검토
- [codex/zero-beta-native-guesses-upstream-20260813 · 382b8ef0](branches/382b8ef00b70.md) — main 이력에 포함
- [feat/advanced-guess · ef1ffa41](branches/ef1ffa415ba6.md) — 동일 끝점 PR 병합 확인
- [feat/advanced-guess · a3ec8907](branches/a3ec89074053.md) — 고유 변경·반영 여부 검토
- [feat/progressive-screening-scf · 84b19bec](branches/84b19bec2959.md) — 동일 끝점 PR 병합 확인
- [feat/scf-cleanup · 3ba40812](branches/3ba40812f00f.md) — 고유 변경·반영 여부 검토
- [feat/scf-cleanup · 8fb04a2c](branches/8fb04a2caede.md) — 고유 변경·반영 여부 검토
- [fix/ecp-cffi-buffer-lifetime · 0c3c29a2](branches/0c3c29a29f91.md) — 동일 끝점 PR 병합 확인
- [fix/ecp-cffi-buffer-lifetime · bbd581ba](branches/bbd581ba16b3.md) — 고유 변경·반영 여부 검토
- [fix/fcheck-bounds-guess-dftlib · 3245c0e5](branches/3245c0e56a52.md) — 고유 변경·반영 여부 검토
- [fix/formaldehyde-gradient-highroot · 2a0caadc](branches/2a0caadc37c9.md) — 고유 변경·반영 여부 검토
- [fix/higher-root-gradient-zvector · a475b429](branches/a475b429ed05.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-davidson-symmetry-block-coverage · 88184e7b](branches/88184e7b78c2.md) — 동일 끝점 PR 병합 확인
- [fix/mrsf-gradient-remaining-highroot-diagnostics · a8f0fd59](branches/a8f0fd598a76.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-gradient-rohf-mob-trial · 86295c92](branches/86295c924d4c.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-h2-zero-closed · f3a40e98](branches/f3a40e988d56.md) — 동일 끝점 PR 병합 확인
- [fix/mrsf-h2-zero-closed · ee4ce62f](branches/ee4ce62fd98f.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-h2-zero-closed · 8733067e](branches/8733067effb5.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-ovov-gradient-sign · 8c23f15d](branches/8c23f15dd92c.md) — 동일 끝점 PR 병합 확인
- [fix/mrsf-ovov-gradient-sign-diagnostic · 47914f5b](branches/47914f5be800.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-reference-scf-stability · 6ac212a0](branches/6ac212a054ee.md) — 동일 끝점 PR 병합 확인
- [fix/mrsf-reference-scf-stability · f70d72f8](branches/f70d72f8d3ce.md) — 고유 변경·반영 여부 검토
- [fix/mrsf-zvector-lhs · a7b8e2f5](branches/a7b8e2f579c3.md) — 고유 변경·반영 여부 검토
- [fix/scf-stability-opt-in · 90bdea5b](branches/90bdea5bf6a8.md) — 동일 끝점 PR 병합 확인
- [fix/tddft-higher-root-gradient · 78502154](branches/78502154ca67.md) — non-merge patch 일치·merge 검토
- [fix/tddft-higher-root-gradient · 6fd72e34](branches/6fd72e34b335.md) — 동일 끝점 PR 병합 확인
- [fix/trah-rstctmo-handoff · 11dd37f7](branches/11dd37f7f359.md) — main 이력에 포함
- [main · 6ab4b867](branches/6ab4b86701bf.md) — main 이력에 포함
- [main · 452af2e5](branches/452af2e5cc0c.md) — main 이력에 포함
- [mrsf_z_vector · 955c13f7](branches/955c13f75a95.md) — 고유 변경·반영 여부 검토
- [openqp-zvector-davidson-stability · f601f83e](branches/f601f83efb12.md) — 고유 변경·반영 여부 검토
- [perf/hf_dft_gradient · 31523f3d](branches/31523f3dafd8.md) — 동일 끝점 PR 병합 확인
- [perf/hf_dft_gradient · 3d6cf9b0](branches/3d6cf9b0f291.md) — 고유 변경·반영 여부 검토
- [perf/mrsf-fock-digestion · fc868854](branches/fc868854932d.md) — 동일 끝점 PR 병합 확인
- [perf/mrsf-zvector · 6f679446](branches/6f6794468a86.md) — 동일 끝점 PR 병합 확인
- [perf/tdhf-davidson-timers · 6064b5b4](branches/6064b5b4f9a5.md) — 고유 변경·반영 여부 검토
- [pr-172 · 0029fd6b](branches/0029fd6b3cb5.md) — 고유 변경·반영 여부 검토
- [refactor/huckel-guess · e3c012f5](branches/e3c012f53720.md) — 동일 끝점 PR 병합 확인
- [revert-63-main · 0a2023a8](branches/0a2023a8f024.md) — 동일 끝점 PR 병합 확인
- [wt/pr363-ci-fix · cc222b55](branches/cc222b55cc06.md) — 고유 변경·반영 여부 검토
- [z-vector · 05cfd48a](branches/05cfd48abc30.md) — 고유 변경·반영 여부 검토
### DFT / grid (15)

- [agent/fix-sg1-radial-grid · 4a9aea0a](branches/4a9aea0a42ba.md) — 고유 변경·반영 여부 검토
- [claude/gallant-ramanujan-aJf6x · 20f82940](branches/20f8294073bc.md) — 고유 변경·반영 여부 검토
- [claude/gallant-ramanujan-aJf6x · b8c5f0a1](branches/b8c5f0a12eff.md) — 고유 변경·반영 여부 검토
- [codex/nc-mrsf-kernel-sign-scale · 9312426d](branches/9312426d459f.md) — 고유 변경·반영 여부 검토
- [codex/xc-grid-gradient-20260828 · a4009ae4](branches/a4009ae4aa8f.md) — 동일 끝점 PR 병합 확인
- [draft/xc-orbit-weight · cdaa1790](branches/cdaa179083c9.md) — 고유 변경·반영 여부 검토
- [feat/coarse-to-fine-xc-grid · d749e29d](branches/d749e29deb37.md) — 동일 끝점 PR 병합 확인
- [feat/dft-grid-default · c60285ed](branches/c60285edcb0d.md) — 동일 끝점 PR 병합 확인
- [feat/dft-xc-grid-reuse · e931030b](branches/e931030b42db.md) — 고유 변경·반영 여부 검토
- [feat/dft-xc-grid-reuse-upstream · 5e106e56](branches/5e106e567fe1.md) — 동일 끝점 PR 병합 확인
- [feature/sg-grids · c62fa90e](branches/c62fa90ee31c.md) — 동일 끝점 PR 병합 확인
- [fix/sg1-grid · 954342ae](branches/954342ae8c54.md) — 동일 끝점 PR 병합 확인
- [perf/tdhf-xc-response-cache · 7a704d10](branches/7a704d103882.md) — 고유 변경·반영 여부 검토
- [perf/xc-numerical-kernel · ed06cd99](branches/ed06cd997fdd.md) — 동일 끝점 PR 병합 확인
- [resolve-241 · 868389cf](branches/868389cf4dd8.md) — 고유 변경·반영 여부 검토
### Symmetry (11)

- [codex/molecular-symmetry-validate · 7c9bc92b](branches/7c9bc92b2bb8.md) — 고유 변경·반영 여부 검토
- [codex/molecular-symmetry-validate · ee20e381](branches/ee20e381152c.md) — 고유 변경·반영 여부 검토
- [draft/grd2-petite-optin · c9af7ecf](branches/c9af7ecfe49a.md) — 동일 끝점 PR 병합 확인
- [draft/inivec-irrep-coverage · 120c10aa](branches/120c10aad6f0.md) — 고유 변경·반영 여부 검토
- [draft/thermo-symmetry-number · 1cfc745d](branches/1cfc745df1ff.md) — 동일 끝점 PR 병합 확인
- [feat/molecular-symmetry · fc4388ad](branches/fc4388ad475a.md) — 고유 변경·반영 여부 검토
- [feat/molecular-symmetry · 6103c85c](branches/6103c85c70e3.md) — 고유 변경·반영 여부 검토
- [feat/petite-no-reorient · c90c8a9e](branches/c90c8a9e5192.md) — 동일 끝점 PR 병합 확인
- [feat/symmetry-default-on · ef745b72](branches/ef745b72f75d.md) — 동일 끝점 PR 병합 확인
- [feat/wake-petite-smoke-tests · 3e6ee0c8](branches/3e6ee0c8119b.md) — 고유 변경·반영 여부 검토
- [pr199 · 16e35679](branches/16e356791566.md) — 고유 변경·반영 여부 검토
### Build / CI / release (69)

- [agent/dftd4-shared-libraries · eec27707](branches/eec27707c79c.md) — 동일 끝점 PR 병합 확인
- [agent/docker-release-hardening · 767e654a](branches/767e654a3638.md) — 동일 끝점 PR 병합 확인
- [agent/integrate-afqmc-build · e64d2931](branches/e64d29316200.md) — 고유 변경·반영 여부 검토
- [agent/v1.3-license-governance · 290780e1](branches/290780e1a6c7.md) — 동일 끝점 PR 병합 확인
- [agent/v130-exclude-private-backends-20260813 · 41c72d2e](branches/41c72d2e9ba7.md) — 동일 끝점 PR 병합 확인
- [backup-windows-prerebase · 0511f7fd](branches/0511f7fd957e.md) — 고유 변경·반영 여부 검토
- [blas-platform-autoselect · 8edb48df](branches/8edb48dfabed.md) — 고유 변경·반영 여부 검토
- [blas-platform-autoselect · bea2aef2](branches/bea2aef20b32.md) — 고유 변경·반영 여부 검토
- [build/docker-openblas-ilp64 · 450644a1](branches/450644a130d5.md) — 동일 끝점 PR 병합 확인
- [build/openblas-ilp64-presets · 09de7dee](branches/09de7deef3c3.md) — 고유 변경·반영 여부 검토
- [chore/allowlist-digest-gate · 2cd8b847](branches/2cd8b84754d6.md) — main 이력에 포함
- [chore/remove-poc-benchmark · 376a13e7](branches/376a13e7d8aa.md) — 동일 끝점 PR 병합 확인
- [chore/strip-dev-artifacts · 6a7dea0d](branches/6a7dea0daefd.md) — 동일 끝점 PR 병합 확인
- [chore/strip-dev-artifacts · 666aa200](branches/666aa200bc9f.md) — 고유 변경·반영 여부 검토
- [ci/run-pytest-suite · 75fa7b9d](branches/75fa7b9d7f51.md) — 동일 끝점 PR 병합 확인
- [ci/wheel-openblas-no-lapack · 767df868](branches/767df8682b81.md) — 동일 끝점 PR 병합 확인
- [cl · 1686aecb](branches/1686aecbaca4.md) — 고유 변경·반영 여부 검토
- [claude-manager · 1890fefd](branches/1890fefdb25d.md) — main 이력에 포함
- [claude/amazing-ramanujan-b881b8 · f2a24c1c](branches/f2a24c1ce970.md) — main 이력에 포함
- [claude/ci-concurrency-group · 05c2056f](branches/05c2056f5283.md) — 고유 변경·반영 여부 검토
- [claude/ci-concurrency-group2 · 827b3876](branches/827b3876a444.md) — 동일 끝점 PR 병합 확인
- [claude/ci-concurrency-group2 · 32da9b46](branches/32da9b467aff.md) — non-merge patch 일치·merge 검토
- [claude/pr-188-review-ZC0GI · 903e686b](branches/903e686ba738.md) — 고유 변경·반영 여부 검토
- [claude/remove-defined-ilp64 · 96334e17](branches/96334e1753ad.md) — 동일 끝점 PR 병합 확인
- [codex/auto-request-mohsen-review · 9e2691ca](branches/9e2691ca39e2.md) — 동일 끝점 PR 병합 확인
- [codex/cheol-ci-codex-review-20260921 · 5668ef51](branches/5668ef511d43.md) — main 이력에 포함
- [codex/cheol-ci-codex-review-20260921 · 6d3f42b2](branches/6d3f42b2279f.md) — 고유 변경·반영 여부 검토
- [codex/claude-all-pr-review · 70a663f0](branches/70a663f0fe01.md) — 동일 끝점 PR 병합 확인
- [codex/docker-ci-cache · 28fe8d66](branches/28fe8d669588.md) — 동일 끝점 PR 병합 확인
- [codex/github-sync-automerge-retry-20260920 · 3035051a](branches/3035051a54dc.md) — main 이력에 포함
- [codex/gitlab-auto-review-20260920 · eb97168c](branches/eb97168c7ad3.md) — main 이력에 포함
- [codex/gitlab-main-migration-20260919 · 48355a72](branches/48355a721b35.md) — main 이력에 포함
- [codex/internal-ci-dedup-20260826 · 6236543c](branches/6236543cdcee.md) — main 이력에 포함
- [codex/internal-gitlab-ci-20260826 · 564c3237](branches/564c3237f0d1.md) — main 이력에 포함
- [codex/macos-blas-lp64-bench · 18829034](branches/18829034b915.md) — 고유 변경·반영 여부 검토
- [codex/macos-blas-lp64-bench-upstream · 759ceaec](branches/759ceaec3d85.md) — 동일 끝점 PR 병합 확인
- [codex/macos-blas-lp64-bench-upstream · d0c91ad5](branches/d0c91ad56ef1.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-final-build-7c73f6fd-20260901 · 7c73f6fd](branches/7c73f6fdb8aa.md) — 고유 변경·반영 여부 검토
- [codex/mrsf-formaldehyde-tpx-build-20260902 · a5ee259b](branches/a5ee259bd140.md) — 고유 변경·반영 여부 검토
- [codex/pr344-review-maintenance-20260816 · 0769abb0](branches/0769abb03d75.md) — 고유 변경·반영 여부 검토
- [codex/pr345-ci-review-maintenance-20260816 · 2a015820](branches/2a0158208d07.md) — 고유 변경·반영 여부 검토
- [codex/reintegrate-devkit-20260921 · 3799e40e](branches/3799e40e2c17.md) — main 이력에 포함
- [codex/release-v1.2.1 · 52507490](branches/525074908253.md) — 동일 끝점 PR 병합 확인
- [codex/reuse-external-build-cache · 3f7da4ad](branches/3f7da4ad90af.md) — 동일 끝점 PR 병합 확인
- [feat/intel-oneapi-support · fa7ae05c](branches/fa7ae05ca407.md) — 동일 끝점 PR 병합 확인
- [feat/regression-registry-and-coverage · 78f78e3c](branches/78f78e3ca98d.md) — 동일 끝점 PR 병합 확인
- [feat/windows-intel-build · 9de4cfbb](branches/9de4cfbbd3cb.md) — 동일 끝점 PR 병합 확인
- [feature/blas-native-policy · dcc49ad2](branches/dcc49ad2f671.md) — 동일 끝점 PR 병합 확인
- [feature/test-on-CI · 90525350](branches/905253505e8e.md) — 고유 변경·반영 여부 검토
- [fix-macos-accelerate-oqp-coords-warnings · b5f72bfb](branches/b5f72bfbfb91.md) — 동일 끝점 PR 병합 확인
- [fix/build-and-docs-followup · fc6a1fe9](branches/fc6a1fe9a728.md) — 동일 끝점 PR 병합 확인
- [fix/claude-review-github-token · e7647d82](branches/e7647d82678f.md) — 동일 끝점 PR 병합 확인
- [fix/claude-review-track-progress · 6e007eb4](branches/6e007eb4ab70.md) — main 이력에 포함
- [fix/lf-license-files-only · b4226229](branches/b42262292d7b.md) — 동일 끝점 PR 병합 확인
- [fix/lf-line-endings · 7a49b678](branches/7a49b678a90a.md) — 동일 끝점 PR 병합 확인
- [fix/macos-wheel-delocate-0.13 · cb5c8509](branches/cb5c8509bbd6.md) — 동일 끝점 PR 병합 확인
- [fix/train-selector-tag-push · 35ee8853](branches/35ee8853e762.md) — main 이력에 포함
- [fix/uhf-stability-and-ci-runtests · 24c8c5c8](branches/24c8c5c8c8a3.md) — 동일 끝점 PR 병합 확인
- [main · d8fcc411](branches/d8fcc4119c73.md) — main 이력에 포함
- [main · 220afefb](branches/220afefb2f72.md) — main 이력에 포함
- [main-preflight · 9d9cc53c](branches/9d9cc53c822b.md) — main 이력에 포함
- [perf/eigen-blas-thread-policy · f0d7f437](branches/f0d7f4374c93.md) — 동일 끝점 PR 병합 확인
- [policy/pr-enforcement · 9c2a8f7e](branches/9c2a8f7e263a.md) — 동일 끝점 PR 병합 확인
- [pr-169 · 0f78a348](branches/0f78a348cf95.md) — 고유 변경·반영 여부 검토
- [pr-250-review · d0446a0c](branches/d0446a0c7fef.md) — 고유 변경·반영 여부 검토
- [pr203-check · fe1440f9](branches/fe1440f9df4d.md) — 고유 변경·반영 여부 검토
- [release-prep · d452be68](branches/d452be6832d4.md) — main 이력에 포함
- [release/v1.3.1 · ba072ed1](branches/ba072ed1f6a6.md) — 동일 끝점 PR 병합 확인
- [update/libtagarray · caa16352](branches/caa1635263c3.md) — 고유 변경·반영 여부 검토
### Docs / API / input (35)

- [agent/clarify-current-gpl · 8f99b590](branches/8f99b59083b9.md) — 동일 끝점 PR 병합 확인
- [agent/export-mo-frequency-formats · 57a70ade](branches/57a70ade8737.md) — 동일 끝점 PR 병합 확인
- [agent/native-openqp-terminology · 8a02e7cd](branches/8a02e7cd2e1e.md) — non-merge patch 일치·merge 검토
- [agent/readable-oqp-input · 2999ee7c](branches/2999ee7c4fdc.md) — 고유 변경·반영 여부 검토
- [chore/split-devkit-and-layout-gate · 7938b613](branches/7938b613ca66.md) — main 이력에 포함
- [claude/all-open-issues-20260821 · 7bc51dbd](branches/7bc51dbd1718.md) — 고유 변경·반영 여부 검토
- [codex/export-mrsf-analysis · 6b35f1ee](branches/6b35f1eecbd0.md) — 동일 끝점 PR 병합 확인
- [codex/list-available-functionalities · 1fd8aeaa](branches/1fd8aeaa99d1.md) — 고유 변경·반영 여부 검토
- [codex/openqp-log-consistency-authors-20260816 · b3854705](branches/b385470512b3.md) — 동일 끝점 PR 병합 확인
- [codex/openqp-log-consistency-authors-20260816 · 0baea2f1](branches/0baea2f17d14.md) — 고유 변경·반영 여부 검토
- [codex/oqp-natural-input · 2d64f666](branches/2d64f6664e50.md) — 동일 끝점 PR 병합 확인
- [codex/pr347-oqup-pronunciation-20260817 · 19335083](branches/193350832e51.md) — 고유 변경·반영 여부 검토
- [codex/readme-foldable-methods-20260829 · 390121dd](branches/390121ddaa55.md) — 동일 끝점 PR 병합 확인
- [codex/readme-foldable-parser-tests-20260829 · f95091ae](branches/f95091ae73ec.md) — 동일 끝점 PR 병합 확인
- [codex/readme-oqp-studio · eec316ed](branches/eec316edae2c.md) — 동일 끝점 PR 병합 확인
- [codex/readme-oqp-studio-mo-20260829 · 6ba58184](branches/6ba58184f2c5.md) — 동일 끝점 PR 병합 확인
- [codex/readme-web-link · 3c8d1eeb](branches/3c8d1eeb0785.md) — main 이력에 포함
- [codex/space-separated-method-basis-20260829 · dff991c2](branches/dff991c29e50.md) — 동일 끝점 PR 병합 확인
- [docs/drop-tight-binding · c21f216b](branches/c21f216b2e66.md) — 동일 끝점 PR 병합 확인
- [docs/readme-ecosystem · c9548bf9](branches/c9548bf9ec41.md) — 고유 변경·반영 여부 검토
- [docs/readme-enhance-upstream · d257246a](branches/d257246aa4a9.md) — 동일 끝점 PR 병합 확인
- [docs/readme-tutorial-deeplinks · c3831638](branches/c3831638cfea.md) — 동일 끝점 PR 병합 확인
- [docs/readme-wiki-updates · 67b379f1](branches/67b379f1ec6b.md) — 동일 끝점 PR 병합 확인
- [docs/simplify-oqp-inputs · ec510664](branches/ec5106648a6a.md) — 동일 끝점 PR 병합 확인
- [draft/input-to-standard · d35754c1](branches/d35754c19199.md) — 고유 변경·반영 여부 검토
- [feat/log-grid-info · 5b9485df](branches/5b9485dfc76c.md) — 동일 끝점 PR 병합 확인
- [feat/pythonic-api-compat · 347e8d9b](branches/347e8d9b15c3.md) — 동일 끝점 PR 병합 확인
- [feat/pythonic-api-compat · 3466f477](branches/3466f477116c.md) — 고유 변경·반영 여부 검토
- [fix/molden-high-angular-momentum · c00ba614](branches/c00ba6147a72.md) — 동일 끝점 PR 병합 확인
- [fix/molden-high-angular-momentum · df8e7266](branches/df8e726632a6.md) — 고유 변경·반영 여부 검토
- [karmachoi-patch-1 · 3f47f2f7](branches/3f47f2f76a57.md) — main 이력에 포함
- [karmachoi/README.md · b378baae](branches/b378baae0337.md) — 별도 이력
- [perf/omp-input · cd957d48](branches/cd957d4831be.md) — 동일 끝점 PR 병합 확인
- [perf/omp-input · df2a2e7f](branches/df2a2e7f4119.md) — 고유 변경·반영 여부 검토
- [python_version_selection · af888b1b](branches/af888b1b24c9.md) — 고유 변경·반영 여부 검토
### Historical documents (4)

- [backup/thymine-before-pr-cleanup-20260917 · 537a4181](branches/537a418104ec.md) — main 이력에 포함
- [gh-pages · d47819b4](branches/d47819b45257.md) — 별도 이력
- [gh-pages · cc8d6bef](branches/cc8d6bef93f3.md) — 별도 이력
- [gh-pages · ffc302f5](branches/ffc302f570af.md) — 별도 이력
### Other / mixed (18)

- [agent/native-openqp-terminology · e8ab06dd](branches/e8ab06dd302a.md) — 동일 끝점 PR 병합 확인
- [claude/lucid-ritchie-aa76a2 · 19381380](branches/193813805d1a.md) — main 이력에 포함
- [codex/fix-fort6-cleanup-20260826 · 293308b7](branches/293308b7c30a.md) — 고유 변경·반영 여부 검토
- [feat/mixedprec-fp32-routing · d5bd9724](branches/d5bd9724e420.md) — main 이력에 포함
- [feat/native-dftd4 · 4a335a32](branches/4a335a3230fe.md) — 고유 변경·반영 여부 검토
- [feat/perf-levels · 365391f6](branches/365391f66c3e.md) — 동일 끝점 PR 병합 확인
- [feat/three-state-ci-benchmark · ba75c728](branches/ba75c7280227.md) — 고유 변경·반영 여부 검토
- [fix/bounds-check-conformance · bd43ba5e](branches/bd43ba5e8126.md) — 동일 끝점 PR 병합 확인
- [fix/fcheck-bounds-aborts · 32a8e121](branches/32a8e1214890.md) — 고유 변경·반영 여부 검토
- [fix/omp-barrier-int1 · 1ccde326](branches/1ccde326aae9.md) — 동일 끝점 PR 병합 확인
- [fix/post-v1.1-test-failures · 87ff6774](branches/87ff677406ff.md) — 동일 끝점 PR 병합 확인
- [fix/states-overlap-ov-exact-diagonal · e5eeed07](branches/e5eeed07f9b9.md) — main 이력에 포함
- [fix/write-xyz-numpy2-regression · 96105f13](branches/96105f13115c.md) — 동일 끝점 PR 병합 확인
- [kk/rm_xints · 8746fd6d](branches/8746fd6d6344.md) — main 이력에 포함
- [main · b9a6ba76](branches/b9a6ba76bede.md) — main 이력에 포함
- [main · f2bc79bf](branches/f2bc79bf3ea8.md) — main 이력에 포함
- [resolve/pr168-conflicts · d002ee3c](branches/d002ee3c6d77.md) — 동일 끝점 PR 병합 확인
- [upstream-main · c7db4d22](branches/c7db4d222a7a.md) — main 이력에 포함

분야 분류는 검색을 돕기 위한 것이다. 115개 이름 불명확 항목에 로컬 Qwen의 구조 검증된 분류를 참고했으며, SHA·개수·병합 판정은 모두 Git/API에서 계산했다. 분류 자체는 과학적 검증이 아니다.