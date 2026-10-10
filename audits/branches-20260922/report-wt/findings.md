# OpenQP 개인·private 개발 브랜치 전수 조사

기준일: 2026-09-22. 조사 대상은 GitHub `karmachoi/openqp`,
`karmachoi/openqp-private`, 현재 GitLab `open-quantum-platform/internal/openqp`,
그리고 Ultra에서 확인한 해당 저장소의 로컬 개발 브랜치다.
GitLab의 과거 `cheol/openqp`는 현재 내부 그룹 저장소이며,
별도 `cheol/openqp-private` 프로젝트는 이번 GitLab 목록에서 발견되지 않았다.

**두 개인 저장소는 단순 중복이 아니다. 개발 이력 대부분은 GitLab에 보존되어
있지만, 아직 GitLab에 없는 최신 REKS 변경과 로컬 변경이 있다. 브랜치 이름이나
마지막 커밋 제목만 보고 삭제·archive하면 안 된다.**

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
