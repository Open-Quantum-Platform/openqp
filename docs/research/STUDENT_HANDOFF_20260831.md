# OpenQP MRSF analytic NAC · NAMD 연구 인수인계

최초 작성: 2026-08-31 12:36 KST · **갱신: 2026-09-07 (NAMD 캠페인 결과 반영)**
중앙 연구 과제: `CHC-20260829-44893D`
이 문서는 학생이 기존 계산을 중복하지 않고, 확인된 결함을 다시 만들지 않으며, 코드·계산·논문을 이어가기 위한 기준 문서이다.
analytic NAC과 NAMD 두 갈래를 한 문서로 합쳐 관리한다.

## 0. 2026-09-07 갱신: 아래 1절의 진단은 대체되었다

8월 31일자 1절은 Uracil의 조기 S0 전이 원인을 `tdc=npi`와 `decoherence=edc`로 지목하고
`tdc=fd`, `tlf=2`, `decoherence=off`로 재현할 것을 지시했다. **이 진단과 처방은 모두 틀렸다.**
9월 5일부터 7일까지의 캠페인에서 확인된 사실은 다음과 같다.

1. **근본 원인은 초기 속도였다.** 2022년 KNU-GAMESS 궤적은 `TVELQM`이 적용되지 않아
   정지 상태에서 출발한다. 영점 에너지가 없다. OpenQP는 Wigner 분포에서 출발하며 평균
   운동 에너지가 1.28 eV이다. 두 계산의 차이 대부분이 이 하나에서 나온다. analytic NAC의
   정확도 문제가 아니다.
2. **`tlf=2`는 처방이 아니라 결함이다.** 축약 Leibniz 전개는 준축퇴 점유 오비탈이 회전할 때
   상태 겹침을 붕괴시킨다. 1절이 지목한 IC49 56.5 fs의 33배 과대 결합이 바로 이 붕괴의
   결과이며, NPI 자체의 결함이 아니다. 정확한 소행렬식(`tlf=0`)이 이제 기본값이고 비용도
   같다. `tlf=2`로 새 계산을 하지 말 것.
3. **DIIS 연속화가 이중 점유 오비탈과 SOMO를 교환하는 경우가 있다.** 이제 검출하고
   `ref_follow`가 SOSCF와 TRAH로 단계적으로 대응하며, 실패하면 Hückel 추정으로 되돌린다.
4. **에너지 표류는 dt 0.5 fs의 적분 오차였다.** 가파른 S1 구간에서 스텝을 놓친다.
   불연속이 감지되면 그 스텝을 부분 스텝으로 다시 밟는 `disc_substeps`가 이를 없앴다.
   표류 중앙값이 +3~8 eV에서 −0.01 eV로 내려갔다.

## 0.05 앞으로의 기본 설정: TDC+NAC + 운동량 반전

**새 MRSF FSSH 계산은 TDC+NAC에 좌절 도약 운동량 반전을 켜고, 간격 컷오프 없이 돌린다.**
캠페인 결과 디렉터리로는 `ht_nac_nogate_reflect`에 해당한다. 입력의 핵심 키워드는 다음과 같다.

```
tdc=npi                 # 겹침으로 시간 미분 결합을 전파
rescale=hop_analytic_nac # 도약 쌍에 대해서만 해석적 d를 계산해 그 방향으로 재조정
frustrated=reflect      # 좌절 도약에서 d 방향 운동량 성분을 반전
thrshe=0.36749322175655 # 간격 컷오프 사실상 해제 (10 eV)
decoherence=edc, edc_c=0.1
tlf=0                   # 정확한 소행렬식. tlf=2를 쓰지 말 것
mo_reuse=true, ref_follow=soscf, ref_switch_rescale=true
disc_rescale=true, disc_tol=0.002, disc_substeps=10
substep=50000, dt 0.5 fs, NVE
```

근거는 세 가지다.

1. **해석적 결합이 유일하게 의미 있는 요인이다.** 겹침만으로 결합을 추정하면 150 fs 이후
   하향 도약이 54회로 폭증하고 바닥상태 점유율이 0.72까지 올라간다. 도약 쌍에 해석적 d를
   넣으면 궤적 50개 중 41개가 달라지고 0.58로 내려가며, 150 fs 이후 거의 변하지 않는
   GAMESS 기록의 모양에 가까워진다.
2. **컷오프는 겹침 결합의 증상을 가리는 장치다.** TDC+NAC에서는 궤적 50개 중 42개를 그대로
   두어 사실상 아무 일도 하지 않는다. 결합을 제대로 계산하면 큰 간격에서 도약 후보가
   애초에 생기지 않으므로 컷오프로 막을 것이 없다. 증상을 가리는 대신 원인을 고친다.
3. **반전은 점유율을 바꾸지 않지만 켜 두는 것이 옳다.** 300 fs 점유율 변화가 정확히
   0.00이고 궤적 14개만 개별적으로 달라진다. 그러나 좌절 도약에서 d 방향 성분을 반전하는
   것이 표준적인 처방이고 운동 에너지를 보존하며 추가 비용이 없다. 결과를 바꾸지 않으면서
   이론적으로 방어 가능한 선택이므로 기본값으로 둔다. 반대로 반전을 껐을 때 결과가
   달라진다면 그 계산은 좌절 도약 처리에 민감하다는 뜻이므로 따로 보고해야 한다.

Full NAC(전 쌍 해석 결합)은 비용이 훨씬 크다. 궤적당 맥미니 14코어로 약 7.5시간,
TDC+NAC의 여러 배다. 앙상블 결과가 나오면 TDC+NAC과 비교해 추가 비용이 정당한지
판단하되, 그때까지 기본은 TDC+NAC + 반전이다.

## 0.1 NAMD 처리 비교 캠페인 (2026-09-05 ~ )

일곱 가지 처리를 각각 Wigner 초기조건 50개, 300 fs, dt 0.5 fs, NVE,
MRSF-TDDFT/BH&HLYP/6-31G*, EDC 0.1로 계산한다. 일곱 앙상블이 **같은 초기조건과 같은
궤적별 난수열**을 쓰므로 평균뿐 아니라 궤적 단위로 짝지어 비교할 수 있다.

| 결과 디렉터리 | 처리 | 결합 | 재조정 | 간격 컷오프 | 좌절 도약 |
| --- | --- | --- | --- | --- | --- |
| `td_npi_iso` | Overlap TDC, cutoff | 겹침(NPI) | 등방 | 10 kcal/mol | 무시 |
| `td_npi_iso_nocut` | Overlap TDC | 겹침(NPI) | 등방 | 없음 | 무시 |
| `baeck_an_iso` | Curvature TDC, cutoff | Baeck–An | 등방 | 10 kcal/mol | 무시 |
| `ht_nac` | TDC+NAC, cutoff | 겹침 + 도약쌍 해석 d | d 방향 | 10 kcal/mol | 무시 |
| `ht_nac_nogate` | TDC+NAC | 위와 같음 | d 방향 | 없음 | 무시 |
| `ht_nac_nogate_reflect` | TDC+NAC, reversal | 위와 같음 | d 방향 | 없음 | d 따라 반전 |
| `full_analytic_nogate_reflect` | Full NAC | 전 쌍 해석 d·v | d 방향 | 없음 | d 따라 반전 |

컷오프는 `thrshe` 키워드이다. 0.0159362 Hartree(10 kcal/mol)이면 켜짐,
0.36749 Hartree(10 eV)이면 사실상 꺼짐이다.

### 바닥상태 점유율 (완성 앙상블만, n=50)

| 처리 | 100 fs | 150 fs | 300 fs |
| --- | --- | --- | --- |
| Overlap TDC | 0.42 | 0.52 | 0.72 |
| Overlap TDC, cutoff | 0.50 | 0.56 | 0.62 |
| Curvature TDC, cutoff | 0.40 | 0.50 | 0.52 |
| TDC+NAC, cutoff | 0.44 | 0.50 | 0.56 |
| TDC+NAC | 0.42 | 0.52 | 0.58 |
| TDC+NAC, reversal | 0.48 | 0.52 | 0.58 |
| Park 2022 GAMESS (v=0) | 0.36 | 0.40 | 0.41 |

n=50에서 부트스트랩 구간이 ±0.13이므로 개별 값끼리는 겹친다. 짝지은 비교를 볼 것.

### 요인별 효과 (같은 초기조건끼리 짝지음)

| 요인 | 비교 | 동일 궤적 | 300 fs 변화 |
| --- | --- | --- | --- |
| 결합 모델 | 겹침 → 겹침 + 해석 NAC | 9/50 | −0.14 |
| 컷오프 | Overlap TDC, 켬 → 끔 | 27/49 | +0.10 |
| 컷오프 | TDC+NAC, 켬 → 끔 | 42/50 | +0.02 |
| 운동량 반전 | TDC+NAC → 반전 추가 | 36/50 | 0.00 |

**결합 모델만이 앙상블을 재편한다.** 여섯 처리가 150 fs까지 0.02 안에서 일치하고 그
이후에만 갈라진다. 150 fs 이후 겹침 전용 처리는 하향 도약을 54회 일으키는 반면 나머지는
4~15회이며, 상향 도약도 궤적당 2.12회 대 0.55회다. 궤적이 오래 왕복하다 바닥상태에
갇힌다. 컷오프는 큰 간격에서 도약을 제안하는 유일한 추정인 겹침에서만 의미가 있다.
운동량 반전은 궤적 14개를 바꾸지만 점유율은 정확히 0만큼 바꾼다.

**결어긋남은 시험되지 않았다.** 일곱 앙상블 모두 EDC 0.1이고 대조군이 없다. PyRAI2MD
레인만 결어긋남이 꺼져 있으나 결합·재조정·좌절 도약 처리가 모두 달라 요인을 분리할 수
없다. 분리하려면 `ht_nac_nogate`에 `decoherence=off`만 바꾼 50궤적 앙상블이 하나 더
필요하며, 현재 자원으로 하루면 된다.

### 진행 상황과 경로

| 레인 | 상태 |
| --- | --- |
| 여섯 앙상블 | 완료 (n=50) |
| Full NAC | 39/50, 오늘 밤 완료 예상 |
| PyRAI2MD 300 fs | 31/50, 하루 더 |

- 세션 기록: `/Volumes/OpenQP.DB/claude/sessions/20260905_uracil_namd_nac_check/`
  (`REPORT_uracil_namd_nac_check.md`에 전체 경위, `analysis/`에 분석, `fig8_campaign/`에 캠페인 스크립트)
- 클러스터 캠페인 루트: `/bighome/cheolho.choi/openqp-uracil-fig8-dd3275d-20260905/`
  (`results/prod300`은 수정 전 폐기 단계이므로 쓰지 말 것)
- PyRAI2MD: `/bighome/cheolho.choi/openqp-uracil-namd-moreuse-1ff8139-20260905/pyrai2md/`
- 논문: Overleaf `6a6e9cf68a2d15fed156be72`. Figure 8은 완성된 여섯 앙상블을 담고 있다.
  `namd_results.tex` 본문은 아직 옛 이름과 옛 수치가 남아 있어 Full NAC 완료 후 정리해야 한다.

### 마무리 절차

1. `fig8_campaign/scripts/audit_completeness.sh` — 이상이 없으면 앙상블별 개수만 출력한다.
2. `analysis/sync_results.sh` — 완주 궤적만 세션 트리로 가져온다.
3. `analysis/final_analysis.py` — 점유율·도약 통계·에너지 표류 표를 만든다.
   50개 미만 앙상블은 의도적으로 표에서 제외한다.
4. `analysis/fig8_complete.py` — 완성 앙상블만으로 Figure 8을 다시 만든다.
5. `analysis/factor_effects.py` — 위 요인별 효과를 재현한다.

### 반드시 피해야 할 함정

- **부분 앙상블의 점유율을 인용하지 말 것.** 벽시계 시간이 활성 상태와 상관된다(r=−0.85).
  S0에 일찍 떨어진 궤적이 몇 시간 먼저 끝나므로 먼저 끝난 것들만 모으면 바닥상태가
  과대평가된다. TDC+NAC의 처음 25개는 S0(300 fs)=0.96이었으나 참값은 0.56이다.
- **종료 코드 0이 완주를 뜻하지 않는다.** 입력의 `nstep`과 실제 도달 스텝을 대조할 것.
  PyRAI2MD가 궤적당 300스텝을 정확히 계산해 놓고 한 줄만 기록한 적이 있다. 연장 시
  `direct`를 새 총 스텝 수로 올리고 `buffer 1`로 둘 것.
- **선점 디렉터리는 작업보다 오래 남는다.** 죽은 작업이 남긴 선점이 궤적 하나를 27시간
  묶어 두었고 아무 실패도 보고되지 않았다. `squeue -t all`과 대조할 것.
- **다시 돌리지 말고 재개할 것.** 20스텝마다 체크포인트와 실행 가능한 재시작 매니페스트를
  쓴다. 재개에는 체크포인트, 매니페스트, 조밀 궤적 파일이 함께 필요하다. 재개된 실행은
  새 매니페스트를 쓰지 않고 자기 `resume.oqp`를 유지하므로 연쇄 재개는 두 이름을 모두
  받아야 한다. 재개의 스레드 폭을 바꾸려면 `input(omp_threads=N)`을 다시 써야 하며
  명령줄만으로는 무시된다.
- **맥에서 나오는 NumPy matmul 경고는 잡음이다.** Apple Accelerate가 SIMD 잔여 레인을
  계산하며 IEEE 플래그를 더럽힌다. 정상 실행도 낸다.
- **자원 규칙.** 에너지 효율상 맥미니 우선, chc2는 로그인 노드라 절대 사용 금지,
  chc4는 물리 76코어의 절반까지. 맥미니의 Full NAC은 14코어 전부를 궤적 하나에 쓴다.
- **실행 중인 스크립트를 제자리에서 고치지 말 것.** bash가 파일에서 이어 읽는다.
  새 파일명으로 배포하고 러너를 다시 띄울 것.
- 알려진 흠: `pyoqp/oqp/library/namd.py`의 `NAMD._enforce_nacme_gate`가 두 번 정의되어
  있다(2325행, 2426행). 두 번째만 유효하며 동작은 같다. 병합 이전부터 있던 것이다.

## 1. 2026-08-31 시점의 기록 (0절로 대체됨, 경위 참고용)


1. 현재 Uracil의 TD 및 HT-NAC 500 fs 결과는 원 논문의 KNU-GAMESS 조건을 재현한 결과가 아니다. 두 계산 모두 시간 미분 결합에 `tdc=npi`를 썼고 `decoherence=edc`도 추가했다. 따라서 현재 Figure 8의 조기 S0 전이는 analytic NAC의 정확도 결론으로 사용할 수 없다.
2. IC49의 56.5 fs 지점에서 원인이 수치적으로 확인되었다. raw TLF overlap에 NPI의 polar-orthogonalization과 matrix logarithm을 적용하면 `T10=-0.0354857 au`가 되지만, KNU 방식의 Hammes-Schiffer–Tully finite difference는 `T10=-0.00108107 au`, centered analytic 값은 `-0.0010681 au`이다. NPI가 약 33배 큰 거짓 S1→S0 후보를 만들었다.
3. HT-NAC은 hop 후보가 생긴 뒤에만 analytic derivative-coupling vector를 계산하여 속도 재조정 방향에 사용한다. 따라서 잘못된 NPI hop 후보 자체를 막지 못한다.
4. 다음 계산은 코드를 먼저 고치는 일이 아니라 입력 조건을 KNU 기준으로 맞추는 일이다. `tdc=fd`, `tlf=2`, `decoherence=off`를 사용한 짧은 paired trajectory 1–2개를 먼저 확인한다.
5. 현재 Slurm에서는 기존 full analytic 및 HT-NAC 계산이 아직 실행 중이다. 새 계산을 제출하기 전에 반드시 `squeue`와 결과 경로를 대조하여 중복을 막는다.
6. PR은 열지 않는다. 코드 수정이 필요하면 별도 worktree/branch에서 최소 수정하고 회귀 검사를 통과시킨 뒤 교수님 승인을 받는다.

## 2. 학생의 첫 30분

### 2.1 규칙과 중앙 상태를 읽기 전용으로 확인

다음 global skill을 **처음부터 끝까지** 읽는다.

- `research-operations-chief`
- `openqp-build`
- `chc-compute`
- `mac-mini-compute`
- `openqp-db`
- `synology-operations`
- `mrsf-tddft-theory` — MRSF 이론이나 resident NAC 수식을 수정할 때 필수
- `chc-scientific-writing-style`, `scientific-language`, `scientific-manuscript-editor` — 논문 수정 시 필수

Ultra에서 중앙 기록을 읽기 전용으로 확인한다. 토큰이나 capability를 문서·로그·메시지에 복사하지 않는다.

```bash
ROC=$(echo /Users/cheolhochoi/.research-ops/runtime/dispatchers/*/research-ops-runtime-dispatch)
"$ROC" status
"$ROC" list tasks
"$ROC" resource list
"$ROC" doctor
```

계산 제출·취소·재시작·공유 경로 쓰기·OpenQP.DB 저장은 별도 권한과 중앙 lease가 필요하다. 학생이 스스로 기존 lease를 인계받았다고 가정하면 안 된다. 만료된 lease도 자동으로 비어 있는 자원이 아니다.

### 2.2 현재 Slurm 상태를 먼저 확인

```bash
ssh chc4
sinfo -N -h -o '%P|%N|%t|%C|%c|%X|%Y|%Z|%m|%f'
squeue -u "$USER" -o '%.18i|%.12P|%.28j|%.8T|%.10M|%.10l|%R'
```

2026-08-31 12:36 KST의 실행 상태는 아래와 같다. 이 표는 시간이 지나면 달라지므로 위 명령이 항상 우선이다.

| 방법 | 현재 결과 | 실행 중 Slurm ID |
|---|---:|---|
| TD/NPI isotropic | 50개 파일, 49개 500 fs 완료, IC101 412 fs | 없음 |
| HT-NAC/NPI | 43개 파일, 43개 500 fs 완료 | `5430434_[1-4]` (`r640`) |
| full analytic | 31개 고유 IC 파일, 0개 500 fs 완료, 최대 224 fs | `5430624_[0-13]`, `5430625_14`, `5430666_[0-11]`, `5430667_[1-7]`, `5430701_0`, `5430702_0` |
| RLZT10 | 0개 | 없음 |
| Baeck–An | 6개 500 fs 완료: IC69, 70, 93, 99, 105, 107 | 당시 `5430603`이 IC46/47 실행 중이었으나 8 h 제한 직전이므로 즉시 재확인 |

full analytic의 Slurm 행 수와 고유 결과 디렉터리 수가 다른 것은 여러 retry와 packed worker가 섞여 있기 때문이다. **Slurm job 수를 trajectory 수로 해석하지 말고**, 각 `run-manifest.txt`, trajectory 마지막 시각, 부모 job ID를 대조한다.

Mac mini의 OpenQP 생산 계산은 모두 종료된 상태로 전달되었다. `chc2`는 사용하지 않는다. 로컬에는 10분 간격 감시 프로그램이 남아 있었다.

```text
PID 27371  /tmp/oqp-uracil-ec04414-artifacts/monitor_mac_fleet.sh 600
PID 30797  /tmp/oqp-uracil-ec04414-artifacts/monitor_population_checkpoints.py 600
PID 33359  /tmp/oqp-uracil-ec04414-artifacts/monitor_population_checkpoints.py 600
```

PID는 달라질 수 있다. 학생은 `ps`로 실제 실행 여부를 확인하고, 기존 감시가 살아 있으면 새 감시를 중복 실행하지 않는다.

## 3. 코드와 계산 위치

### 3.1 주 계산 코드

초기 지정 private commit은 다음이다.

```text
repository: myprivate/codex/nac-performance
commit: 8eb9a8d9730a3b821861083268a08e71bd327e9f
```

현재 Uracil 생산 계산에 사용한 코드는 그 commit의 후손이다.

```text
commit: ec04414dc2226c67cd1b59355060faa2d5c0cd4b
remote branch recorded in checkout: origin/agent/ht-nac-20260830
checkout: /bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830/source/openqp-ec04414dc2226c67cd1b59355060faa2d5c0cd4b
build: /bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830/build/attempt-001/venv
later executable used by full retry manifest: /bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830/build/attempt-002/venv/bin/openqp
attempt-002 executable SHA256: 15602f30028b97ff43d2653debde774a309caa9705d537c73ce42c25de615cf3
```

`ec04414` checkout은 detached HEAD이며 깨끗하다. **여기서 직접 수정하지 않는다.** 학생 전용 Git worktree와 branch를 만든다.

주요 코드:

```text
pyoqp/oqp/library/namd.py
  _state_overlap              overlap 및 TDC 계산
  _update_analytic_nac        resident analytic d와 v·d 계산
  _compute_tdc                fd 및 npi 정의
  _hop_triggered_analytic_rescale  HT-NAC 속도 재조정
pyoqp/oqp/library/nac_utils.py
  canonical_state_overlap
  hst_derivative_coupling     (S-S^T)/(2 dt)
pyoqp/oqp/library/nac_analytic.py
  analytic_nac
pyoqp/oqp/library/single_point.py
  NACME.nacme                 TLF overlap 및 FD coupling
```

`ec04414`까지 포함된 주요 변경은 analytic TDC, NAC 방향 속도 재조정, temporal Z-vector predictor, gradient/NAC Z-vector 결합, HT-NAC deferred candidate 계약, XC grid gradient 분리이다. 코드를 새로 구현하기 전에 이 계보를 그대로 사용한다.

### 3.2 Baeck–An 별도 코드

```text
commit: 27dbada1062f3c374faac3bf6d1f40cc63adea02
checkout: /bighome/cheolho.choi/openqp-uracil-baeck-an-27dbada-20260830/source/openqp-27dbada1062f3c374faac3bf6d1f40cc63adea02
campaign: /bighome/cheolho.choi/openqp-uracil-baeck-an-27dbada-20260830
```

이 역시 detached HEAD이며 깨끗하다. 이 구현은 한 step 지연된 Baeck–An coupling magnitude에 phase-tracked overlap sign을 붙이고, 초기 두 step은 NPI로 시작한다. `tdc=baeck_an,rescale=isotropic,nacme_check=off`이다. 360개 회귀 중 360 passed, OpenMM 전용 1개 skip이 기록되어 있다.

학생용 private GitHub branch `codex/student-handoff-20260831`에는 이 exact
`27dbada1062f3c374faac3bf6d1f40cc63adea02` 계보가 merge되어 있다. 따라서 별도
bundle을 다시 가져오지 말고 이 branch에서 학생 전용 worktree를 만든다.

Baeck–An은 full derivative-coupling vector가 아니라 에너지 곡률에 근거한 time-derivative coupling의 크기 근사이다. analytic NAC 정확도와 동일한 물리량이라고 쓰면 안 된다.

### 3.3 계산 campaign

```text
root: /bighome/cheolho.choi/openqp-uracil-namd-ec04414-20260830
input generator: scripts/generate_uracil_inputs.py
selection: scripts/selection.json
build manifest: scripts/build_manifest.json
production inputs: inputs/
production results: results/
Slurm/retry manifests: manifest/
current population analysis: analysis/existing-500fs-td-ht-preliminary-v3/
```

초기 조건은 seed `20260830`으로 고른 50개이다. S1 시작 6개는 14, 40, 42, 61, 74, 107이고, 나머지 44개는 S2 시작이다. `scripts/selection.json`을 기준으로 삼고 목록을 다시 뽑지 않는다.

현재 generator의 네 설정은 다음과 같다.

| 내부 이름 | 현재 설정 | 해석 |
|---|---|---|
| `baseline_npi_iso` | `tdc=npi`, `rescale=isotropic`, EDC | 현재 비교용이지만 KNU 기준 아님 |
| `ht_nac` | `tdc=npi`, `rescale=hop_analytic_nac`, EDC | 잘못된 NPI 후보 뒤 analytic 방향만 사용 |
| `full_analytic` | `tdc=analytic`, `rescale=analytic_nac`, EDC | full analytic 생산 계산 |
| `rlzt10` | analytic + linear Z predictor, exact every 10, EDC | 아직 생산 결과 없음 |

공통 입력은 `dt=0.5 fs`, `substep=50000`, `edc_c=0.1`이다. 논문 원 계산과의 비교에서는 EDC를 끈다.

### 3.4 빌드 정보

CHC Linux build manifest:

```text
GCC/GFortran 12.3.0
Python 3.11.3
Intel MKL 2023.1.0 ILP64 sequential
BLAS/LAPACK integer: 8 byte
USE_LIBINT=OFF
OpenMP=ON
source bundle SHA256: 31b93b55acfe0d2411f8574db0b966087d0af11af08f80900047cfdd91b48b83
```

OpenQP는 bare CMake로 시작하지 않는다. package build와 persistent external cache를 사용한다. Mac 빌드가 새로 필요할 때는 아래 값을 바꾸지 않는다.

```bash
CMAKE_ARGS='-DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 -DUSE_LIBINT=OFF -DENABLE_OPENMP=ON -DLINALG_LIB=auto -DENABLE_OPENTRAH=OFF -DOQP_REUSE_EXTERNALS=ON' \
  /opt/homebrew/bin/python3.11 -m pip install . --no-deps --force-reinstall
```

Libint가 빌드되기 시작하거나, 동일 ABI cache가 있는데 Libxc를 다시 빌드하면 중지하고 원인을 확인한다. 변경 없는 commit·실행 파일·입력·node class·thread 배치에 대해 검사를 반복하지 않는다. 전자구조 계산을 1 core로 시험하지 않는다.

## 4. 원 논문과 Woojin 자료

### 4.1 원자료

Synology에서는 읽기 전용으로 조사한다.

```text
Wigner source:
/volume1/synohome/jin/NX-2-B19/Uracil/final_output_300K

Woojin production:
/volume1/synohome/jin/Projects/0_Done/Uracil/NAMD_300K_10kcal_new

population tables:
/volume1/synohome/jin/Projects/0_Done/Uracil/NAMD_300K_10kcal_new/population_analysis

IC49 input/log:
/volume1/synohome/jin/Projects/0_Done/Uracil/NAMD_300K_10kcal_new/traj/md.49.new.uracil.10kcal.inp
/volume1/synohome/jin/Projects/0_Done/Uracil/NAMD_300K_10kcal_new/traj/md.49.new.uracil.10kcal.log

archived KNU-GAMESS source:
/volume1/synohome/jin/Projects/0_Done/DHA/GAMESS/knu-gamess/source/namd.src
SHA256: 9fbfa1479ffabd1ca1d37f8d0b30fbcd7b4e60c4f7cc1873d37b2806a633369c
archived branch: Konstantin/MRSF-NAMD-REKS
archived HEAD: 8092424faedf4871c28be9cf9d8ca67414266907
```

2021 실행 파일의 정확한 source commit은 보존된 Git checkout만으로 단정할 수 없다. 원 로그는 `GAMESS VERSION = 14 FEB 2018 (R1)`, host `i87k8`, 6 processors를 기록한다.

첨부 SI:

```text
/Users/cheolhochoi/Library/CloudStorage/Dropbox/Exchange/jz2c01694_si_001-2.pdf
```

Figure S6a(page S-20)는 MRSF/BH&HLYP/6-31G*, 88개 성공 trajectory(초기 S2 78개, S1 10개), 핵 dt 0.5 fs, 전자 substep 10^-5 fs, NVE이다. decoherence correction은 기재되어 있지 않다. Figure S7은 torsional coordinate 비교의 기준이다.

### 4.2 KNU time-derivative coupling

KNU-GAMESS는 TLF(2) overlap과 Hammes-Schiffer–Tully 식을 사용한다.

\[
T_{ij}(t+\Delta t/2)=\frac{S_{ij}(t,t+\Delta t)-S_{ji}(t,t+\Delta t)}{2\Delta t}.
\]

KNU `MRSFOV`는 diagonal state overlap을 건너뛰고, `NACVFD`가 위 식을 적용한다. retained-state overlap을 polar-orthogonalize하지 않으며 matrix logarithm도 사용하지 않는다.

따라서 원 논문을 재현할 기준 입력은 다음이다.

```text
TD reference: tdc=fd, tlf=2, decoherence=off, rescale=isotropic
HT-NAC:       tdc=fd, tlf=2, decoherence=off, rescale=hop_analytic_nac
Full analytic:tdc=analytic, decoherence=off, rescale=analytic_nac
Baeck-An:     별도 27dbada 구현; decoherence=off 조건을 별도 확인
```

네 방법의 비교에서는 공통 SCF, Davidson, gradient, 핵 적분, RNG, 초기 조건을 같게 한다. TD와 HT의 차이는 hop 이후 속도 재조정 방향뿐이어야 한다. Full analytic은 coefficient propagation과 hop probability에도 \(\dot{\mathbf R}\cdot\mathbf d_{ij}\)를 사용한다.

## 5. 확인된 NPI 결함

IC49, current HT, 56.5 fs의 raw retained-state overlap은 다음이다.

```text
 0.026248881  -0.010838567  -0.015653126   0.008249603
-0.055531473  -0.146177222   0.158574974   0.073965082
-0.009270643   0.110538644  -0.077338506  -0.030524060
 0.073092593  -0.069540967  -0.083798296  -0.031516829
```

singular values는 `0.270910906, 0.140010333, 0.020542510, 0.013489465`이다. retained-state 공간 밖으로 큰 population이 빠졌기 때문에 이 행렬은 unitary step propagator에 가깝지 않다. 그런데 현재 NPI는

\[
Q=S(S^T S)^{-1/2},\qquad T=\log(Q)/\Delta t
\]

를 적용한다. 이 지점에서 `det(Q)=-1`이고 \(\pi\) branch가 생겨 real antisymmetric logarithm 가정이 성립하지 않는다.

수치 비교:

```text
NPI T10                    -0.0354857 au
KNU HST/FD T10             -0.00108107 au
centered analytic T10      -0.0010681 au
endpoint analytic T10      -0.0022023 au
```

현재 HT trajectory의 상태 이력은 S2→S1 at 17.5 fs, S1→S0 candidate at 56.5 fs이다. 56.0→56.5 fs population은 `(0.000836,0.970351,0.027178,0.001635)`에서 `(0.096583,0.590085,0.256817,0.056515)`로 급변했고 RNG는 `0.017942`였다.

Woojin IC49는 같은 시각에 S2에 남지만, 이 한 경로는 이미 14.5–17.0 fs 사이 hop 이력이 달라졌으므로 pointwise ensemble 증거로 쓰면 안 된다. 원 경로와 현재 경로의 RNG stream도 다르다.

### 현재 population의 정량적 차이

현재 paired 26-trajectory 분석:

| time/fs | TD/NPI S0 | HT/NPI S0 |
|---:|---:|---:|
| 40 | 0.1154 | 0.1538 |
| 50 | 0.1923 | 0.1923 |
| 60 | 0.2692 | 0.3462 |
| 75 | 0.4231 | 0.6154 |
| 100 | 0.5000 | 0.6538 |

Figure S6a에서 읽은 대략적인 S0 population은 40 fs에서 거의 0, 50 fs에서 0.03–0.05, 60 fs에서 0.10–0.15, 75 fs에서 0.25–0.30, 100 fs에서 약 0.45이다. 표본 수 26 대 88의 차이는 불확실성에 포함해야 하지만, 초기 차이는 NPI와 EDC 조건 불일치의 영향을 먼저 제거해야 한다.

## 6. 다음 계산의 정확한 순서

### 단계 A — 기존 작업 보존과 재조정

1. `squeue`와 결과 디렉터리를 다시 읽는다.
2. 각 실행 중 job의 `run-manifest.txt`, 마지막 trajectory 시각, log 증가, CPU affinity를 확인한다.
3. 성공·실패·timeout·cancelled를 분리한다. 이름에 `.failed-<jobid>`가 있어도 dense trajectory가 500 fs에 도달한 경우가 있으므로 파일 내용을 검사한다.
4. 기존 결과를 덮어쓰지 않는다. retry는 새 attempt 디렉터리와 부모 job ID를 사용한다.
5. 현재 NPI/EDC 결과는 `diagnostic-only`로 보존한다. 삭제하지 않는다.

### 단계 B — 1–2개 paired trajectory 확인

새 worktree/branch와 새 calculation root를 만든 뒤 IC49와 정상 overlap을 보이는 IC 하나만 사용한다. 252개 진단 matrix나 대규모 trajectory를 만들지 않는다.

비교할 두 입력:

```text
TD-FD: tdc=fd, tlf=2, decoherence=off, rescale=isotropic
HT-FD: tdc=fd, tlf=2, decoherence=off, rescale=hop_analytic_nac
```

동일 geometry, velocity, initial state, seed, `rng_stream`, `dt=0.5 fs`, `substep=50000`을 사용한다. 처음부터 7–8 OpenMP threads로 실행한다. 이미 통과한 동일 build/node-class 확인은 반복하지 않는다.

확인할 값:

- raw overlap \(S\), singular values, determinant
- state assignment와 phase sign
- FD \(T_{ij}\), centered analytic \(\dot{R}\cdot d_{ij}\)
- coefficient와 population의 합
- 각 hop probability와 random number
- candidate hop, energy 허용/거절, velocity-rescaling 허용/거절
- S0↔S1 양방향 후보와 back transition
- total-energy drift

IC49의 56.5 fs에서 FD T10이 약 `-0.00108 au`이고 NPI의 `-0.0355 au` spike가 사라지는지가 첫 판정 기준이다. HT-FD가 candidate를 만들지 않으면 analytic rescaling은 호출되지 않는 것이 정상이다.

### 단계 C — 네 방법의 50-trajectory 생산 계산

짧은 paired 확인이 통과한 뒤에만 같은 50개 IC로 실행한다.

1. `TD`: TLF(2) HST finite-difference TDC + isotropic rescaling
2. `Baeck–An`: energy-curvature TDC approximation + isotropic rescaling
3. `HT-NAC`: TD-FD hop probability + hop 시에만 exact analytic NAC 방향 rescaling
4. `Full analytic`: 매 step exact analytic \(\dot{R}\cdot d_{ij}\) + analytic NAC 방향 rescaling

RLZT는 네 방법의 주 비교에 바로 넣지 않는다. 먼저 full analytic exact를 기준으로 Z-vector prediction의 \(d\), \(h\), \(v\cdot d\), hop decision 오차를 정량화한 뒤 SI의 가속화 방법으로 다룬다.

자원 배치:

- Slurm의 모든 호환 node class를 실시간 확인한다. `r630`, `r640`, `ryzn`, `trd`, `xeon`을 한 종류로 가정하지 않는다.
- 독립 trajectory는 7–8 OpenMP threads/job을 기본으로 하고 물리 core를 빈틈없이 채운다. BLAS threads는 1이다.
- `r630` 28 physical core: 4 × 7-thread worker.
- `r640` 40 physical core: 5 × 8-thread worker. r630 packed launcher를 그대로 쓰지 않는다.
- `ryzn16–40` 일부는 affinity상 8 physical core만 노출했다. `ryzn103/108` 등 16-core 노드와 분리한다.
- Mac M4 mini 14 physical core: host당 2 × 7-thread worker. 단, method 수정 전에는 재시작하지 않는다.
- `chc2`는 사용하지 않는다.
- 예측 시간이 짧은 TD, HT, Baeck–An부터 제출하고 full analytic을 마지막에 채운다. 다만 이미 실행 중인 full analytic을 순서 변경만을 위해 취소하지 않는다.

### 단계 D — 분석

population은 모든 electronic state를 한 그림에 표시한다. 50 fs마다 임시 population figure를 만들되, 서로 다른 방법은 동일한 완료 trajectory 집합과 동일한 시간 범위로 비교한다. 최종 500 fs 그림에는 500 fs를 완주한 동일 IC만 쓴다.

추가 분석:

- Figure S6a와 같은 adiabatic population
- Figure S7과 같은 ring-puckering/torsional coordinates
- 40 fs 이전의 빠른 decay channel 비율
- first-hop time 및 S0 도달 시간 분포
- forward/back hop 수, frustrated hop 수
- analytic 방향과 isotropic 방향에 따른 product channel 비율
- 각 방법의 wall time/step, electronic calls, SCF/Davidson/gradient/NAC 시간

성능 비교에서는 공통 SCF, Davidson, active-state gradient 시간을 기준으로 분리하고, resident analytic NAC 자체의 추가 시간과 step 전체 시간을 둘 다 제시한다. time-derivative TLF는 scalar \(T_{ij}\), analytic NAC은 full \(3N\) vector라는 차이를 반드시 쓴다.

## 7. 논문

Overleaf:

```text
remote: https://git@git.overleaf.com/6a6e9cf68a2d15fed156be72
current remote main at handoff: c01989e438d973c713775e54dfe5dbad7a3d753b
commit message: Replace uracil TD population panel with 500 fs ensemble
```

학생은 최신 `origin/main`에서 자기 이름이 포함된 새 worktree/branch를 만든다. 다른 사람의 worktree를 공유하지 않는다. push 전에는 remote SHA를 다시 읽고 fast-forward만 허용한다. plain force push는 금지한다.

코드와 이 인수인계 문서를 함께 보존한 private GitHub branch:

```text
repository: https://github.com/karmachoi/openqp-private
branch: codex/student-handoff-20260831
analytic/HT source: ec04414dc2226c67cd1b59355060faa2d5c0cd4b
Baeck-An source: 27dbada1062f3c374faac3bf6d1f40cc63adea02
```

현재 논문의 3.4–3.6은 네 NAMD 방법, population, torsional motion으로 재구성하는 방향이다. 그러나 현재 NPI/EDC population figure는 조건이 잘못 맞춰졌으므로 최종 수치로 사용하지 않는다. 새 TD-FD/HT-FD 결과로 교체하고, 표본 수와 confidence interval을 함께 쓴다.

논문의 중심 주장:

1. analytic NAC은 coordinate finite-difference scan과 step/phase 의존성을 제거한다.
2. CI에서 \(h=(E_1-E_0)d\)와 2D branching plane을 직접 제공한다.
3. full \(3N\) vector를 numerical spatial TLF보다 낮은 비용으로 제공한다.
4. MRSF branching plane은 외부 MRCISD 기준과 정합한다.
5. Uracil NAMD에서는 analytic coupling의 정확도뿐 아니라 hop probability와 velocity-rescaling 방향이 population 및 product channel에 미치는 영향을 보인다.

`analytic–TLF comparison`을 정확도 주장의 중심으로 쓰지 않는다. 이것은 구현의 수치적 확인이다. H2에서 FCI와 직접 비교하는 것은 정확도 평가로 의미가 있다. MRCISD comparison은 MRSF의 CI/branching-plane 정확도를 평가한다.

기존 그림 및 생성 script는 Overleaf 저장소의 `figures/`와 `scripts/`에 있다. 주요 그림은 H2 dissociation NAC, branching-plane vector agreement, MECI double cone, Berry phase, MECI 최적화, Uracil population/torsion, 성능 및 Z-vector predictor이다. figure를 바꾸면 source data와 manifest도 함께 보존한다.

## 8. 아직 해결되지 않은 일

1. TD-FD/decoherence-off 짧은 paired 확인이 아직 없다.
2. HT-NAC의 조기 S0 전이가 FD에서 사라지는지 아직 직접 확인하지 않았다.
3. Baeck–An은 50개 중 6개만 500 fs 완료된 상태로 확인되었고, 현 생산 방식은 시간이 길다.
4. full analytic은 500 fs 완주 trajectory가 아직 없다. retry lineage와 고유 IC를 먼저 정리해야 한다.
5. RLZT10 생산 결과가 없다. Z-vector predictor는 exact analytic 기준 오차를 먼저 제시해야 한다.
6. Woojin 88-trajectory 조건과 현재 50-trajectory 선택의 표본 차이를 통계적으로 다뤄야 한다.
7. 원 2021 KNU-GAMESS 실행 파일의 정확한 Git commit은 확정되지 않았다.
8. final product branching ratio와 40 fs 이내 특수 decay channel 분석이 남아 있다.
9. 현재 500 fs NPI/EDC figure와 그 해석은 수정되어야 한다.
10. 완료 계산은 OpenQP.DB에 입력·출력·manifest·실행 파일 hash와 함께 저장해야 한다. 실패와 불완전 결과는 성공으로 표시하지 않는다.

## 9. 하지 말아야 할 일

- 현재 NPI/EDC 결과를 원 Figure S6 재현 또는 analytic NAC 정확도 증거로 발표하지 않는다.
- 불량 overlap을 unitary step propagator로 강제한 NPI 결과를 그대로 사용하지 않는다.
- 1-core full analytic 확인으로 시간을 낭비하지 않는다.
- 동일 build와 동일 node class에 대해 import/build 검사를 반복하지 않는다.
- `.inp` 대신 현재 campaign의 `.oqp` 형식을 사용한다.
- 50개 trajectory가 같은 파일명이나 scratch 디렉터리를 쓰게 하지 않는다.
- 실패 디렉터리를 덮어쓰거나 삭제하지 않는다.
- source checkout 또는 Synology 원자료를 직접 수정하지 않는다.
- `chc2`를 사용하지 않는다.
- 허가 없이 job 취소·재제출, Overleaf/GitHub push, PR, registry 수정, OpenQP.DB 쓰기를 하지 않는다.
- 교수님의 명시적 승인 없이 upstream fast-forward 또는 branch 변경을 하지 않는다.

## 10. 학생이 교수님께 첫날 보고할 형식

```text
1. 읽은 source commit/build hash:
2. live Slurm job ID와 고유 trajectory 수:
3. TD/HT/full/BA 각각 500 fs 완료·부분·실패 수:
4. IC49 TD-FD와 HT-FD 확인 결과:
   raw overlap singular values, FD T10, analytic T10, hop candidate/accept/reject
5. 원 Figure S6 대비 40/50/60/75/100 fs S0 population:
6. 새로 제출할 정확한 calculation identities와 사용 node/thread 배치:
7. 남은 과학적 문제와 필요한 승인:
```

핵심은 계산 수를 늘리는 것이 아니라, 먼저 KNU와 동일한 TLF(2) finite-difference TDC 및 NVE 조건을 맞춘 뒤 analytic NAC가 **hop probability**와 **속도 재조정 방향**에 미치는 영향을 분리하는 것이다.
