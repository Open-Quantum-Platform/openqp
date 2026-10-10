# fix/formaldehyde-gradient-highroot

SHA: `2a0caadc37c9da6b5a62952bd57dbe1542b3b16a`  
판정: **고유 변경·반영 여부 검토**  
분야: SCF / response (검색용 분류)  
마지막 commit: 2026-05-27T14:27:49+09:00 / cheolhochoi  
제목: test: document MRSF XC gradient density handoff gap

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / fix/formaldehyde-gradient-highroot](https://github.com/karmachoi/openqp/tree/2a0caadc37c9da6b5a62952bd57dbe1542b3b16a)
- [github-private / fix/formaldehyde-gradient-highroot](https://github.com/karmachoi/openqp-private/tree/2a0caadc37c9da6b5a62952bd57dbe1542b3b16a)
- [gitlab / fix/formaldehyde-gradient-highroot](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/2a0caadc37c9da6b5a62952bd57dbe1542b3b16a)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `987ae53522ae4e11c3c7dcd904b04faa60f9e2a2`
- 앞선 커밋 8 / 뒤처진 커밋 749
- non-merge patch: main과 일치 0, 다름 8
- merge/empty 등 patch 비교 제외: 0
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 6; 그중 현재 main과 동일 0, 다름 6

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.


## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `2a0caadc37c9` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | test: document MRSF XC gradient density handoff gap |
| `854a7b952916` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | Revert "fix: pass MRSF transition density to XC gradient" |
| `ca26976f1f59` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | fix: pass MRSF transition density to XC gradient |
| `67243ca2e693` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | Revert "fix: preserve beta MO channel in MRSF Z-vector" |
| `217a3dc75e5c` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | fix: preserve beta MO channel in MRSF Z-vector |
| `cf57549d7425` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | test: document MRSF SPC operator mapping |
| `ec9697f59464` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | fix: preserve MRSF channel-7 ball density |
| `464c084116f6` | 2026-05-27 | 다름 | commit 보존 또는 비교 제외 | test: add gradient state-character diagnostics |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `source/modules/tdhf_mrsf_z_vector.F90` | 다름 |
| A | `tests/test_gradient_state_character_diagnostics.py` | 다름 |
| A | `tests/test_mrsf_channel7_ball_consistency.py` | 다름 |
| A | `tests/test_mrsf_spc_operator_consistency.py` | 다름 |
| A | `tests/test_mrsf_xc_gradient_density.py` | 다름 |
| A | `tools/gradient_state_character.py` | 다름 |
