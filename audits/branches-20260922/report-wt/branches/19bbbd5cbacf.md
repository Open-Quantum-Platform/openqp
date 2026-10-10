# feat/dftbplus-external-backend

SHA: `19bbbd5cbacf36eb311ce732ba5016fdda5e8d43`  
판정: **고유 변경·반영 여부 검토**  
분야: DFTB / xTB / QM/MM (검색용 분류)  
마지막 commit: 2026-05-27T22:13:19+09:00 / cheolhochoi  
제목: test: add DFTB+ parameter env fallback

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- `local-05 / feat/dftbplus-external-backend`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `199f526b0b0acdf1ba967409866b20f4e7d3ab81`
- 앞선 커밋 13 / 뒤처진 커밋 756
- non-merge patch: main과 일치 0, 다름 12
- merge/empty 등 patch 비교 제외: 1
- GitLab 어느 브랜치에서도 끝점 도달 가능: False
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 2
- 공통 조상 이후 변경 파일: 4; 그중 현재 main과 동일 0, 다름 4

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #155: Add optional external DFTB+ backend](https://github.com/Open-Quantum-Platform/openqp/pull/155) — state=closed; merged=아니오; merge commit in main=False

## 같은 이름의 PR — 끝점이 다르므로 전체 반영 증거가 아님

- [upstream #143: feat: add external DFTB+ backend](https://github.com/Open-Quantum-Platform/openqp/pull/143) — closed; PR head `0504a84112f1`; merged=2026-05-26T12:08:02Z
- [upstream #141: Add optional DFTB+ external backend](https://github.com/Open-Quantum-Platform/openqp/pull/141) — closed; PR head `e5130256d28d`; merged=2026-05-26T09:35:44Z
- [upstream #139: Add optional external DFTB+ backend](https://github.com/Open-Quantum-Platform/openqp/pull/139) — closed; PR head `007cd5a6c177`; merged=아니오
- [upstream #137: Add optional external DFTB+ backend](https://github.com/Open-Quantum-Platform/openqp/pull/137) — closed; PR head `a989e05cd923`; merged=아니오

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `19bbbd5cbacf` | 2026-05-27 | 다름 | 미발견 | test: add DFTB+ parameter env fallback |
| `c8474bc5381a` | 2026-05-27 | 다름 | 미발견 | test: add DFTB+ availability skip probe |
| `0504a84112f1` | 2026-05-26 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into feat/dftbplus-external-backend |
| `4a5db42c7a5a` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | docs: condense DFTB README entry |
| `e5130256d28d` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | fix: fail fast when DFTB+ gradients are unavailable |
| `b838848bdbfa` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | Restore DFTB+ schema and fixtures after upstream sync |
| `a96139ec856e` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | Add DFTB+ optimization validation reference |
| `11e0046f0d1f` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | fix: transpose DFTB+ column-major forces |
| `5e42351cd6fe` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | Fix DFTB optimizer without SciPy |
| `139cd4e4c2b3` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | Add DFTB+ capability matrix and optimize path |
| `e2d2923171ee` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | fix: validate DFTB+ gradient runner with real output |
| `ac373802540e` | 2026-05-26 | 다름 | commit 보존 또는 비교 제외 | Preserve optional DFTB+ work directories |
| `5e19b00a75fb` | 2026-05-25 | 다름 | commit 보존 또는 비교 제외 | Allow DFTB+ ground-state gradient input |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `README.md` | 다름 |
| M | `pyoqp/oqp/library/dftbplus.py` | 다름 |
| M | `tests/test_dftbplus_backend.py` | 다름 |
| M | `tests/test_geometric_optimizer.py` | 다름 |
