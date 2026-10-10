# codex/reuse-external-build-cache

SHA: `3f7da4ad90af0cdd894dabd8feaf37c5116d14d5`  
판정: **동일 끝점 PR 병합 확인**  
분야: Build / CI / release (검색용 분류)  
마지막 commit: 2026-06-12T18:05:23+09:00 / Cheol Ho Choi  
제목: build: hash resolved BLAS provider into externals cache key

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / codex/reuse-external-build-cache](https://github.com/karmachoi/openqp/tree/3f7da4ad90af0cdd894dabd8feaf37c5116d14d5)
- [gitlab / codex/reuse-external-build-cache](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/3f7da4ad90af0cdd894dabd8feaf37c5116d14d5)
- `local-11 / pr203-cache-key`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `6687782f22b3a55c2efba8eebef5f968f03bf614`
- 앞선 커밋 2 / 뒤처진 커밋 713
- non-merge patch: main과 일치 0, 다름 2
- merge/empty 등 patch 비교 제외: 0
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 2; 그중 현재 main과 동일 0, 다름 2

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #204: ci: key the externals cache by full compiler version](https://github.com/Open-Quantum-Platform/openqp/pull/204) — state=closed; merged=2026-06-12T09:34:34Z; merge commit in main=True

## 같은 이름의 PR — 끝점이 다르므로 전체 반영 증거가 아님

- [upstream #203: [codex] Cache bundled external builds](https://github.com/Open-Quantum-Platform/openqp/pull/203) — closed; PR head `d29e7111234b`; merged=2026-06-12T08:44:41Z

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `3f7da4ad90af` | 2026-06-12 | 다름 | commit 보존 또는 비교 제외 | build: hash resolved BLAS provider into externals cache key |
| `cc8b878c14c7` | 2026-06-12 | 다름 | commit 보존 또는 비교 제외 | ci: key the externals cache by full compiler version |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.github/workflows/CI.yml` | 다름 |
| M | `external/CMakeLists.txt` | 다름 |

## 이 끝점을 사용 중인 로컬 worktree

- `/Users/cheolhochoi/Documents/claude/openqp-pr203` — exists=True, status 항목 0; 
