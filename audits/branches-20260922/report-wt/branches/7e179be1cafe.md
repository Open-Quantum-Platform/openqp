# perf/int2-pair-decode-opt

SHA: `7e179be1cafee9b22ebf949190cd579f5cb8007b`  
판정: **동일 끝점 PR 병합 확인**  
분야: ERI / integrals (검색용 분류)  
마지막 commit: 2026-05-27T04:57:37+09:00 / cheolhochoi  
제목: Merge branch 'main' into perf/int2-pair-decode-opt

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-private / perf/int2-pair-decode-opt](https://github.com/karmachoi/openqp-private/tree/7e179be1cafee9b22ebf949190cd579f5cb8007b)
- [gitlab / github-private/perf/int2-pair-decode-opt](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/7e179be1cafee9b22ebf949190cd579f5cb8007b)
- `local-05 / resolve-pr150-conflicts`

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `d33f9cce9af06b5e52f85c310245015c87bb7848`
- 앞선 커밋 5 / 뒤처진 커밋 752
- non-merge patch: main과 일치 2, 다름 0
- merge/empty 등 patch 비교 제외: 3
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 2; 그중 현재 main과 동일 0, 다름 2

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #150: perf: precompute int2 shell-pair decode map](https://github.com/Open-Quantum-Platform/openqp/pull/150) — state=closed; merged=2026-05-26T20:32:15Z; merge commit in main=True

## 같은 이름의 PR — 끝점이 다르므로 전체 반영 증거가 아님

- [upstream #145: perf: precompute int2 shell-pair decode map](https://github.com/Open-Quantum-Platform/openqp/pull/145) — closed; PR head `392e95a54c20`; merged=아니오

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `7e179be1cafe` | 2026-05-27 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into perf/int2-pair-decode-opt |
| `392e95a54c20` | 2026-05-26 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into perf/int2-pair-decode-opt |
| `f4cfad69b733` | 2026-05-26 | merge/empty/비교 제외 | commit 보존 또는 비교 제외 | Merge branch 'main' into perf/int2-pair-decode-opt |
| `c05c18f6ba80` | 2026-05-26 | 일치 | commit 보존 또는 비교 제외 | perf: precompute int2 shell-pair decode map |
| `b5ed10175121` | 2026-05-26 | 일치 | commit 보존 또는 비교 제외 | fix: flatten int2 OpenMP shell-pair workshare |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `source/integrals/int2.F90` | 다름 |
| M | `tests/test_int2_openmp_workshare.py` | 다름 |
