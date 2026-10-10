# cl

SHA: `1686aecbaca467f2acbd09e61123d2a3bfc13b94`  
판정: **고유 변경·반영 여부 검토**  
분야: Build / CI / release (검색용 분류)  
마지막 commit: 2026-05-31T20:43:25+09:00 / cheolhochoi  
제목: ci: skip Claude auto review for fork PRs

[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)

## 같은 끝점을 가리키는 모든 원본

- [github-personal / cl](https://github.com/karmachoi/openqp/tree/1686aecbaca467f2acbd09e61123d2a3bfc13b94)
- [github-private / cl](https://github.com/karmachoi/openqp-private/tree/1686aecbaca467f2acbd09e61123d2a3bfc13b94)
- [gitlab / cl](https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/1686aecbaca467f2acbd09e61123d2a3bfc13b94)

## 현재 GitLab main과의 비교

- 기준 main: `d8fcc4119c73a2be3c1c7e5f49a934f4d7ccad92`
- 공통 조상: `795deaf25fc61c3200376b7d07b064db3847a94a`
- 앞선 커밋 3 / 뒤처진 커밋 744
- non-merge patch: main과 일치 0, 다름 3
- merge/empty 등 patch 비교 제외: 0
- GitLab 어느 브랜치에서도 끝점 도달 가능: True
- GitLab 전체에서 동일 patch를 못 찾은 커밋: 0
- 공통 조상 이후 변경 파일: 26; 그중 현재 main과 동일 6, 다름 20

앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.

## 이 끝점과 정확히 일치하는 PR/MR

- [upstream PR #164: Remove defined optimizer and make ILP64 BLAS/LAPACK mandatory](https://github.com/Open-Quantum-Platform/openqp/pull/164) — state=closed; merged=아니오; merge commit in main=False

## 추가 커밋 전체

| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |
| --- | --- | --- | --- | --- |
| `1686aecbaca4` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | ci: skip Claude auto review for fork PRs |
| `5d0f1f4a05a4` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | chore: remove static ILP64 config test |
| `b2b55bc7ad02` | 2026-05-31 | 다름 | commit 보존 또는 비교 제외 | Remove defined optimizer and require ILP64 BLAS |

## 변경 파일 전체

| 상태 | 경로 | 현재 main과 내용 |
| --- | --- | --- |
| M | `.github/workflows/CI.yml` | 다름 |
| M | `.github/workflows/claude.yml` | 다름 |
| M | `.gitlab-ci.yml` | 다름 |
| M | `CMakeLists.txt` | 다름 |
| M | `Dockerfile` | 다름 |
| M | `README.md` | 다름 |
| M | `cmake/oqp_functions.cmake` | 다름 |
| M | `cmake/patches/OpenTrustRegion.CMakeLists.txt` | 동일 |
| M | `external/CMakeLists.txt` | 다름 |
| M | `pyoqp/CMakeLists.txt` | 다름 |
| M | `pyoqp/README.md` | 다름 |
| D | `pyoqp/oqp/library/libdlfind.py` | 동일 |
| M | `pyoqp/oqp/library/libscipy.py` | 다름 |
| M | `pyoqp/oqp/library/runfunc.py` | 다름 |
| M | `pyoqp/oqp/molecule/oqpdata.py` | 다름 |
| M | `pyoqp/oqp/utils/file_utils.py` | 다름 |
| M | `pyoqp/oqp/utils/input_checker.py` | 다름 |
| D | `pyoqp/patch.sh` | 동일 |
| M | `pyoqp/requirements.txt` | 동일 |
| M | `pyoqp/setup.py` | 다름 |
| M | `pyproject.toml` | 다름 |
| M | `source/CMakeLists.txt` | 다름 |
| M | `source/dftlib/CMakeLists.txt` | 동일 |
| M | `source/integrals/CMakeLists.txt` | 동일 |
| M | `source/mathlib/CMakeLists.txt` | 다름 |
| M | `tests/test_geometric_optimizer.py` | 다름 |
