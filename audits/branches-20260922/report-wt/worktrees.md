# 로컬 worktree 상태

tracked/untracked 경로만 확인했으며 파일을 변경하지 않았다. 상태는 조사 시점의 스냅샷이다.

## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-reintegrate-devkit-20260921

- 원본 `local-01`, branch `refs/heads/codex/reintegrate-devkit-20260921`, HEAD `3799e40e2c171b81c4d21ad5b998a98b502a6fbc`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-worktrees/cheol-ci-codex-review-20260921

- 원본 `local-02`, branch `refs/heads/codex/cheol-ci-codex-review-20260921`, HEAD `6d3f42b2279f385d4d6af39b39883924a920bb20`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-worktrees/gitlab-auto-review-20260920

- 원본 `local-03`, branch `refs/heads/codex/gitlab-auto-review-20260920`, HEAD `eb97168c7ad32be04996409256daa0a8806394a7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-worktrees/github-main-inbound-merge-20260920

- 원본 `local-03`, branch `refs/heads/codex/github-main-inbound-merge-20260920`, HEAD `5668ef511d43762a6e2487cb69643e3927fc2550`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-worktrees/github-sync-automerge-retry-20260920

- 원본 `local-03`, branch `refs/heads/codex/github-sync-automerge-retry-20260920`, HEAD `3035051a54dc2bcfdf7d16ac0cab7ad1c8d70bd8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP

- 원본 `local-04`, branch `refs/heads/main`, HEAD `0000000000000000000000000000000000000000`
- 존재: True; status 항목: 4

```text
?? audits/
?? outputs/
?? tmp/
?? work/
```
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-analytic-ht-20260903

- 원본 `local-04`, branch `refs/heads/claude/analytic-ht-transition-dipole-20260903`, HEAD `060af26aa1ba442512b789ca103570e7690aa53c`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-claude-mrsf-hessian-handoff-20260902

- 원본 `local-04`, branch `refs/heads/codex/claude-mrsf-hessian-handoff-20260902`, HEAD `74ff3cd78f735fb3ace5e5fd59991acd463730cf`
- 존재: True; status 항목: 11

```text
 M .github/workflows/CI.yml
 M pyoqp/oqp/runtime.py
 M pyoqp/oqp/utils/mpi_utils.py
D  pyoqp/requirements.txt
 M pyproject.toml
 M source/modules/tdhf_hessian_response.F90
 M source/modules/tdhf_mrsf_hessian_amplitude.F90
 M source/tdhf_mrsf_lib.F90
 M tests/test_pcm_canonical_runtime_path.py
?? pyoqp/oqp/runtime.py.bak
?? source/modules/tdhf_hessian_response.F90.bak
```
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-fock-basis-fix-20260901

- 원본 `local-04`, branch `refs/heads/codex/mrsf-fock-basis-fix-20260901`, HEAD `6f1e569c7a8195054b382d4da29023136b78bd2a`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-formaldehyde-debug-20260902

- 원본 `local-04`, branch `refs/heads/codex/mrsf-formaldehyde-debug-20260902`, HEAD `46c9b92950f9e5eeb0f13df7d5d8eeec2e99ab20`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-formaldehyde-tpx-build-20260902

- 원본 `local-04`, branch `refs/heads/codex/mrsf-formaldehyde-tpx-build-20260902`, HEAD `a5ee259bd1407c4cdb9d41233639d876383cdbf1`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-hessian-20260901

- 원본 `local-04`, branch `refs/heads/codex/mrsf-tddft-hessian-20260901`, HEAD `879b6b60953b4b1f88c573047b47648541e1b031`
- 존재: True; status 항목: 40

```text
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN.freq.molden
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN.hess.json
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN_NOXC.freq.molden
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN_NOXC.hess.json
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN_NOXC.inp
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN_NOXC_scf_rohf_bhhlyp_6-31g.molden
?? examples/HESS/H2O_BHHLYP-MRSF_ANALYTIC_HESSIAN_scf_rohf_bhhlyp_6-31g.molden
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN.freq.molden
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN.hess.json
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN.json
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_H0005.freq.molden
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_H0005.hess.json
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_H0005.inp
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_H0005.json
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_H0005_num_hess/
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_H0005_scf_rohf_bhhlyp_6-31g.molden
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_NOXC.freq.molden
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_NOXC.hess.json
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_NOXC.inp
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_NOXC.json
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_NOXC_num_hess/
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_NOXC_scf_rohf_bhhlyp_6-31g.molden
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_num_hess/
?? examples/HESS/H2O_BHHLYP-MRSF_NUMERICAL_HESSIAN_scf_rohf_bhhlyp_6-31g.molden
?? examples/HESS/H2O_HF-MRSF_ANALYTIC_HESSIAN.freq.molden
?? examples/HESS/H2O_HF-MRSF_ANALYTIC_HESSIAN_scf_rohf_hf_6-31g.molden
?? examples/HESS/H2O_HF-MRSF_NUMERICAL_HESSIAN.freq.molden
?? examples/HESS/H2O_HF-MRSF_NUMERICAL_HESSIAN.hess.json
?? examples/HESS/H2O_HF-MRSF_NUMERICAL_HESSIAN.json
?? examples/HESS/H2O_HF-MRSF_NUMERICAL_HESSIAN_num_hess/
?? examples/HESS/H2O_HF-MRSF_NUMERICAL_HESSIAN_scf_rohf_hf_6-31g.molden
?? examples/HESS/H2O_SVWN-MRSF_ANALYTIC_HESSIAN.freq.molden
?? examples/HESS/H2O_SVWN-MRSF_ANALYTIC_HESSIAN.hess.json
?? examples/HESS/H2O_SVWN-MRSF_ANALYTIC_HESSIAN_scf_rohf_svwn_6-31g.molden
?? examples/HESS/H2O_SVWN-MRSF_NUMERICAL_HESSIAN.freq.molden
?? examples/HESS/H2O_SVWN-MRSF_NUMERICAL_HESSIAN.hess.json
?? examples/HESS/H2O_SVWN-MRSF_NUMERICAL_HESSIAN.json
?? examples/HESS/H2O_SVWN-MRSF_NUMERICAL_HESSIAN_num_hess/
?? examples/HESS/H2O_SVWN-MRSF_NUMERICAL_HESSIAN_scf_rohf_svwn_6-31g.molden
?? examples/HESS/hess.status
```
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-hessian-xc-perf-20260904

- 원본 `local-04`, branch `refs/heads/claude/mrsf-hessian-xc-perf-20260904`, HEAD `fdc60421be088b1a4c77c0d9d6d2ac290350c62d`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-s0-diagnostics-20260901

- 원본 `local-04`, branch `refs/heads/codex/mrsf-s0-diagnostics-20260901`, HEAD `718bcc3329ca27bad43f7953d151ec6d89307dde`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-s0-screening-fix-20260901

- 원본 `local-04`, branch `refs/heads/codex/mrsf-s0-screening-fix-20260901`, HEAD `5fcabbd3fea6bc3bb75ecc18a2a99ec6d87b0eb7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-mrsf-seam-residual-20260904

- 원본 `local-04`, branch `refs/heads/claude/mrsf-hessian-seam-residual-20260904`, HEAD `5f20bba43fb77046b30ff38cf6b6e60d48b48d9c`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-nac-analytic-20260903

- 원본 `local-04`, branch `detached`, HEAD `492887b801be7d4512ba622777ec9a3568f3f009`
- 존재: True; status 항목: 1

```text
 M pyoqp/oqp/utils/input_checker.py
```
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-pr391-closed-shell-hessian-opt-20260902

- 원본 `local-04`, branch `refs/heads/codex/private-pr391-closed-shell-hessian-opt-20260902`, HEAD `c081e7b8404f3128a326391772364ffee97991f6`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-readme-foldable-methods-20260829

- 원본 `local-04`, branch `refs/heads/codex/readme-foldable-methods-20260829`, HEAD `390121ddaa551e711cfecae48848a0b46970adc8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-readme-foldable-tests-20260829

- 원본 `local-04`, branch `refs/heads/codex/readme-foldable-parser-tests-20260829`, HEAD `f95091ae73ecbf6e21949ef98a23215767a5ad90`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-readme-oqp-studio

- 원본 `local-04`, branch `refs/heads/codex/readme-oqp-studio`, HEAD `eec316edae2cbe38a1fada30441a0b3281b848ef`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-readme-oqp-studio-mo-20260829

- 원본 `local-04`, branch `refs/heads/codex/readme-oqp-studio-mo-20260829`, HEAD `6ba58184f2c51a1ac4185bcf48d2684a6c7bb083`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-rks-hessian-speed-20260904

- 원본 `local-04`, branch `refs/heads/claude/tddft-hessian-relaxed-density-fix-20260904`, HEAD `9795db633e0eb4c7e03efea9115d41ca0f10effb`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-s1-freq-20260907-01a07a3a

- 원본 `local-04`, branch `refs/heads/codex/mrsf-s1-freq-20260907-01a07a3a`, HEAD `30d8ed9c7cc1350ec460ba7c936c971f4730dcb6`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-space-separated-input-20260829

- 원본 `local-04`, branch `refs/heads/codex/space-separated-method-basis-20260829`, HEAD `dff991c29e507a21d4806b78833bb1cfecc70c50`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/ChatGPT/OpenQP-tddft-hessian-20260831

- 원본 `local-04`, branch `refs/heads/codex/import-tddft-hessian-20260831`, HEAD `576be1a7a3d0707f5ae16b63cec20676407bc133`
- 존재: True; status 항목: 1

```text
?? .calc/
```
## /Users/cheolhochoi/Library/Caches/openqp/benchmarks/pr391-speed-576be1a7/source

- 원본 `local-04`, branch `refs/heads/codex/pr391-speed-benchmark-20260902`, HEAD `576be1a7a3d0707f5ae16b63cec20676407bc133`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/mrsf-curvature-origin-20260920-01a07f50

- 원본 `local-04`, branch `refs/heads/codex/mrsf-curvature-origin-20260920-01a07f50`, HEAD `dbd2214df861913e592e59d3b885f2f07488fda3`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/mrsf-curvature-postfix-20260920-01a07f50

- 원본 `local-04`, branch `refs/heads/codex/mrsf-curvature-postfix-20260920-01a07f50`, HEAD `a0b9d20ab724a6198448e40c095205742e60da88`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/oqp-builds/mrsf-final-4fb453d1-20260901/source

- 원본 `local-04`, branch `detached`, HEAD `4fb453d11f61a960a33515d6415dd931f7abf19b`
- 존재: True; status 항목: 2

```text
?? tests/mrsf_h2/h2_mrsftdhf_hessian_rohf.freq.molden
?? tests/mrsf_h2/h2_mrsftdhf_hessian_rohf.hess.json
```
## /Users/cheolhochoi/oqp-builds/mrsf-final-7c73f6fd-20260901/source

- 원본 `local-04`, branch `refs/heads/codex/mrsf-final-build-7c73f6fd-20260901`, HEAD `7c73f6fdb8aafd2e99744b491df985b406abd75e`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/oqp-builds/mrsf-final-99f5df16-20260901/source

- 원본 `local-04`, branch `detached`, HEAD `99f5df1675f1e818b5eb6f3d2bd7f46f05889dcb`
- 존재: True; status 항목: 2

```text
?? tests/mrsf_h2/h2_mrsftdhf_hessian_rohf.freq.molden
?? tests/mrsf_h2/h2_mrsftdhf_hessian_rohf.hess.json
```
## /Users/cheolhochoi/oqp-builds/mrsf-hessian-5f48567e/source

- 원본 `local-04`, branch `detached`, HEAD `5f48567e9e40c3869290ecd79cb9e7b0a2664dde`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/oqp-builds/mrsf-ir-raman-956d2c87-20260901/source

- 원본 `local-04`, branch `detached`, HEAD `956d2c8732837e2ad3f8b08e014b2becac27f990`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp

- 원본 `local-05`, branch `refs/heads/codex/openqp-dftb-cmake-hook`, HEAD `2b08607666ca813eefe231490fe66973d86a870b`
- 존재: True; status 항목: 106

```text
 M include/oqp.h
 M pyoqp/oqp/__init__.py
 M pyoqp/oqp/library/liboqp.py
 M pyoqp/oqp/library/libscipy.py
 M pyoqp/oqp/library/oqp_engine.py
 M pyoqp/oqp/library/single_point.py
 M pyoqp/oqp/molecule/oqpdata.py
 M pyoqp/oqp/openqp.py
 M pyoqp/oqp/runtime.py
 M pyoqp/oqp/utils/input_checker.py
 M source/modules/get_states_overlap.F90
 M source/modules/tdhf_mrsf_z_vector.F90
 M source/modules/tdhf_sf_z_vector.F90
 M source/modules/tdhf_z_vector.F90
 M source/tdhf_lib.F90
 M source/tdhf_sf_lib.F90
 M source/zvector_common.F90
 M tests/test_davidson_solver_stability.py
 M tests/test_openqp_dftb_cmake_integration.py
 M tests/test_oqp_optimizer.py
 M tests/test_zvector_solver_stability.py
?? ".dockerignore 2"
?? ".github/workflows/claude 2.yml"
?? "GRAD_SCREENING_NOTES 2.md"
?? "MRSF_ZVECTOR_PERF_NOTES 2.md"
?? "cmake/FindDDX 2.cmake"
?? "cmake/patches/ddx-v0.8.0-ilp64 2.patch"
?? "cmake/patches/libtagarray-v0.0.6-default-integer-shapes 2.patch"
?? "docs/coarse_to_fine_xc_grid 2.md"
?? "docs/dft_xc_reuse 2.md"
?? "docs/perf_levels 2.md"
?? "docs/progressive_screening 2.md"
?? "docs/release-packaging 2.md"
?? "examples/HESS/H2O_RHF-DFT_ANA_HESS 2.inp"
?? "examples/HESS/H2O_RHF-DFT_ANA_HESS 2.json"
?? "examples/HESS/H2O_RHF-DFT_ANA_HESS.hess 2.json"
?? "examples/ISPHER/H2O_BHHLYP-MRSFTDDFT_GRADIENT_F_SHELL_ISPHER 2.inp"
?? "examples/ISPHER/H2O_BHHLYP-MRSFTDDFT_GRADIENT_ISPHER 2.inp"
?? "examples/ISPHER/H2O_RHF-BHHLYP_ANALYTIC_HESSIAN_ISPHER 2.inp"
?? "examples/ISPHER/H2O_ROHF-MRSF-EKT_EA_ISPHER 2.inp"
?? "examples/ISPHER/H2O_ROHF-MRSF-EKT_IP_ISPHER 2.inp"
?? "examples/ISPHER/HBr_RHF-BHHLYP_ECP_GRADIENT_ISPHER 2.inp"
?? "examples/NMR/CH3_ROHF-DFT-PBE0-GIAO-NMR 2.inp"
?? "examples/NMR/CH3_ROHF-DFT-PBE0-GIAO-NMR 2.json"
?? "examples/NMR/CH3_ROHF-GIAO-NMR 2.inp"
?? "examples/NMR/CH3_ROHF-GIAO-NMR 2.json"
?? "examples/NMR/CH3_UHF-DFT-PBE0-GIAO-NMR 2.inp"
?? "examples/NMR/CH3_UHF-DFT-PBE0-GIAO-NMR 2.json"
?? "examples/NMR/CH3_UHF-GIAO-NMR 2.inp"
?? "examples/NMR/CH3_UHF-GIAO-NMR 2.json"
?? "examples/NMR/H2O_RHF-DFT-PBE0-GIAO-NMR 2.inp"
?? "examples/NMR/H2O_RHF-DFT-PBE0-GIAO-NMR 2.json"
?? "examples/NMR/H2O_RHF-GIAO-NMR 2.inp"
?? "examples/NMR/H2O_RHF-GIAO-NMR 2.json"
?? "examples/NMR/H2O_RHF-NMR 2.inp"
?? "examples/NMR/H2O_RHF-NMR 2.json"
?? "examples/OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP 2.inp"
?? "examples/OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP 2.json"
?? "examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP 2.inp"
?? "examples/OPT/HCN_RHF-DFT_IRC_OQP 2.inp"
?? "examples/OPT/HCN_RHF-DFT_NEB_OQP 2.inp"
?? "examples/OPT/HCN_RHF-DFT_NEB_OQP_product 2.xyz"
?? "examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER 2.inp"
?? "examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER 2.json"
?? "examples/PCM/OH_ROHF-HF_DDPCM_ENERGY 2.inp"
?? "examples/PCM/OH_ROHF-HF_DDPCM_ENERGY 2.json"
?? "examples/PERF/H2O_DFT_perf1 2.inp"
?? "examples/PERF/H2O_MRSF_perf1 2.inp"
?? "examples/PERF/H2O_MRSF_perf2_override 2.inp"
?? "examples/PROP/H2O_BHHLYP_PROPERTIES 2.inp"
?? "examples/PROP/H2O_CATION_UHF_PROPERTIES 2.inp"
?? "examples/PROP/H2O_CATION_UHF_PROPERTIES 2.json"
?? "examples/PROP/H2O_RHF_PROPERTIES 2.inp"
?? "examples/PROP/H2O_RHF_PROPERTIES 2.json"
?? "examples/QUANTUM/h2 2.inp"
?? "examples/SOC/CH3Br-BHHLYP-SOC 2.inp"
?? "examples/SOC/CH3Br-BHHLYP-SOC 2.json"
?? "examples/SOC/H2O_BHHLYP_SOC 2.inp"
?? "examples/SOC/H2O_BHHLYP_SOC 2.json"
?? "examples/XAS/README 2.md"
?? "examples/XAS/delta-CHP-MRSF/HCN_CHP-MRSF.inp 2.example"
?? "examples/XAS/delta-CHP-MRSF/HCN_CHP-MRSF.json 2.example"
?? "examples/other/h2o_rhf_3-21g_minao 2.inp"
?? "examples/other/h2o_rhf_3-21g_minao 2.json"
?? "examples/other/h2o_rhf_3-21g_modhuckel 2.inp"
?? "examples/other/h2o_rhf_3-21g_modhuckel 2.json"
?? "examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp 2.inp"
?? "examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp 2.json"
?? "examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp 2.inp"
?? "examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp 2.json"
?? pyoqp/oqp/library/openqp_dftb.py
?? "pyoqp/oqp/library/scf_selector_model 2.py"
?? "pyoqp/oqp/quantum/fcidump 2.py"
?? source/openqp_dftb_bridge.F90
?? source/response_solver_common.F90
?? "tests/smoke_response_symmetry_validation 2.py"
?? "tests/test_ddx_scf_integration_seam 2.py"
?? "tests/test_mrsf_pure_hf_spc_scale 2.py"
?? "tests/test_nmr_giao_para_live 2.py"
?? "tests/test_openqp_dftb_cmake_integration 2.py"
?? "tests/test_openqp_dftb_cmake_integration 3.py"
?? "tests/test_openqp_dftb_cmake_integration 4.py"
?? tests/test_openqp_dftb_schema_hooks.py
?? "tests/test_perf_levels 2.py"
?? tmp/
?? work/
```
## /Users/cheolhochoi/.codex/worktrees/164d/openqp

- 원본 `local-05`, branch `detached`, HEAD `e9ff159da9624975e939602f788a8eb4b533372d`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/.codex/worktrees/7021/openqp

- 원본 `local-05`, branch `refs/heads/codex/molecular-symmetry-validate`, HEAD `ee20e381152c638786e0d95218dee0b2d6537dc7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/.codex/worktrees/9d05/openqp

- 원본 `local-05`, branch `detached`, HEAD `d5bd9724e420bf1515c223547dc6ee8bbebee370`
- 존재: True; status 항목: 8

```text
 M pyoqp/oqp/library/single_point.py
 M pyoqp/oqp/molecule/oqpdata.py
 M pyoqp/oqp/utils/input_checker.py
 M source/scf.F90
 M tests/test_single_point_scf_fallback.py
?? pyoqp/oqp/library/scf_selector_model.py
?? tests/test_scf_active_manager.py
?? tests/test_scf_selector_manager_modes.py
```
## /Users/cheolhochoi/.codex/worktrees/b8e4/openqp

- 원본 `local-05`, branch `refs/heads/codex/soc-namd-options`, HEAD `7b1bb7597dea35c81cd0e65f59d4a97151d25a7b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/.codex/worktrees/d435/openqp

- 원본 `local-05`, branch `detached`, HEAD `e9ff159da9624975e939602f788a8eb4b533372d`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/.codex/worktrees/db66/openqp

- 원본 `local-05`, branch `refs/heads/codex/native-trah-default`, HEAD `5f87db071b419280ba21fbb1e177422d5c894683`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp-ovexact-fix

- 원본 `local-05`, branch `refs/heads/feat/dftb-conventional-default`, HEAD `db391506c82090739d024210d378ee0cff605643`
- 존재: True; status 항목: 4

```text
?? openqp_soc_test_tmp_2026-07-26_10-51-18/
?? openqp_soc_test_tmp_2026-07-26_10-51-20/
?? openqp_soc_test_tmp_2026-07-26_10-55-32/
?? openqp_soc_test_tmp_2026-07-26_10-55-35/
```
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp-reks22

- 원본 `local-05`, branch `refs/heads/feature/reks22-scf`, HEAD `53d1d7e8f1ab5dedb9efcbb4a013c27e0930eaf0`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp-xtb-pyoqp

- 원본 `local-05`, branch `refs/heads/feat/xtb-log-overhaul`, HEAD `022016fdecf8dd7baf9f7d13093743d444d61133`
- 존재: True; status 항목: 391

```text
M  .github/workflows/CI.yml
M  README.md
A  docs/native_optimization_recovery.md
A  examples/DFT/H2O_RHF-DFT_ENERGY.oqp
A  examples/DFT/H2O_RHF-DFT_GRADIENT.oqp
A  examples/DFT/H2O_ROHF-DFT_ENERGY.oqp
A  examples/DFT/H2O_ROHF-DFT_GRADIENT.oqp
A  examples/DFT/H2O_UHF-DFT_ENERGY.oqp
A  examples/DFT/H2O_UHF-DFT_GRADIENT.oqp
M  examples/DFTB/C2_LC-MRSF-TDDFTB_ERF-TUNED_ENERGY.inp
A  examples/DFTB/C2_LC-MRSF-TDDFTB_ERF-TUNED_ENERGY.oqp
A  examples/DFTB/CH2_MRSF-TDDFTB_ENERGY.oqp
A  examples/DFTB/CH2_SF-TDDFTB_ENERGY.oqp
M  examples/DFTB/H2CH2_MRSF-TDDFTB_GRADIENT.inp
A  examples/DFTB/H2CH2_MRSF-TDDFTB_GRADIENT.oqp
A  examples/DFTB/H2O_TDDFTB_ENERGY.oqp
A  examples/DFTB/H2_DFTB_ENERGY.oqp
A  examples/DFTB/README.md
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_ENERGY.inp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_ENERGY.oqp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_EXPLICIT_ENERGY.inp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_EXPLICIT_ENERGY.oqp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_S1_GRADIENT.inp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_S1_GRADIENT.oqp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_S1_OPTIMIZE.inp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_S1_OPTIMIZE.oqp
A  examples/DTCAM-TB/BUTADIENE_DTCAMTB_S1_OPTIMIZE_GEOMETRIC.oqp
A  examples/DTCAM-TB/ETHYLENE_DTCAMTB_S1S0S2_MECI_BAEKA.oqp
A  examples/DTCAM-TB/ETHYLENE_DTCAMTB_S1S0_MECI.inp
A  examples/DTCAM-TB/ETHYLENE_DTCAMTB_S1S0_MECI.oqp
A  examples/DTCAM-TB/README.md
A  examples/ECP/C2H4_BHHLYP-MRSFTDDFT_Energy.oqp
A  examples/ECP/C2H4_BHHLYP-MRSFTDDFT_Grad.oqp
A  examples/ECP/Custom_Basis-Set/C2H4_Custom_Basis-Set-MRSFTDDFT_Energy.oqp
A  examples/ECP/HBr_BHHLYP-MRSFTDDFT_ENERGY.oqp
A  examples/ECP/HBr_RHF-DFT_ENERGY.oqp
A  examples/ECP/HBr_RHF-DFT_GRADIENT.oqp
A  examples/ECP/Mg_RHF-DFT_ENERGY.oqp
A  examples/ECP/NaCl_BHHLYP-MRSFTDDFT_ENERGY.oqp
A  examples/ECP/NaCl_RHF-DFT_GRADIENT.oqp
M  examples/HESS/H2O_BHHLYP-MRSFTDDFT_NUM_HESS.hess.json
A  examples/HESS/H2O_BHHLYP-MRSFTDDFT_NUM_HESS.oqp
M  examples/HESS/H2O_RHF-DFT_ANA_HESS.hess.json
A  examples/HESS/H2O_RHF-DFT_ANA_HESS.oqp
M  examples/HESS/H2O_RHF-DFT_NUM_HESS.hess.json
A  examples/HESS/H2O_RHF-DFT_NUM_HESS.oqp
A  examples/HF/H2O_RHF-HF_ENERGY.oqp
A  examples/HF/H2O_RHF-HF_GRADIENT.oqp
A  examples/HF/H2O_ROHF-HF_ENERGY.oqp
A  examples/HF/H2O_ROHF-HF_GRADIENT.oqp
A  examples/HF/H2O_UHF-HF_ENERGY.oqp
A  examples/HF/H2O_UHF-HF_GRADIENT.oqp
A  examples/ISPHER/H2O_BHHLYP-MRSFTDDFT_GRADIENT_F_SHELL_ISPHER.oqp
A  examples/ISPHER/H2O_BHHLYP-MRSFTDDFT_GRADIENT_ISPHER.oqp
A  examples/ISPHER/H2O_RHF-BHHLYP_ANALYTIC_HESSIAN_ISPHER.oqp
A  examples/ISPHER/H2O_ROHF-MRSF-EKT_EA_ISPHER.oqp
A  examples/ISPHER/H2O_ROHF-MRSF-EKT_IP_ISPHER.oqp
A  examples/ISPHER/HBr_RHF-BHHLYP_ECP_GRADIENT_ISPHER.oqp
A  examples/MP2/h2o_ump2_6-31g.oqp
A  examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.oqp
A  examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.oqp
A  examples/NMR/CH3_ROHF-DFT-PBE0-GIAO-NMR.oqp
A  examples/NMR/CH3_ROHF-GIAO-NMR.oqp
A  examples/NMR/CH3_UHF-DFT-PBE0-GIAO-NMR.oqp
A  examples/NMR/CH3_UHF-GIAO-NMR.oqp
A  examples/NMR/H2O_RHF-DFT-PBE0-GIAO-NMR.oqp
A  examples/NMR/H2O_RHF-GIAO-NMR.oqp
A  examples/NMR/H2O_RHF-NMR.oqp
R  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECI_GEOMETRIC.inp -> examples/OPT/C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.inp
R  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECP_GEOMETRIC.json -> examples/OPT/C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.json
A  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.oqp
A  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECI.oqp
D  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECI_GEOMETRIC.json
R  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECP_GEOMETRIC.inp -> examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECP_OQP.inp
A  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECP_OQP.json
A  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECP_OQP.oqp
M  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MEP.inp
M  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MEP.json
A  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MEP.oqp
A  examples/OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.oqp
A  examples/OPT/H2O_BHHLYP-MRSFTDDFT_OPTIMIZE.oqp
A  examples/OPT/H2O_RHF-DFT_OPTIMIZE.oqp
D  examples/OPT/H2O_RHF-DFT_OPTIMIZE_GEOMETRIC.inp
D  examples/OPT/H2O_RHF-DFT_OPTIMIZE_GEOMETRIC.json
A  examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.json
A  examples/OPT/H2O_RHF-DFT_OPTIMIZE_OQP.oqp
R  examples/OPT/HCN_BHHLYP-MRSFTDDFT_TS_GEOMETRIC.inp -> examples/OPT/HCN_BHHLYP-MRSFTDDFT_TS_OQP.inp
R  examples/OPT/HCN_BHHLYP-MRSFTDDFT_TS_GEOMETRIC.json -> examples/OPT/HCN_BHHLYP-MRSFTDDFT_TS_OQP.json
A  examples/OPT/HCN_BHHLYP-MRSFTDDFT_TS_OQP.oqp
D  examples/OPT/HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.constraints
D  examples/OPT/HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.inp
D  examples/OPT/HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.json
R  examples/OPT/HCN_RHF-DFT_IRC_GEOMETRIC.inp -> examples/OPT/HCN_RHF-DFT_CONSTRAINED_OQP.inp
A  examples/OPT/HCN_RHF-DFT_CONSTRAINED_OQP.json
A  examples/OPT/HCN_RHF-DFT_CONSTRAINED_OQP.oqp
D  examples/OPT/HCN_RHF-DFT_IRC_GEOMETRIC.json
A  examples/OPT/HCN_RHF-DFT_IRC_OQP.json
A  examples/OPT/HCN_RHF-DFT_IRC_OQP.oqp
A  examples/OPT/HCN_RHF-DFT_NEB_OQP.json
A  examples/OPT/HCN_RHF-DFT_NEB_OQP.oqp
D  examples/OPT/HCN_RHF-DFT_TS_GEOMETRIC.json
R  examples/OPT/HCN_RHF-DFT_TS_GEOMETRIC.inp -> examples/OPT/HCN_RHF-DFT_TS_OQP.inp
A  examples/OPT/HCN_RHF-DFT_TS_OQP.json
A  examples/OPT/HCN_RHF-DFT_TS_OQP.oqp
A  examples/OQP_INPUT/H2O_MRSF_S0_OPT.oqp
A  examples/OQP_INPUT/h2o.xyz
A  examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.oqp
A  examples/PCM/OH_ROHF-HF_DDPCM_ENERGY.oqp
A  examples/PERF/H2O_DFT_perf1.oqp
A  examples/PERF/H2O_MRSF_perf1.oqp
A  examples/PERF/H2O_MRSF_perf2_override.oqp
A  examples/PROP/H2O_BHHLYP_PROPERTIES.oqp
A  examples/PROP/H2O_CATION_UHF_PROPERTIES.oqp
A  examples/PROP/H2O_RHF_PROPERTIES.oqp
A  examples/QMMM/2E4E_RHF-DFT-QMMM_energy.oqp
A  examples/QMMM/H2CO-water_BHHLYP-MRSF-NAMD-QMMM.oqp
A  examples/QMMM/H2CO-water_BHHLYP-SOC-NAMD-QMMM.oqp
A  examples/QMMM/ala-dipeptide_BHHLYP-QMMM-MD-RCD.oqp
A  examples/QMMM/ala.oqp
M  examples/QMMM/run.inp
A  examples/QMMM/run.oqp
A  examples/QUANTUM/h2.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_adiis.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_basic.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_diis-reset.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_ediis.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_mom.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_pfon.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_soscf.oqp
A  examples/SCF/h2o_rhf_6-31g_pbe_vdiis.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_adiis.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_basic.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_diis-reset.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_ediis.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_mom.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_pfon.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_soscf.oqp
A  examples/SCF/h2o_rohf_6-31g_pbe_vshift.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_adiis.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_basic.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_diis-reset.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_ediis.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_mom.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_pfon.oqp
A  examples/SCF/h2o_uhf-s_6-31g_pbe_vdiis.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_adiis.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_basic.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_diis-reset.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_ediis.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_mom.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_pfon.oqp
A  examples/SCF/h2o_uhf-t_6-31g_pbe_vdiis.oqp
A  examples/SCF/h2o_uhf_6-31g_pbe_soscf.oqp
A  examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_ENERGY.oqp
A  examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_GRADIENT.oqp
A  examples/SOC/CH3Br-BHHLYP-SOC.oqp
A  examples/SOC/H2O_BHHLYP_SOC.oqp
A  examples/TDDFT/H2O_B3LYP5-TDDFT_ENERGY.oqp
A  examples/TDDFT/H2O_B3LYP5-TDDFT_GRADIENT.oqp
A  examples/TDHF/H2O_TDHF_ENERGY.oqp
A  examples/TDHF/H2O_TDHF_GRADIENT.oqp
A  examples/TRAH/H2O_RHF-DFT_GRAD_TRAH.oqp
A  examples/TRAH/H2O_ROHF-DFT_GRAD_TRAH.oqp
A  examples/TRAH/H2O_UHF-DFT_GRAD_TRAH.oqp
A  examples/TRAH/h2o_rohf_mrsf-t_6-31g_prop.oqp
A  examples/UMRSF-TDDFT/C4H6_BHHLYP_UMRSFTDDFT_ENERGY.oqp
A  examples/XAS/MRSF/HCN_MRSF.oqp
A  examples/XAS/delta-CHP-MRSF/HCN_CHP-MRSF_v2.oqp
A  examples/XAS/delta-CHP-MRSF/HCN_DFT.oqp
A  examples/geometries/BrH-a6a9d6f6a6f3.xyz
A  examples/geometries/C2-d0affeda5215.xyz
A  examples/geometries/C2H4-25af03d76164.xyz
A  examples/geometries/C2H4-41911fceceb6.xyz
A  examples/geometries/C2H4-5ec6cdbfac9d.xyz
A  examples/geometries/C2H4-c7ec226122c3.xyz
A  examples/geometries/C4H6-0a0dabfd2a8a.xyz
A  examples/geometries/C4H6-e783f4d884c7.xyz
A  examples/geometries/CH2-e41258125cee.xyz
A  examples/geometries/CH2O-2bc62dda4b8a.xyz
A  examples/geometries/CH3-0b060bed6a1d.xyz
A  examples/geometries/CH3Br-82fdf47b3e2c.xyz
A  examples/geometries/CHN-682f5d7d753a.xyz
A  examples/geometries/CHN-825a4a518dc8.xyz
A  examples/geometries/CHN-cef1bbffc7f3.xyz
A  examples/geometries/CHN-d957be87de22.xyz
A  examples/geometries/ClNa-0ed89ad56c59.xyz
A  examples/geometries/H2-069e977a383f.xyz
A  examples/geometries/H2O-0381125c86f2.xyz
A  examples/geometries/H2O-1bdadad53a3a.xyz
A  examples/geometries/H2O-781cba5cd72b.xyz
A  examples/geometries/H2O-7c17dce6a2e9.xyz
A  examples/geometries/H2O-95fa99c614ed.xyz
A  examples/geometries/H2O-a5185d3fce44.xyz
A  examples/geometries/H4O2-6c06eba36ab5.xyz
A  examples/geometries/HO-0e90d87efea2.xyz
A  examples/geometries/Mg-21a9b71503dd.xyz
A  examples/other/h2o-2_rhf_cc-pvtz_hf.oqp
A  examples/other/h2o-2_rohf_cc-pvtz_hf.oqp
A  examples/other/h2o-2_uhf-s_cc-pvtz_hf.oqp
A  examples/other/h2o_nacme_rohf_mrsf-s_6-31g_bhhlyp.oqp
A  examples/other/h2o_rhf_3-21g_minao.oqp
A  examples/other/h2o_rhf_3-21g_modhuckel.oqp
A  examples/other/h2o_rhf_3-21g_sap.oqp
A  examples/other/h2o_rhf_6-31g_b3lypv5.oqp
A  examples/other/h2o_rhf_6-31g_bhhlyp.oqp
A  examples/other/h2o_rhf_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rhf_6-31g_hf.oqp
A  examples/other/h2o_rhf_6-31g_m06-2x.oqp
A  examples/other/h2o_rhf_6-31g_pbe.oqp
A  examples/other/h2o_rhf_6-31g_slater.oqp
A  examples/other/h2o_rhf_cc-pvtz_b3lypv5.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_b3lypv5.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_bhhlyp.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_hf.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_m06-2x.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_pbe.oqp
A  examples/other/h2o_rhf_rpa-s_6-31g_slater.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_b3lypv5.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_bhhlyp.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_hf.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_m06-2x.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_pbe.oqp
A  examples/other/h2o_rhf_rpa-t_6-31g_slater.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_b3lypv5.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_bhhlyp.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_hf.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_m06-2x.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_pbe.oqp
A  examples/other/h2o_rhf_tda-s_6-31g_slater.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_b3lypv5.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_bhhlyp.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_hf.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_m06-2x.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_pbe.oqp
A  examples/other/h2o_rhf_tda-t_6-31g_slater.oqp
A  examples/other/h2o_rohf-dft_energy_init_scf.oqp
A  examples/other/h2o_rohf-dft_energy_init_scf_basis_library.oqp
A  examples/other/h2o_rohf_6-31g_b3lypv5.oqp
A  examples/other/h2o_rohf_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rohf_6-31g_hf.oqp
A  examples/other/h2o_rohf_6-31g_m06-2x.oqp
A  examples/other/h2o_rohf_6-31g_pbe.oqp
A  examples/other/h2o_rohf_6-31g_slater.oqp
A  examples/other/h2o_rohf_cc-pvtz_b3lypv5.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_b3lypv5.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dt-bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dt-vee.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-aee.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-b3lyp.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-stg.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-tune.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-vaee.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-vee.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-xi.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-xiv.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_hf.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_m06-2x.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_pbe.oqp
A  examples/other/h2o_rohf_mrsf-q_6-31g_slater.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_b3lypv5.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp-spc-coco.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp-spc-coov.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp-spc-ovov.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dt-bhhlyp-spc.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dt-bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dt-vee.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-aee.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-b3lyp.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-stg.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-tune.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-vaee.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-vee.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-xi.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-xiv.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_hf.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_m06-2x.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_pbe.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_prop.oqp
A  examples/other/h2o_rohf_mrsf-s_6-31g_slater.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_b3lypv5.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_dt-bhhlyp-spc.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_dt-bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_dt-vee.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_dtcam-b3lyp.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_hf.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_m06-2x.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_pbe.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_slater.oqp
A  examples/other/h2o_rohf_mrsf-t_6-31g_stg1x.oqp
A  examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_sf_6-31g_b3lypv5.oqp
A  examples/other/h2o_rohf_sf_6-31g_bhhlyp.oqp
A  examples/other/h2o_rohf_sf_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_rohf_sf_6-31g_dt-bhhlyp.oqp
A  examples/other/h2o_rohf_sf_6-31g_dtcam-b3lyp.oqp
A  examples/other/h2o_rohf_sf_6-31g_hf.oqp
A  examples/other/h2o_rohf_sf_6-31g_m06-2x.oqp
A  examples/other/h2o_rohf_sf_6-31g_pbe.oqp
A  examples/other/h2o_rohf_sf_6-31g_slater.oqp
A  examples/other/h2o_uhf-s_6-31g_b3lypv5.oqp
A  examples/other/h2o_uhf-s_6-31g_bhhlyp.oqp
A  examples/other/h2o_uhf-s_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_uhf-s_6-31g_hf.oqp
A  examples/other/h2o_uhf-s_6-31g_m06-2x.oqp
A  examples/other/h2o_uhf-s_6-31g_pbe.oqp
A  examples/other/h2o_uhf-s_6-31g_slater.oqp
A  examples/other/h2o_uhf-s_cc-pvtz_b3lypv5.oqp
A  examples/other/h2o_uhf-t_6-31g_b3lypv5.oqp
A  examples/other/h2o_uhf-t_6-31g_bhhlyp.oqp
A  examples/other/h2o_uhf-t_6-31g_cam-b3lyp.oqp
A  examples/other/h2o_uhf-t_6-31g_hf.oqp
A  examples/other/h2o_uhf-t_6-31g_m06-2x.oqp
A  examples/other/h2o_uhf-t_6-31g_pbe.oqp
A  examples/other/h2o_uhf-t_6-31g_slater.oqp
M  pyoqp/README.md
A  pyoqp/oqp/library/baeka.py
M  pyoqp/oqp/library/liboqp.py
M  pyoqp/oqp/library/libscipy.py
UU pyoqp/oqp/library/namd.py
M  pyoqp/oqp/library/neb_utils.py
UU pyoqp/oqp/library/openqp_dftb.py
M  pyoqp/oqp/library/oqp_engine.py
M  pyoqp/oqp/library/oqp_irc.py
M  pyoqp/oqp/library/oqp_neb.py
UU pyoqp/oqp/library/qmmm_driver.py
M  pyoqp/oqp/library/runfunc.py
M  pyoqp/oqp/library/set_basis.py
M  pyoqp/oqp/library/single_point.py
M  pyoqp/oqp/molecule/molecule.py
M  pyoqp/oqp/molecule/oqpdata.py
UU pyoqp/oqp/openqp.py
M  pyoqp/oqp/pyoqp.py
A  pyoqp/oqp/utils/dftb_trace.py
M  pyoqp/oqp/utils/file_utils.py
UU pyoqp/oqp/utils/input_checker.py
M  pyoqp/oqp/utils/input_parser.py
A  pyoqp/oqp/utils/json_utils.py
A  pyoqp/oqp/utils/oqp_input.py
UU pyoqp/oqp/utils/oqp_tester.py
M  pyoqp/oqp/utils/qmmm.py
M  pyoqp/oqp/utils/regression.py
A  pyoqp/oqp/utils/state_labels.py
M  pyoqp/setup.py
M  pyproject.toml
M  source/modules/get_states_overlap.F90
M  source/modules/soc_mrsf.F90
M  source/modules/tdhf_mrsf_energy.F90
M  source/modules/tdhf_mrsf_gradient.F90
M  source/modules/tdhf_mrsf_z_vector.F90
M  source/scf.F90
M  source/tdhf_sf_lib.F90
M  tests/test_advanced_guess_examples.py
M  tests/test_analytic_hessian.py
M  tests/test_analytic_hessian_bindings.py
M  tests/test_davidson_solver_stability.py
A  tests/test_dftb_model_default.py
A  tests/test_dftb_parameter_defaults.py
A  tests/test_dftb_trace.py
M  tests/test_geometric_optimizer.py
A  tests/test_input_geometry_references.py
A  tests/test_json_td_layout.py
A  tests/test_legacy_example_conversion.py
M  tests/test_mrsf_ekt_scaffold.py
M  tests/test_neb_input.py
A  tests/test_openmm_optional_backend.py
UU tests/test_openqp_api.py
A  tests/test_openqp_dftb_abi1.py
M  tests/test_openqp_dftb_api.py
A  tests/test_openqp_dftb_logging.py
A  tests/test_oqp_input.py
A  tests/test_oqp_input_schema_manifest.py
A  tests/test_oqp_neb_native.py
M  tests/test_oqp_optimizer.py
A  tests/test_set_basis_tagged_geometry.py
M  tests/test_single_point_scf_fallback.py
M  tests/test_soc_namd_qmmm_production.py
A  tests/test_state_labels.py
M  tests/test_symmetry_metadata.py
A  tools/convert_legacy_examples.py
```
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp-xtb-qmmm

- 원본 `local-05`, branch `refs/heads/feat/xtb-qmmm-fullespf`, HEAD `b217f46d57853ed2ae046a05363dc6f7bac4e336`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp/.claude/worktrees/determined-pasteur-f56cf9

- 원본 `local-05`, branch `detached`, HEAD `f2a24c1ce97077e2b286e0d4baf121526cc40161`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp/.claude/worktrees/elastic-wozniak-9b6f5b

- 원본 `local-05`, branch `detached`, HEAD `f2a24c1ce97077e2b286e0d4baf121526cc40161`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-24/openqp/.claude/worktrees/objective-poitras-2113a3

- 원본 `local-05`, branch `detached`, HEAD `f2a24c1ce97077e2b286e0d4baf121526cc40161`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-05-31/openqp-solvent-backend-spike

- 원본 `local-05`, branch `refs/heads/feat/mrsf-pcm-spike`, HEAD `89b5e312ee085ba45d84fd496e2ea401614eb3b1`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-06-27/openqp-pr234-review

- 원본 `local-05`, branch `detached`, HEAD `7e411732d9e3029165bc747b4e4cda5fc1e3eef2`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-06-29/openqp-claude-all-pr

- 원본 `local-05`, branch `refs/heads/codex/claude-all-pr-review`, HEAD `70a663f0fe01a9d62e64238ee2caec578a3eff4b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-07-07/mrsf-tddftb-manuscript-continuation/work/openqp-pr-dftb-native

- 원본 `local-05`, branch `refs/heads/codex/openqp-dftb-native`, HEAD `929ea903b7a9ba7a45a261dab55a8ef56dc2bcb6`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-07-07/new-chat/work/openqp-upstream-head-20260707

- 원본 `local-05`, branch `detached`, HEAD `8c30c85d68a19ad107d6a5f9f7abbfe671f958a8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-02/openqp-mo-frequency-export

- 원본 `local-05`, branch `refs/heads/agent/export-mo-frequency-formats`, HEAD `57a70ade873754bbea04178ae19a1c4ee7bfce50`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-03/openqp-finite-water-droplet

- 원본 `local-05`, branch `refs/heads/agent/finite-water-droplet-boundary`, HEAD `f5c755db59a2da8490785a34cbdb38cdb1c0e343`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-03/openqp-odp

- 원본 `local-05`, branch `refs/heads/agent/odp-umbrella-native`, HEAD `717ecfbb631227d2c1cff72da1ea851b05e94c3b`
- 존재: True; status 항목: 1

```text
?? work/
```
## /Users/cheolhochoi/Documents/Codex/2026-08-08/openqp-dlc-native-default

- 원본 `local-05`, branch `refs/heads/codex/dlc-native-default`, HEAD `b88b7fb0e8e9c5cb843c3b9b73c645bfc87835ea`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-12/openqp-xms-caspt2-molcas-benchmark

- 원본 `local-05`, branch `refs/heads/codex/xms-caspt2-molcas-benchmark-20260812`, HEAD `9bd4d64f588059111b9e258a11c922b34f8f5473`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-12/openqp-xms-caspt2-molcas-upstream

- 원본 `local-05`, branch `refs/heads/codex/xms-caspt2-molcas-upstream-20260812`, HEAD `28a6e2f41eee319e6b8756e3cc587ca768b5da39`
- 존재: True; status 항목: 1

```text
?? benchmarks/
```
## /Users/cheolhochoi/Documents/claude/gradient-fix/mrsf-nmr-clean

- 원본 `local-06`, branch `refs/heads/feat/mrsf-nmr-gate5a`, HEAD `71e064cc4e01792137af59cc53fac51f513b10be`
- 존재: True; status 항목: 10

```text
 M source/atomic_structure.F90
 M tests/test_opentrustregion_linalg_config.py
?? build-abi-off/
?? pyoqp/oqp/__pycache__/
?? pyoqp/oqp/library/__pycache__/
?? pyoqp/oqp/molden/__pycache__/
?? pyoqp/oqp/molecule/__pycache__/
?? pyoqp/oqp/periodic_table/__pycache__/
?? pyoqp/oqp/utils/__pycache__/
?? tests/__pycache__/
```
## /Users/cheolhochoi/Documents/claude/gradient-fix/openqp

- 원본 `local-07`, branch `refs/heads/main`, HEAD `1890fefdb25dea17bb99dd38215b8e4e85b5a8fe`
- 존재: True; status 항목: 5

```text
 M source/tdhf_mrsf_lib.F90
?? ".github/workflows/claude 2.yml"
?? .venv/
?? build-romrsf/
?? build/
```
## /private/tmp/ekt_a64

- 원본 `local-07`, branch `detached`, HEAD `a64cffddf833d7bd3d11ad536193741a6712662c`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /private/tmp/gpu-xc-rebase

- 원본 `local-07`, branch `refs/heads/feat/gpu-xc-response-rebase`, HEAD `7f9eda08b20886137eadf8c7466158097196e057`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Users/cheolhochoi/Documents/claude/gradient-fix/ekt-mrsf-wt

- 원본 `local-07`, branch `refs/heads/research/ekt-mrsf-wtilde-derivation`, HEAD `9e2b94f7612092b98133d9eba7c02367744f7ab8`
- 존재: True; status 항목: 3

```text
?? build/
?? pyoqp/PyOpenQP.egg-info/
?? pyoqp/build/
```
## /Users/cheolhochoi/Documents/claude/gradient-fix/gpu-metc-wt

- 원본 `local-07`, branch `refs/heads/feat/gpu-metc-regression-test`, HEAD `fc92929bbdbbdf6c4afdcfffb851afad7dc63019`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/gradient-fix/gpu-workspace-wt

- 원본 `local-07`, branch `refs/heads/feat/gpu-workspace-manager`, HEAD `6413215eeb41ba63452ba641661acc9713f6af23`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/gradient-fix/nmr-wt

- 원본 `local-07`, branch `refs/heads/feat/mrsf-nmr-gate5a`, HEAD `4faaa117420e925b3a5a68f4afdeb7b349790d6b`
- 존재: True; status 항목: 3

```text
?? .venv
?? examples/HF/H2O_RHF-HF_ENERGY_scf_rhf_hf_6-31gs.molden
?? tests/__pycache__/
```
## /Users/cheolhochoi/Documents/claude/gradient-fix/openqp-pcm

- 원본 `local-07`, branch `refs/heads/feat/solvent-backend-spike-private-integrated`, HEAD `903f0eb4cd8d7dcef7f9cb7ddc20d970fb0493e4`
- 존재: True; status 항목: 1

```text
?? build/
```
## /Users/cheolhochoi/Documents/claude/openqp-private-opeqp-GPU

- 원본 `local-07`, branch `refs/heads/rot-seam-v12`, HEAD `9636a98b85dc9dbeea649629018757161d8a52f9`
- 존재: True; status 항목: 25

```text
?? ".git 2"
?? ".github 2/"
?? ".gitignore 2"
?? ".gitlab-ci 2.yml"
?? "CMakeLists 2.txt"
?? "Dockerfile 2"
?? "GPU_PORTING_UNIFIED_HANDOFF 2.md"
?? "LICENSE 2"
?? "MANIFEST 2.in"
?? "README 2.md"
?? "TROUBLESHOOTING 2.md"
?? "basis_sets 2/"
?? build-mac/
?? "cmake 2/"
?? "examples 2/"
?? "external 2/"
?? "include 2/"
?? install-mac/
?? "pyoqp 2/"
?? pyoqp/OpenQP.egg-info/
?? pyoqp/PyOpenQP.egg-info/
?? "pyproject 2.toml"
?? "source 2/"
?? "tests 2/"
?? "tools 2/"
```
## /Users/cheolhochoi/Documents/claude/openqp-hessian

- 원본 `local-08`, branch `refs/heads/claude/practical-hypatia-4tcBa`, HEAD `947781d78bead03da35952e5a5fae47fa434971c`
- 존재: True; status 항목: 2

```text
M  pyoqp/oqp/molecule/molecule.py
?? openqp_test_tmp_2026-05-30_18-44-09/
```
## /Users/cheolhochoi/Documents/claude/openqp-private

- 원본 `local-09`, branch `refs/heads/main`, HEAD `6ab4b86701bf972b56acc7c79bd40338a51b13ac`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-private-analytic-nac

- 원본 `local-09`, branch `refs/heads/feat/mrsf-analytic-nac`, HEAD `ecda4b84db966ff6804df378c6499ee8678f23f4`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-review

- 원본 `local-10`, branch `refs/heads/fix/oqp-timer-parser-blockers`, HEAD `901c3e2bc723c2a1bc8496fd265ba84d80d44aac`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp

- 원본 `local-11`, branch `refs/heads/feat/routec-dfjk-bridge`, HEAD `d9a752d4002a0a6dc4f5a19d1313fc087390f5a1`
- 존재: True; status 항목: 26

```text
 M source/modules/tdhf_mrsf_z_vector.F90
?? ".github/workflows/claude 2.yml"
?? __pycache__/
?? build-fix/
?? build-routec/
?? build/
?? fd_allstates.py
?? fd_classify.py
?? fd_grad_check.py
?? fd_one.py
?? fd_tables.py
?? h2o_rhf_3-21g_sad.sad.pyscf
?? h2o_rhf_3-21g_sap.sap.pyscf
?? iterate.sh
?? openqp_test_tmp_2026-05-30_13-59-21/
?? pyoqp/oqp/__pycache__/
?? pyoqp/oqp/library/__pycache__/
?? pyoqp/oqp/molden/__pycache__/
?? pyoqp/oqp/molecule/__pycache__/
?? pyoqp/oqp/periodic_table/__pycache__/
?? pyoqp/oqp/utils/__pycache__/
?? pyscf_wfn.json
?? "source/integrals/int2_routec 2.F90"
?? tests/__pycache__/
?? "tests/test_mrsf_pure_hf_spc_scale 2.py"
?? "tests/test_tda_gradient_zvector_hessian 2.py"
```
## /private/tmp/claude-501/-Users-cheolhochoi-Documents-claude-openqp/db4bb5d6-64f3-4637-9986-c4afd166f2d0/scratchpad/wt-241

- 원본 `local-11`, branch `refs/heads/resolve-241`, HEAD `868389cf4dd87eaa172eee9b1779c3ab50931e90`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Users/cheolhochoi/.codex/worktrees/openqp-h2-pr

- 원본 `local-11`, branch `refs/heads/fix/mrsf-h2-zero-closed`, HEAD `ee4ce62fd98f0f1697ea6351984a237789f9f95b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-delocate-wt

- 원본 `local-11`, branch `refs/heads/fix/macos-wheel-delocate-0.13`, HEAD `cb5c8509bbd69c42568a88c46afd977873d04ed3`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-intel-wt

- 원본 `local-11`, branch `refs/heads/feat/intel-oneapi-support`, HEAD `fa7ae05ca40719e9af7babfc9ed6aa4bb39aa374`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-numnac

- 원본 `local-11`, branch `refs/heads/bugfix-num-NAC`, HEAD `04e58a47ac4c89ed498e06ff00586671476ff45c`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-pr203

- 원본 `local-11`, branch `refs/heads/pr203-cache-key`, HEAD `3f7da4ad90af0cdd894dabd8feaf37c5116d14d5`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-stability-wt

- 원본 `local-11`, branch `refs/heads/fix/scf-stability-opt-in`, HEAD `90bdea5bf6a8953a36f4e1db3ba6e1b29f9558f8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp-windows-intel

- 원본 `local-11`, branch `refs/heads/fix/lf-license-files-only`, HEAD `b42262292d7bc7085f9f3889c3edfbf57d33f6bf`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp/.claude/worktrees/brave-villani-adc645

- 원본 `local-11`, branch `refs/heads/feat/dft-xc-grid-reuse`, HEAD `e931030b42dbc697cf98100bf60cfd5271838839`
- 존재: True; status 항목: 1

```text
?? build/
```
## /Users/cheolhochoi/Documents/claude/openqp/.claude/worktrees/condescending-maxwell-3d664c

- 원본 `local-11`, branch `detached`, HEAD `d749e29deb37e258561d5b1fcd81504f96f195fd`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp/.claude/worktrees/dft-xc-upstream

- 원본 `local-11`, branch `refs/heads/feat/dft-xc-grid-reuse-upstream`, HEAD `5e106e567fe16b050294408d5a27fdc0577822d1`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp/.claude/worktrees/gifted-montalcini-51ed09

- 원본 `local-11`, branch `refs/heads/chore/trim-example-reference-jsons`, HEAD `ed86fa809f0d46852d184a10c61e5bf64c425001`
- 존재: True; status 항목: 244

```text
M  .github/workflows/CI.yml
M  CMakeLists.txt
M  cmake/FindDDX.cmake
MM examples/DFT/H2O_RHF-DFT_ENERGY.json
MM examples/DFT/H2O_RHF-DFT_GRADIENT.json
MM examples/DFT/H2O_RHF-DFT_OPTIMIZE.json
MM examples/DFT/H2O_ROHF-DFT_ENERGY.json
MM examples/DFT/H2O_ROHF-DFT_GRADIENT.json
M  examples/DFT/H2O_UHF-DFT_ENERGY.inp
MM examples/DFT/H2O_UHF-DFT_ENERGY.json
M  examples/DFT/H2O_UHF-DFT_GRADIENT.inp
MM examples/DFT/H2O_UHF-DFT_GRADIENT.json
MM examples/ECP/C2H4_BHHLYP-MRSFTDDFT_Energy.json
MM examples/ECP/C2H4_BHHLYP-MRSFTDDFT_Grad.json
MM examples/ECP/Custom_Basis-Set/C2H4_Custom_Basis-Set-MRSFTDDFT_Energy.json
MM examples/ECP/HBr_BHHLYP-MRSFTDDFT_ENERGY.json
MM examples/ECP/HBr_RHF-DFT_ENERGY.json
MM examples/ECP/HBr_RHF-DFT_GRADIENT.json
MM examples/ECP/Mg_RHF-DFT_ENERGY.json
MM examples/ECP/NaCl_BHHLYP-MRSFTDDFT_ENERGY.json
MM examples/ECP/NaCl_RHF-DFT_GRADIENT.json
MM examples/HESS/H2O_BHHLYP-MRSFTDDFT_NUM_HESS.json
MM examples/HESS/H2O_RHF-DFT_NUM_HESS.json
MM examples/HF/H2O_RHF-HF_ENERGY.json
MM examples/HF/H2O_RHF-HF_GRADIENT.json
MM examples/HF/H2O_ROHF-HF_ENERGY.json
MM examples/HF/H2O_ROHF-HF_GRADIENT.json
M  examples/HF/H2O_UHF-HF_ENERGY.inp
MM examples/HF/H2O_UHF-HF_ENERGY.json
M  examples/HF/H2O_UHF-HF_GRADIENT.inp
MM examples/HF/H2O_UHF-HF_GRADIENT.json
MM examples/MRSF-TDDFT/C2H4_BHHLYP-MRSFTDDFT_MECI.json
MM examples/MRSF-TDDFT/C2H4_BHHLYP-MRSFTDDFT_MEP.json
MM examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.json
MM examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.json
MM examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_OPTIMIZE.json
MM examples/NMR/CH3_ROHF-DFT-PBE0-GIAO-NMR.json
MM examples/NMR/CH3_ROHF-GIAO-NMR.json
MM examples/NMR/CH3_UHF-DFT-PBE0-GIAO-NMR.json
MM examples/NMR/CH3_UHF-GIAO-NMR.json
MM examples/NMR/H2O_RHF-DFT-PBE0-GIAO-NMR.json
MM examples/NMR/H2O_RHF-GIAO-NMR.json
MM examples/NMR/H2O_RHF-NMR.json
MM examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECI.json
MM examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECI_GEOMETRIC.json
MM examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MECP_GEOMETRIC.json
MM examples/OPT/C2H4_BHHLYP-MRSFTDDFT_MEP.json
MM examples/OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.json
MM examples/OPT/H2O_BHHLYP-MRSFTDDFT_OPTIMIZE.json
MM examples/OPT/H2O_RHF-DFT_OPTIMIZE.json
MM examples/OPT/H2O_RHF-DFT_OPTIMIZE_GEOMETRIC.json
MM examples/OPT/HCN_BHHLYP-MRSFTDDFT_TS_GEOMETRIC.json
MM examples/OPT/HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.json
MM examples/OPT/HCN_RHF-DFT_IRC_GEOMETRIC.json
MM examples/OPT/HCN_RHF-DFT_TS_GEOMETRIC.json
M  examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.inp
MM examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.json
M  examples/PCM/OH_ROHF-HF_DDPCM_ENERGY.inp
MM examples/PCM/OH_ROHF-HF_DDPCM_ENERGY.json
MM examples/SCF/h2o_rhf_6-31g_pbe_adiis.json
MM examples/SCF/h2o_rhf_6-31g_pbe_basic.json
MM examples/SCF/h2o_rhf_6-31g_pbe_diis-reset.json
MM examples/SCF/h2o_rhf_6-31g_pbe_ediis.json
MM examples/SCF/h2o_rhf_6-31g_pbe_mom.json
MM examples/SCF/h2o_rhf_6-31g_pbe_pfon.json
MM examples/SCF/h2o_rhf_6-31g_pbe_soscf.json
MM examples/SCF/h2o_rhf_6-31g_pbe_vdiis.json
MM examples/SCF/h2o_rohf_6-31g_pbe_adiis.json
MM examples/SCF/h2o_rohf_6-31g_pbe_basic.json
MM examples/SCF/h2o_rohf_6-31g_pbe_diis-reset.json
MM examples/SCF/h2o_rohf_6-31g_pbe_ediis.json
MM examples/SCF/h2o_rohf_6-31g_pbe_mom.json
MM examples/SCF/h2o_rohf_6-31g_pbe_pfon.json
MM examples/SCF/h2o_rohf_6-31g_pbe_soscf.json
MM examples/SCF/h2o_rohf_6-31g_pbe_vshift.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_adiis.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_basic.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_diis-reset.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_ediis.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_mom.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_pfon.json
MM examples/SCF/h2o_uhf-s_6-31g_pbe_vdiis.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_adiis.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_basic.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_diis-reset.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_ediis.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_mom.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_pfon.json
MM examples/SCF/h2o_uhf-t_6-31g_pbe_vdiis.json
MM examples/SCF/h2o_uhf_6-31g_pbe_soscf.json
MM examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_ENERGY.json
MM examples/SF-TDDFT/H2O_BHHLYP-SFTDDFT_GRADIENT.json
MM examples/SOC/CH3Br-BHHLYP-SOC.json
MM examples/SOC/H2O_BHHLYP_SOC.json
MM examples/TDDFT/H2O_B3LYP5-TDDFT_ENERGY.json
MM examples/TDDFT/H2O_B3LYP5-TDDFT_GRADIENT.json
MM examples/TDHF/H2O_TDHF_ENERGY.json
MM examples/TDHF/H2O_TDHF_GRADIENT.json
MM examples/TRAH/H2O_RHF-DFT_GRAD_TRAH.json
MM examples/TRAH/H2O_ROHF-DFT_GRAD_TRAH.json
MM examples/TRAH/H2O_UHF-DFT_GRAD_TRAH.json
MM examples/TRAH/h2o_rohf_mrsf-t_6-31g_prop.json
MM examples/UMRSF-TDDFT/C4H6_BHHLYP_UMRSFTDDFT_ENERGY.json
MM examples/XAS/MRSF/HCN_MRSF.json
MM examples/XAS/delta-CHP-MRSF/HCN_CHP-MRSF_v2.json
MM examples/XAS/delta-CHP-MRSF/HCN_DFT.json
MM examples/other/h2o-2_rhf_cc-pvtz_hf.json
MM examples/other/h2o-2_rohf_cc-pvtz_hf.json
MM examples/other/h2o-2_uhf-s_cc-pvtz_hf.json
MM examples/other/h2o_nacme_rohf_mrsf-s_6-31g_bhhlyp.json
MM examples/other/h2o_rhf_3-21g_minao.json
MM examples/other/h2o_rhf_3-21g_modhuckel.json
MM examples/other/h2o_rhf_3-21g_sap.json
MM examples/other/h2o_rhf_6-31g_b3lypv5.json
MM examples/other/h2o_rhf_6-31g_bhhlyp.json
MM examples/other/h2o_rhf_6-31g_cam-b3lyp.json
MM examples/other/h2o_rhf_6-31g_hf.json
MM examples/other/h2o_rhf_6-31g_m06-2x.json
MM examples/other/h2o_rhf_6-31g_pbe.json
MM examples/other/h2o_rhf_6-31g_slater.json
MM examples/other/h2o_rhf_cc-pvtz_b3lypv5.json
MM examples/other/h2o_rhf_rpa-s_6-31g_b3lypv5.json
MM examples/other/h2o_rhf_rpa-s_6-31g_bhhlyp.json
MM examples/other/h2o_rhf_rpa-s_6-31g_cam-b3lyp.json
MM examples/other/h2o_rhf_rpa-s_6-31g_hf.json
MM examples/other/h2o_rhf_rpa-s_6-31g_m06-2x.json
MM examples/other/h2o_rhf_rpa-s_6-31g_pbe.json
MM examples/other/h2o_rhf_rpa-s_6-31g_slater.json
MM examples/other/h2o_rhf_rpa-t_6-31g_b3lypv5.json
MM examples/other/h2o_rhf_rpa-t_6-31g_bhhlyp.json
MM examples/other/h2o_rhf_rpa-t_6-31g_cam-b3lyp.json
MM examples/other/h2o_rhf_rpa-t_6-31g_hf.json
MM examples/other/h2o_rhf_rpa-t_6-31g_m06-2x.json
MM examples/other/h2o_rhf_rpa-t_6-31g_pbe.json
MM examples/other/h2o_rhf_rpa-t_6-31g_slater.json
MM examples/other/h2o_rhf_tda-s_6-31g_b3lypv5.json
MM examples/other/h2o_rhf_tda-s_6-31g_bhhlyp.json
MM examples/other/h2o_rhf_tda-s_6-31g_cam-b3lyp.json
MM examples/other/h2o_rhf_tda-s_6-31g_hf.json
MM examples/other/h2o_rhf_tda-s_6-31g_m06-2x.json
MM examples/other/h2o_rhf_tda-s_6-31g_pbe.json
MM examples/other/h2o_rhf_tda-s_6-31g_slater.json
MM examples/other/h2o_rhf_tda-t_6-31g_b3lypv5.json
MM examples/other/h2o_rhf_tda-t_6-31g_bhhlyp.json
MM examples/other/h2o_rhf_tda-t_6-31g_cam-b3lyp.json
MM examples/other/h2o_rhf_tda-t_6-31g_hf.json
MM examples/other/h2o_rhf_tda-t_6-31g_m06-2x.json
MM examples/other/h2o_rhf_tda-t_6-31g_pbe.json
MM examples/other/h2o_rhf_tda-t_6-31g_slater.json
MM examples/other/h2o_rohf-dft_energy_init_scf.json
MM examples/other/h2o_rohf-dft_energy_init_scf_basis_library.json
MM examples/other/h2o_rohf_6-31g_b3lypv5.json
MM examples/other/h2o_rohf_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_6-31g_cam-b3lyp.json
MM examples/other/h2o_rohf_6-31g_hf.json
MM examples/other/h2o_rohf_6-31g_m06-2x.json
MM examples/other/h2o_rohf_6-31g_pbe.json
MM examples/other/h2o_rohf_6-31g_slater.json
MM examples/other/h2o_rohf_cc-pvtz_b3lypv5.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_b3lypv5.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_cam-b3lyp.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dt-bhhlyp.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dt-vee.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-aee.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-b3lyp.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-stg.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-tune.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-vaee.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-vee.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-xi.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_dtcam-xiv.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_hf.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_m06-2x.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_pbe.json
MM examples/other/h2o_rohf_mrsf-q_6-31g_slater.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_b3lypv5.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp-spc-coco.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp-spc-coov.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp-spc-ovov.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_cam-b3lyp.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dt-bhhlyp-spc.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dt-bhhlyp.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dt-vee.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-aee.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-b3lyp.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-stg.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-tune.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-vaee.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-vee.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-xi.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_dtcam-xiv.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_hf.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_m06-2x.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_pbe.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_prop.json
MM examples/other/h2o_rohf_mrsf-s_6-31g_slater.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_b3lypv5.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_cam-b3lyp.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_dt-bhhlyp-spc.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_dt-bhhlyp.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_dt-vee.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_dtcam-b3lyp.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_hf.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_m06-2x.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_pbe.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_slater.json
MM examples/other/h2o_rohf_mrsf-t_6-31g_stg1x.json
MM examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_sf_6-31g_b3lypv5.json
MM examples/other/h2o_rohf_sf_6-31g_bhhlyp.json
MM examples/other/h2o_rohf_sf_6-31g_cam-b3lyp.json
MM examples/other/h2o_rohf_sf_6-31g_dt-bhhlyp.json
MM examples/other/h2o_rohf_sf_6-31g_dtcam-b3lyp.json
MM examples/other/h2o_rohf_sf_6-31g_hf.json
MM examples/other/h2o_rohf_sf_6-31g_m06-2x.json
MM examples/other/h2o_rohf_sf_6-31g_pbe.json
MM examples/other/h2o_rohf_sf_6-31g_slater.json
MM examples/other/h2o_uhf-s_6-31g_b3lypv5.json
MM examples/other/h2o_uhf-s_6-31g_bhhlyp.json
MM examples/other/h2o_uhf-s_6-31g_cam-b3lyp.json
MM examples/other/h2o_uhf-s_6-31g_hf.json
MM examples/other/h2o_uhf-s_6-31g_m06-2x.json
MM examples/other/h2o_uhf-s_6-31g_pbe.json
MM examples/other/h2o_uhf-s_6-31g_slater.json
MM examples/other/h2o_uhf-s_cc-pvtz_b3lypv5.json
MM examples/other/h2o_uhf-t_6-31g_b3lypv5.json
M  examples/other/h2o_uhf-t_6-31g_bhhlyp.inp
MM examples/other/h2o_uhf-t_6-31g_bhhlyp.json
M  examples/other/h2o_uhf-t_6-31g_cam-b3lyp.inp
MM examples/other/h2o_uhf-t_6-31g_cam-b3lyp.json
M  examples/other/h2o_uhf-t_6-31g_hf.inp
MM examples/other/h2o_uhf-t_6-31g_hf.json
M  examples/other/h2o_uhf-t_6-31g_m06-2x.inp
MM examples/other/h2o_uhf-t_6-31g_m06-2x.json
MM examples/other/h2o_uhf-t_6-31g_pbe.json
MM examples/other/h2o_uhf-t_6-31g_slater.json
M  external/CMakeLists.txt
MM pyoqp/oqp/molecule/molecule.py
M  pyoqp/oqp/pyoqp.py
M  pyoqp/oqp/utils/oqp_tester.py
```
## /Users/cheolhochoi/Documents/claude/openqp/.claude/worktrees/jolly-bell-068a3f

- 원본 `local-11`, branch `detached`, HEAD `4d8afd6e122fe5435f263eb2cf5ffa7b41d9d3dc`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/openqp/.claude/worktrees/vigorous-elion-dd912c

- 원본 `local-11`, branch `refs/heads/fix/ddx-linux-link-and-ci-enable`, HEAD `ed6c0905e551d0db4ee2b9204b073b75cbc7f32f`
- 존재: True; status 항목: 3

```text
M  .github/workflows/CI.yml
M  CMakeLists.txt
M  cmake/FindDDX.cmake
```
## /Users/cheolhochoi/Documents/openqp-dftd4-fork

- 원본 `local-11`, branch `refs/heads/feat/native-dftd4`, HEAD `4a335a3230fe1e2794a42f28eb12d1cd767ac8be`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/openqp-dftd4-native

- 원본 `local-11`, branch `refs/heads/feat/dftd4-native`, HEAD `4d8afd6e122fe5435f263eb2cf5ffa7b41d9d3dc`
- 존재: True; status 항목: 7

```text
 M Dockerfile
 M external/CMakeLists.txt
 M include/oqp.h
 M pyoqp/oqp/__init__.py
 M pyoqp/oqp/library/single_point.py
 M source/CMakeLists.txt
?? source/dftd4_interface.F90
```
## /Users/cheolhochoi/Documents/openqp-pr367-fix

- 원본 `local-11`, branch `refs/heads/claude/ci-spin-cluster-20260821`, HEAD `83c48ad56b54698c9a7847b5631cd06d5cba35a9`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/openqp-uhf-ci-fix

- 원본 `local-11`, branch `refs/heads/fix/uhf-stability-and-ci-runtests`, HEAD `24c8c5c8c8a329ea64d959db3274bfdf1a861cd5`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/simplify-oqp-inputs-20260819

- 원본 `local-11`, branch `refs/heads/docs/simplify-oqp-inputs`, HEAD `ec5106648a6a955fd93cf9d8a2b40d3090cc7d3a`
- 존재: True; status 항목: 0
## /Volumes/External_Storage/claude/sessions/20260612_143000_pr199_ispher_review/openqp-impl

- 원본 `local-11`, branch `refs/heads/feat/rotaxis-direct-pure`, HEAD `20f504d410c2b63f8f1163f813fc5b245455c761`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Volumes/External_Storage/claude/sessions/20260612_143000_pr199_ispher_review/openqp-pr199

- 원본 `local-11`, branch `detached`, HEAD `fb7bd53ffca483ce5307cd3822c5b29f6ab82b35`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Volumes/External_Storage/claude/sessions/20260612_143000_pr199_ispher_review/openqp-ryscart

- 원본 `local-11`, branch `refs/heads/bench/rys-cart`, HEAD `91626817b92326a3ab1836af7e3e9b427015dea0`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Volumes/External_Storage/claude/sessions/20260617_154541_strip-dev-artifacts/openqp-cleanup

- 원본 `local-11`, branch `refs/heads/chore/strip-dev-artifacts`, HEAD `666aa200bc9f30e317fd2aca6942f91c521968ef`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Volumes/External_Storage/claude/sessions/20260628_perf_audit_levels/wt-171

- 원본 `local-11`, branch `refs/heads/fix/trah-rstctmo-handoff`, HEAD `11dd37f7f35988a0b624c34115d6085a1df50f61`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Volumes/External_Storage/claude/sessions/20260628_perf_audit_levels/wt-perf

- 원본 `local-11`, branch `refs/heads/feat/perf-levels`, HEAD `365391f66c3e6e32e4e08de6716803a73dcb50b2`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Volumes/OpenQP.DB/claude/sessions/20260819_081024_tutorials-oqp-format/openqp-upstream

- 원본 `local-11`, branch `detached`, HEAD `20ea3f4f252cea188dbff2385578ee39b84c264d`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_082115_readme-drop-tight-binding/wt-readme

- 원본 `local-11`, branch `refs/heads/docs/drop-tight-binding`, HEAD `c21f216b2e66e78dfe72528e03b8149cd9a08133`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-349

- 원본 `local-11`, branch `refs/heads/wt/pr349-conflict`, HEAD `1116c0568fdc0de0597efd2fa6ef1c62b82619cf`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-350

- 원본 `local-11`, branch `refs/heads/wt/pr350-conflict`, HEAD `11965fab7fa961848717a76fca6af31bc87e21db`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-351

- 원본 `local-11`, branch `refs/heads/wt/pr351`, HEAD `fa52fd31229cef143e3fe3af70e79e32c03d05de`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-354

- 원본 `local-11`, branch `refs/heads/wt/pr354`, HEAD `2ea620513530c277feeb0432bfb6d3a5f7c99d65`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-build

- 원본 `local-11`, branch `detached`, HEAD `0bd0050d03a9bc74cdba2d4761bfc5180e5d9e67`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-main

- 원본 `local-11`, branch `detached`, HEAD `8fa5a1c6f24249eacfee31a9e81487af3c5fd33f`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-readme

- 원본 `local-11`, branch `refs/heads/fix/readme-sa-casscf-analytic`, HEAD `c1f378d9786d2ee726ae63a548cd37b6e734f925`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-route

- 원본 `local-11`, branch `refs/heads/fix/scnevpt2-route-contract`, HEAD `5f7dc6468bf708507d74980362f31b54a2b27d56`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260819_pr347_353/wt-sim

- 원본 `local-11`, branch `detached`, HEAD `4e4c716f9bf8013e5dc6f2e42308cd461f5658dc`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260822_160000_pr368-ci-fix/wt-363

- 원본 `local-11`, branch `refs/heads/wt/pr363-ci-fix`, HEAD `cc222b55cc0685c5076f6d95ecce67be9a09a89e`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260822_160000_pr368-ci-fix/wt-368

- 원본 `local-11`, branch `refs/heads/wt/pr368-fix`, HEAD `26a2e0b185f1321b4b70398530669444e57c858b`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260919_002020_claude-codex-review-fix/wt-claudeyml

- 원본 `local-11`, branch `refs/heads/fix/claude-review-track-progress`, HEAD `6e007eb4ab70f2ddd4d2838dba37e9028eac95ff`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260919_002020_claude-codex-review-fix/wt-split

- 원본 `local-11`, branch `refs/heads/chore/split-devkit-and-layout-gate`, HEAD `7938b613ca66ff357c8fe00f692cb1b9a190e828`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/sessions/20260731_nac_audit/openqp_upstream

- 원본 `local-12`, branch `refs/heads/fix/fcheck-bounds-aborts`, HEAD `32a8e1214890981bb064c11d7afdc58bcd1836b1`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/claude/sessions/20260731_nac_audit/repo

- 원본 `local-13`, branch `refs/heads/fix/fcheck-bounds-guess-dftlib`, HEAD `3245c0e56a5262946b77334e889466c28b027245`
- 존재: True; status 항목: 10

```text
?? .venv/
?? .venv_main_g13/
?? build/
?? build_main_dbg_g13/
?? build_pip/
?? build_pip_debug/
?? build_pip_debug_g13/
?? install_main_dbg_g13/
?? pyoqp/PyOpenQP.egg-info/
?? pyoqp/build/
```
## /Users/cheolhochoi/Documents/claude/sessions/20260731_nac_audit/repo_nac

- 원본 `local-13`, branch `refs/heads/nac-lagrangian`, HEAD `81d41e6200fba0349680cece1c0983eb441094b1`
- 존재: True; status 항목: 1

```text
?? tmp/
```
## /Users/cheolhochoi/Documents/openqp-private

- 원본 `local-14`, branch `refs/heads/feat/mrsf-ensemble-reference`, HEAD `2804a0ec3e5de40eba5d183c3304ac6cc9524654`
- 존재: True; status 항목: 8

```text
 M include/oqp.h
 M pyoqp/oqp/library/single_point.py
 M pyoqp/oqp/utils/mrsf_reference.py
 M source/modules/tdhf_mrsf_energy.F90
 M source/tagarray_driver.F90
 M source/tdhf_mrsf_lib.F90
?? CHECKPOINT_MRSF_DETERMINANT_UNION.md
?? tools/_mrsf_response_smoke/
```
## /private/tmp/openqp-pr199-check

- 원본 `local-14`, branch `detached`, HEAD `16e356791566b828704d619f669698ee3dc2ebcb`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /private/tmp/openqp-pr201-check

- 원본 `local-14`, branch `detached`, HEAD `6bb207ccee746143550e70e535cc8d56021f8f1c`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /private/tmp/openqp-pr202-check

- 원본 `local-14`, branch `detached`, HEAD `33b706dd9d0f06ec898f827b9ce4a3d43dc1eec1`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Users/cheolhochoi/.codex/worktrees/0be9/openqp-private

- 원본 `local-14`, branch `refs/heads/claude/ptc-tenno-umu7ko`, HEAD `9041f1c7f14d28a5566728087c289d5d379e4c89`
- 존재: True; status 항목: 2

```text
 M tests/ptc_mrsf/prototype/pes_tenno.dat
 M tests/ptc_mrsf/prototype/tc_h2_pes_tenno.F90
```
## /Users/cheolhochoi/.codex/worktrees/8dee/openqp-private

- 원본 `local-14`, branch `refs/heads/codex/review-qmrsf-dual-pathways`, HEAD `2c26faabec16bdef235bb979db95448e3985e9ed`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/.codex/worktrees/b1c8/openqp-private

- 원본 `local-14`, branch `refs/heads/codex/macos-blas-lp64-bench-upstream`, HEAD `d0c91ad56ef123aa32a81cdbdd1c6f2a4f6662d5`
- 존재: True; status 항목: 12

```text
?? _build_local_accelerate_lp64/
?? _build_local_openblas_ilp64/
?? _install_local_accelerate_lp64/
?? _install_local_openblas_ilp64/
?? _install_local_openblas_lp64/
?? _probe_ci_key_fix/
?? _probe_ilp64_netlib/
?? _probe_ilp64_openblas/
?? _probe_ilp64_openblas_explicit/
?? _probe_lp64_accelerate/
?? _probe_lp64_openblas/
?? blas_matrix_results/
```
## /Users/cheolhochoi/.codex/worktrees/b6b8/openqp-private

- 원본 `local-14`, branch `refs/heads/codex/pr199-ispher-fixes`, HEAD `fd7e8fb412d7e3c508159e72ce097da8bb76eea4`
- 존재: True; status 항목: 14

```text
 M docs/plans/2026-06-09-ispher-spherical-harmonics.md
 M pyoqp/oqp/library/set_basis.py
 M pyoqp/oqp/library/single_point.py
 M pyoqp/oqp/molecule/oqpdata.py
 M source/integrals/int2.F90
 M source/integrals/int_rotaxis.F90
 M source/integrals/int_rys.F90
 M tests/test_ispher_runtime_keyword.py
?? pyoqp/oqp/utils/basis_convention.py
?? source/integrals/int2_rys_pure_generated.F90
?? source/integrals/int_rotaxis_pure_generated.F90
?? tools/bench_rys_pure_paths.py
?? tools/generate_int2_rys_pure_kernels.py
?? tools/generate_int_rotaxis_pure_d.py
```
## /Users/cheolhochoi/.codex/worktrees/ce2c/openqp-private

- 원본 `local-14`, branch `detached`, HEAD `d5bd9724e420bf1515c223547dc6ee8bbebee370`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/.codex/worktrees/mrsf-h2

- 원본 `local-14`, branch `refs/heads/fix/mrsf-h2-zero-closed`, HEAD `8733067effb5912dc1c20703abedc5b8c9bdf241`
- 존재: True; status 항목: 5

```text
 M tests/mrsf_h2/h2_mrsfcis_uhf.inp
?? tests/mrsf_h2/h2_rohf_stab.inp
?? tests/mrsf_h2/h2_rohf_stabauto.inp
?? tests/mrsf_h2/h2_rohf_swap.inp
?? tests/mrsf_h2/h2_rohf_trah.inp
```
## /Users/cheolhochoi/Documents/Codex/2026-08-12/qmrsf-orbital-rotation-14-6/work/openqp-c4dd-readonly

- 원본 `local-14`, branch `detached`, HEAD `c4dd34f895a26ec1b4377960bf32cad40dd49da2`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-26/openqp-fort6-cleanup-20260826

- 원본 `local-14`, branch `refs/heads/codex/fix-fort6-cleanup-20260826`, HEAD `293308b7c30aa7663f4c4b049eb6efe14f467a0e`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-26/openqp-ssc-h2co-benchmark

- 원본 `local-14`, branch `refs/heads/codex/ssc-h2co-benchmark-20260826`, HEAD `0b070671a4a47b34432d94928b4d4360871b2c7b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-26/openqp-x2c-soc-science

- 원본 `local-14`, branch `refs/heads/codex/x2c-soc-science-20260826`, HEAD `e51cf8a016bf15d4412c12b183b6a5f49649f2a2`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-26/openqp-x2c-ssc-pes

- 원본 `local-14`, branch `refs/heads/codex/x2c-ssc-pes-20260826`, HEAD `40d151a7853bf1163a97bbd7c7cd30fabecf0e55`
- 존재: True; status 항목: 3

```text
?? pyoqp/oqp/library/__pycache__/
?? tests/__pycache__/
?? tmp/
```
## /Users/cheolhochoi/Documents/Codex/2026-08-27/openqp-ssc-zfs-pes-calc-20260827

- 원본 `local-14`, branch `refs/heads/codex/ssc-zfs-pes-calc-20260827`, HEAD `40d151a7853bf1163a97bbd7c7cd30fabecf0e55`
- 존재: True; status 항목: 2

```text
?? campaigns/
?? tmp/
```
## /Users/cheolhochoi/Documents/Codex/2026-08-31/openqp-nac-student-handoff-20260831

- 원본 `local-14`, branch `refs/heads/codex/student-handoff-20260831`, HEAD `f78eb044268ac00aa7416b9a8665449d94e10a72`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-09-01/qmrsf-be-c-operator-audit

- 원본 `local-14`, branch `refs/heads/codex/qmrsf-be-c-operator-audit-20260901`, HEAD `f7572f7667ddfb902db7b8bb220adaf2d9735b34`
- 존재: True; status 항목: 10

```text
?? tools/qmrsf_operator_audit/CP4_QMRSF_SPEC.md
?? tools/qmrsf_operator_audit/SPECTRAL_PROJECTOR_COVARIANT_QMRSF.md
?? tools/qmrsf_operator_audit/audit_pair_channel_topology.py
?? tools/qmrsf_operator_audit/develop_covariant_pair_kernel.py
?? tools/qmrsf_operator_audit/develop_spectral_projector_qmrsf.py
?? tools/qmrsf_operator_audit/output/covariant_pair_kernel_development.json
?? tools/qmrsf_operator_audit/output/fixed_manifold_cp4_verification.json
?? tools/qmrsf_operator_audit/output/pair_channel_topology_audit.json
?? tools/qmrsf_operator_audit/output/spectral_projector_qmrsf.json
?? tools/qmrsf_operator_audit/verify_fixed_manifold_cp4.py
```
## /Users/cheolhochoi/Documents/Codex/2026-09-06/tda-nmr-origin-regression

- 원본 `local-14`, branch `refs/heads/codex/tda-nmr-origin-regression-20260906`, HEAD `98930fd73cf1530f39ad59cbfa5bfe6491198427`
- 존재: True; status 항목: 1

```text
?? runs/
```
## /Users/cheolhochoi/Documents/Codex/2026-09-07/mrsf-nmr-build-final

- 원본 `local-14`, branch `detached`, HEAD `c59c74700d1105e86c8a1cfdd59bd1c29241d23e`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-09-07/mrsf-nmr-csf-response

- 원본 `local-14`, branch `refs/heads/codex/mrsf-nmr-csf-response-20260907`, HEAD `fff81acc29a90568b4962f5c0b5e12bf4c80547d`
- 존재: True; status 항목: 1

```text
?? runs/
```
## /Users/cheolhochoi/Documents/Codex/2026-09-08/nmr-integral-response-validation

- 원본 `local-14`, branch `refs/heads/codex/nmr-integral-response-validation-20260908`, HEAD `edd7301654eb7700c2792640d58dcf0193d5a8db`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-09-08/nmr-metric-validation

- 원본 `local-14`, branch `refs/heads/codex/nmr-metric-validation-20260908`, HEAD `a7c6d07b66e97a12a67e85bbb54a474c728e8b23`
- 존재: True; status 항목: 1

```text
?? tests/test_mrsf_nmr_metric_covariance.py
```
## /Users/cheolhochoi/Documents/openqp-digestion

- 원본 `local-14`, branch `refs/heads/perf/mrsf-fock-digestion`, HEAD `fc868854932de3c5aaf6b78ef99b2c370308441d`
- 존재: True; status 항목: 1

```text
?? build-ddx/
```
## /Users/cheolhochoi/Documents/openqp-dk-dft

- 원본 `local-14`, branch `refs/heads/codex/qmrsf-dk-original2-validation`, HEAD `3e5c61b3be96aac8f1b2cedea04f35ecf1c43e85`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/openqp-dk-fullspace

- 원본 `local-14`, branch `refs/heads/feat/dk-fullspace-singles`, HEAD `bc75f319d8afa20f5330112cad22ebd218834078`
- 존재: True; status 항목: 1

```text
 M source/modules/tdhf_qmrsf_dk.F90
```
## /Users/cheolhochoi/Documents/openqp-dk-worktree

- 원본 `local-14`, branch `refs/heads/feat/qmrsf-dk-live`, HEAD `70c9825af713ff902023c16e30133a2662bee0f5`
- 존재: True; status 항목: 2

```text
?? tools/qmrsf_pathways_proto/stageB/h4_quintet_dk.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/qmrsf_dk_full_live.dat
```
## /Users/cheolhochoi/Documents/openqp-doublet

- 원본 `local-14`, branch `refs/heads/feat/doublet-quartet-mrsf`, HEAD `3ee72f82fc5276fe4cba15d67c2a048166825e36`
- 존재: True; status 항목: 3

```text
 M docs/doublet_quartet/aodirect_two_density_backfock.py
 M source/modules/qmrsf_doublet_assemble.F90
 M source/modules/tdhf_qmrsf_doublet.F90
```
## /Users/cheolhochoi/Documents/openqp-mecp-fix

- 원본 `local-14`, branch `refs/heads/fix/mecp-converging-objectives`, HEAD `b019d4caa92100aafbe6d5b84be63ba81b9a54f0`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/openqp-private-qmrsf-pathways

- 원본 `local-14`, branch `refs/heads/feat/qmrsf-dual-pathways`, HEAD `8025f9167c318700f76a7c41aa3a23b6dfba04b1`
- 존재: True; status 항목: 132

```text
 M tools/qmrsf_pathways_proto/fortran/qmrsf_icpt2_full_live.dat
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_631g_dk.inp
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_631g_dk.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_631g_dkg.inp
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_631g_dkg.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_631g_icpt2.inp
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_631g_icpt2.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_sto3g_dk.inp
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_sto3g_dk.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_sto3g_dkg.inp
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_sto3g_dkg.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_sto3g_icpt2.inp
?? tools/qmrsf_pathways_proto/stageB/bench_CBD_sto3g_icpt2.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_H4_631g_dk.inp
?? tools/qmrsf_pathways_proto/stageB/bench_H4_631g_dk.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_H4_631g_dkg.inp
?? tools/qmrsf_pathways_proto/stageB/bench_H4_631g_dkg.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_H4_631g_icpt2.inp
?? tools/qmrsf_pathways_proto/stageB/bench_H4_631g_icpt2.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_H4_sto3g_dk.inp
?? tools/qmrsf_pathways_proto/stageB/bench_H4_sto3g_dk.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_H4_sto3g_dkg.inp
?? tools/qmrsf_pathways_proto/stageB/bench_H4_sto3g_dkg.qmrsf_dk.json
?? tools/qmrsf_pathways_proto/stageB/bench_H4_sto3g_icpt2.inp
?? tools/qmrsf_pathways_proto/stageB/bench_H4_sto3g_icpt2.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D2h_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D2h_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D2h_ccpvtz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D2h_ccpvtz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D4h_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D4h_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D4h_ccpvtz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_loos_D4h_ccpvtz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scal_CBD_631g.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scal_CBD_631g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scal_CBD_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scal_CBD_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scal_CBD_sto3g.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scal_CBD_sto3g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scal_H4_631g.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scal_H4_631g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scal_H4_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scal_H4_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scal_H4_sto3g.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scal_H4_sto3g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_00.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_00.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_01.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_01.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_02.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_02.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_03.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_03.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_04.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_04.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_05.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_05.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_06.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_06.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_07.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_07.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_08.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_08.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_09.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_09.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_10.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_10.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_11.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_11.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_12.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_12.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_13.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanLIN_13.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_00.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_00.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_01.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_01.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_02.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_02.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_03.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_03.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_04.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_04.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_05.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_05.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_06.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_06.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_07.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_07.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_08.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_08.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_09.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_09.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_10.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_10.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_11.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_11.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_12.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_12.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_13.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_13.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_14.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_14.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_15.inp
?? tools/qmrsf_pathways_proto/stageB/bench_scanSR_15.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_tf_CBDrect_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_tf_CBDrect_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_tf_CBDsquare_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_tf_CBDsquare_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_tf_H4_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_tf_H4_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_tf_TME_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_tf_TME_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/bench_tf_TMM_ccpvdz.inp
?? tools/qmrsf_pathways_proto/stageB/bench_tf_TMM_ccpvdz.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/h4_quintet_icpt2.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/h4_quintet_icpt2_631g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/h6_quintet_icpt2.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/oracle_ao.npz
?? tools/qmrsf_pathways_proto/stageB/poly_TME_631g.inp
?? tools/qmrsf_pathways_proto/stageB/poly_TME_631g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/poly_TME_sto3g.inp
?? tools/qmrsf_pathways_proto/stageB/poly_TME_sto3g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/poly_TMM_631g.inp
?? tools/qmrsf_pathways_proto/stageB/poly_TMM_631g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/poly_TMM_sto3g.inp
?? tools/qmrsf_pathways_proto/stageB/poly_TMM_sto3g.qmrsf.json
?? tools/qmrsf_pathways_proto/stageB/qmrsf_cact_live.dat
?? tools/qmrsf_pathways_proto/stageB/qmrsf_cfull_live.dat
?? tools/qmrsf_pathways_proto/stageB/qmrsf_dk_full_live.dat
?? tools/qmrsf_pathways_proto/stageB/qmrsf_icpt2_full_live.dat
?? tools/qmrsf_pathways_proto/stageB/qmrsf_icpt2_live.dat
```
## /Users/cheolhochoi/Documents/openqp-sqp

- 원본 `local-14`, branch `refs/heads/feat/meci-sqp`, HEAD `aef07199a6e00047e83c09a78a070dad356a99f4`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/openqp-ssc-fix

- 원본 `local-14`, branch `refs/heads/codex/fix-ssc-reference-contraction`, HEAD `0b070671a4a47b34432d94928b4d4360871b2c7b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/analnac-merged-regression-20260919

- 원본 `local-14`, branch `refs/heads/codex/analnac-merged-regression-20260919`, HEAD `7b3ba5190f8b3059a1f61d54fbde31fbf14b08ec`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/mrsf-response-diagnostic-20260914-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/mrsf-response-diagnostic-20260914-01a07f50`, HEAD `fcdf2091bc71c440159f0e29ac9f31570fda8e0b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/namd-local-continuation-20260914-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/namd-local-continuation-20260914-01a07f50`, HEAD `e3cdf5268a77ac045ebb1d7f9f97abc887c6313a`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/namd-selected-nac-20260908-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/namd-selected-nac-20260908-01a07f50`, HEAD `406f5a0a6576a2eecb9b53171e74150aededf263`
- 존재: True; status 항목: 8

```text
 M include/oqp.h
 M pyoqp/oqp/library/nac_analytic.py
 M pyoqp/oqp/library/namd.py
 M source/modules/mrsf_nac_driver.F90
 M source/modules/mrsf_nac_metric_data.F90
 M tests/test_mrsf_nac_fortran_driver.py
 M tests/test_mrsf_nac_tagarray_lifetimes.py
 M tests/test_namd_analytic_directional.py
```
## /Users/cheolhochoi/openqp-worktrees/static-analytic-nac-input-20260910-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/static-analytic-nac-input-20260910-01a07f50`, HEAD `1f508e4f9f2ebe3880cfd18a2e6451ff4ac9c8b7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/thymine-five-20260919-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/thymine-five-20260919-01a07f50`, HEAD `eb8877debe501a6d43c501aea0bab86467a43668`
- 존재: True; status 항목: 2

```text
?? campaign-thymine-five/plot_populations.py
?? campaign-thymine-five/recover_tracking.py
```
## /Users/cheolhochoi/openqp-worktrees/thymine-ic0103-scf-20260920

- 원본 `local-14`, branch `refs/heads/codex/thymine-ic0103-scf-20260920`, HEAD `b6f3ebf5cdfab01c784ba13512fd4357a122d0f3`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/thymine-namd-20260917-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/thymine-namd-20260917-01a07f50`, HEAD `3a9157eb70ef05d42ce6a726f3af50b5daa92885`
- 존재: True; status 항목: 1

```text
?? campaigns/
```
## /Users/cheolhochoi/openqp-worktrees/trah-strict-dt-20260914-01a07f50

- 원본 `local-14`, branch `refs/heads/codex/trah-strict-dt-20260914-01a07f50`, HEAD `4ec64a5d3ffd9aca427fe9f187842b6dff633cf3`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260905_mrsf_nmr_root_cause/openqp-tda-nmr

- 원본 `local-14`, branch `refs/heads/claude/tda-nmr-cgo-20260905`, HEAD `194a0504befcef17c64e7ef8f13dde3f463d0c61`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260905_uracil_namd_nac_check/openqp-analnac

- 원본 `local-14`, branch `refs/heads/claude/namd-analytic-nac-handoff-20260907`, HEAD `406f5a0a6576a2eecb9b53171e74150aededf263`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260905_uracil_namd_nac_check/openqp-selnac

- 원본 `local-14`, branch `refs/heads/claude/namd-selected-nac-20260908`, HEAD `1f508e4f9f2ebe3880cfd18a2e6451ff4ac9c8b7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp

- 원본 `local-15`, branch `refs/heads/codex/dftb-s0-continuity`, HEAD `57872b5f6248bb139aab3202b138361652ed3324`
- 존재: True; status 항목: 20

```text
 M pyoqp/oqp/library/namd.py
 M pyoqp/oqp/library/openqp_dftb.py
 M pyoqp/oqp/library/set_basis.py
 M pyoqp/oqp/molecule/molecule.py
 M pyoqp/oqp/molecule/oqpdata.py
 M source/bragg_slater.F90
 M source/dftlib/libxc.F90
 M source/integrals/grd2_rys.F90
 M source/integrals/mod_1e_primitives.F90
 M source/io/messages.F90
 M source/modules/apply_basis.F90
 M source/modules/dk_scalar.F90
 M source/modules/namd.F90
 M source/modules/qmmm.F90
 M source/modules/soc_mrsf.F90
 M source/modules/tdhf_mrsf_ekt.F90
 M tests/test_mrsf_ekt_scaffold.py
 M tests/test_tda_gradient_zvector_hessian.py
 M tools/libxc/gen_all.py
?? legal/
```
## /private/tmp/claude-501/-Users-cheolhochoi-Documents-claude/1960767a-ebce-431c-b7c6-3241a3675e79/scratchpad/ea_cluster/wt-upstream

- 원본 `local-15`, branch `detached`, HEAD `8324dac17c7eb87853d71635f37fb70e77da9829`
- 존재: False; status 항목: 0
- Git 표시: gitdir file points to non-existent location
## /Users/cheolhochoi/.codex/worktrees/openqp-readme-web-link

- 원본 `local-15`, branch `refs/heads/codex/readme-web-link`, HEAD `3c8d1eeb0785a3b58311f636a912a69d7ce09257`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Applications/openqp-s3r-seam-continuity-20260814-019ff5e3

- 원본 `local-15`, branch `refs/heads/codex/fix-s3r-seam-continuity-20260814-019ff5e3`, HEAD `0bacd2d647947776fa1999c7ff97ccf2cd9e88a5`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/benzene_openqp_symmetry_check/upstream-main-worktree

- 원본 `local-15`, branch `detached`, HEAD `34aece888ef69e0004e460a6f9c57dcfe5501634`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-ci-concurrency

- 원본 `local-15`, branch `refs/heads/claude/ci-concurrency-group2`, HEAD `32da9b467affc08b7ced6ec0b715053458323ca8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-dk-seam

- 원본 `local-15`, branch `refs/heads/feat/qmrsf-dk-covariant-seam`, HEAD `c4dd34f895a26ec1b4377960bf32cad40dd49da2`
- 존재: True; status 항목: 1

```text
?? tests/test_qmrsf_dk_doublet_invariance.py
```
## /Users/cheolhochoi/clone/openqp-main-ref

- 원본 `local-15`, branch `detached`, HEAD `fde78570274049d60e1cc43421fa4c01a58bbbf7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-mrsf-window-rebase

- 원본 `local-15`, branch `refs/heads/claude/mrsf-davidson-window-rebased-20260822`, HEAD `a31ae4cd4c576eba58d05e6f392eb7d3c3010ad3`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-perf-grad

- 원본 `local-15`, branch `refs/heads/perf/hf_dft_gradient`, HEAD `3d6cf9b0f29162d408aefd5a4b8e59b9b4188b1d`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-perf-int

- 원본 `local-15`, branch `refs/heads/perf/omp-input`, HEAD `df2a2e7f4119f3652db39dcd68b3839a136837ce`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-pr363-verify

- 원본 `local-15`, branch `refs/heads/claude/pr363-verify-20260822`, HEAD `7bc51dbd171869e0f2f8de9f3233ab32764947a6`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-scf

- 원본 `local-15`, branch `refs/heads/feat/dft-grid-default`, HEAD `c60285edcb0d276a659be35173cc4883bc9c8f5a`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/openqp-symmetry

- 원본 `local-15`, branch `refs/heads/feat/molecular-symmetry`, HEAD `6103c85c70e3d10af9bda20f1773b348f08616b8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Application 작업/openqp-v1.2.1-release

- 원본 `local-15`, branch `refs/heads/codex/release-v1.2.1`, HEAD `5250749082536345ff07e1d193979cd98c0055b9`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Code Checking/openqp-abi1-compat

- 원본 `local-15`, branch `refs/heads/agent/openqp-dftb-abi-fallback`, HEAD `6820c4b314e3c4deed6e55d855fd4eca05e036fc`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Code Checking/openqp-dftb-ui

- 원본 `local-15`, branch `refs/heads/agent/dftb-input-logging`, HEAD `346084274e503184abce31db9a3e7283c5d32d02`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-14/qmrsf-sr-code-integration-20260814/work/openqp-qmrsf-sr-canonical-019fff32

- 원본 `local-15`, branch `refs/heads/codex/qmrsf-sr-canonical-roks-20260814-019fff32`, HEAD `daf0f444877228d37f00cb42dc8f6aef8f1d32c0`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-15/openqp-casscf-analytic-gradient-claude

- 원본 `local-15`, branch `refs/heads/claude/casscf-analytic-gradient-20260815`, HEAD `20d235721e66ecb555e37bf8b8dae837fa00db22`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-15/openqp-casscf-numgrad-readme

- 원본 `local-15`, branch `refs/heads/codex/casscf-numgrad-readme-20260815`, HEAD `ad2fe53a0d2c662852287b5cd1d649787072b9f1`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-log-consistency-authors/work/openqp-log-consistency-authors

- 원본 `local-15`, branch `refs/heads/codex/openqp-log-consistency-authors-20260816`, HEAD `0baea2f17d14ba2a8a26b52f2e0464a249d5bcff`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-mp2-analytic-gradient-codex

- 원본 `local-15`, branch `refs/heads/codex/mp2-analytic-gradient-20260816`, HEAD `e7082ac056bbafba4cdc85f9bfb6f898a333695b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-mp2-analytic-gradient-specialist

- 원본 `local-15`, branch `refs/heads/codex/mp2-analytic-gradient-specialist-20260816`, HEAD `8b7ad53e5d7ec8a6f94af7ee34c9490d5b34dec5`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-pr344-review-maintenance

- 원본 `local-15`, branch `refs/heads/codex/casscf-analytic-gradient-review-v2-20260816`, HEAD `3bdfb60dbd9c94ffdb0d8c487d8fcdd9a2ebcc78`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-pr345-ci-review-maintenance

- 원본 `local-15`, branch `refs/heads/codex/pr345-ci-review-maintenance-20260816`, HEAD `2a0158208d0767e78d97c59d1a627a14a1df7965`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-16/openqp-sa-casscf-zvector-claude

- 원본 `local-15`, branch `refs/heads/claude/sa-casscf-zvector-20260816`, HEAD `95274875095c894465d16c95e049e7f4f2c6604f`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-17/openqp-caspt2-analytic-gradient-claude

- 원본 `local-15`, branch `refs/heads/claude/caspt2-analytic-gradient-20260817`, HEAD `1b6872969208106f854309e96719aea0910f2d4e`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-17/openqp-nevpt2-analytic-gradient-claude

- 원본 `local-15`, branch `refs/heads/claude/nevpt2-analytic-gradient-20260817`, HEAD `45aac8f4611f1a3968ffa0b3d52695024acd18e8`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-17/openqp-pr347-oqup-pronunciation

- 원본 `local-15`, branch `refs/heads/codex/pr347-oqup-pronunciation-20260817`, HEAD `193350832e514c812db2bb0eb909eb0f6602ca7d`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-17/openqp-sa-casscf-analytic-gradient-claude

- 원본 `local-15`, branch `refs/heads/claude/sa-casscf-analytic-gradient-20260817`, HEAD `de565b60a805e8875734ae7a27da02f59712bb49`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/Codex/2026-08-17/openqp-sa-casscf-zvector-pr-codex

- 원본 `local-15`, branch `refs/heads/codex/sa-casscf-zvector-pr-20260817`, HEAD `076db6cfe604572f4a1afe27186514fc50e101da`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/논문작성/.codex-worktrees/openqp-namd-no-fort6-20260827

- 원본 `local-15`, branch `refs/heads/codex/namd-no-fort6-20260827`, HEAD `7001c65c6d10cfc46926dd9b3ac982204ef76afb`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/논문작성/.codex-worktrees/openqp-xc-grid-gradient-pr-20260828

- 원본 `local-15`, branch `refs/heads/codex/xc-grid-gradient-20260828`, HEAD `a4009ae4aa8f5106b156fe284337e4c3d8a66003`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/Documents/서류작업/tmp/openqp_pr292_dual.u2JcYf

- 원본 `local-15`, branch `refs/heads/agent/clarify-current-gpl`, HEAD `8f99b59083b9cbfb8c333e25552b66bd5dd71e08`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-bethe-surface-20260917

- 원본 `local-15`, branch `refs/heads/claude/bethe-surface-scattering-20260917`, HEAD `1f047d562c6b61cc0e0f4d864e9fb8f13125f195`
- 존재: True; status 항목: 2

```text
A  pyoqp/oqp/analysis/scattering_ints.py
A  tests/scattering/test_finite_q_integrals.py
```
## /Users/cheolhochoi/openqp-worktrees/all-open-issues-20260821

- 원본 `local-15`, branch `refs/heads/claude/all-open-issues-20260821`, HEAD `7bc51dbd171869e0f2f8de9f3233ab32764947a6`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/caspt2-ci-phase-20260822

- 원본 `local-15`, branch `refs/heads/claude/caspt2-ci-phase-ff`, HEAD `96e0cc0707d10f39108eb1184fea54c8fad9381a`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/claude-berry-20260919

- 원본 `local-15`, branch `detached`, HEAD `56ebbd39183861d1fc8d5bce12f9f06990d4ab41`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/claude-static-nac-1f508e4f9

- 원본 `local-15`, branch `detached`, HEAD `1f508e4f9f2ebe3880cfd18a2e6451ff4ac9c8b7`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/espf-boundary-20260821

- 원본 `local-15`, branch `refs/heads/claude/espf-boundary-swscale-20260821`, HEAD `8b8a123c18e4dd769d20346ecde32df13d512012`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/irrep-selection-20260821

- 원본 `local-15`, branch `refs/heads/claude/ci-irrep-selection-20260821`, HEAD `a532e2a6c4cf1e00aa784b03ceba7d7ba1b76782`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/mrsf-2pa-20260917

- 원본 `local-15`, branch `refs/heads/claude/mrsf-2pa-sos-20260917`, HEAD `41b3c7d1d4c5e53660c84c3959f95c500b3ccdeb`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/qmrsf-full-xc-clean-20260918-01a0aa9f

- 원본 `local-15`, branch `refs/heads/codex/qmrsf-full-xc-clean-20260918-01a0aa9f`, HEAD `7d8ea0e4a7b1ce5f1a93554d94801319f130f3a6`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/qmrsf-h4-scan-20260918-h4agent

- 원본 `local-15`, branch `refs/heads/codex/qmrsf-h4-scan-20260918-h4agent`, HEAD `31768813cde4e963d05e51d65035c7d59a679c64`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/qmrsf-paper-audit-20260918-01a0aa9f

- 원본 `local-15`, branch `refs/heads/codex/qmrsf-paper-audit-20260918-01a0aa9f`, HEAD `3c7d475e85361c945bfe33de084b1916516ef16e`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/qmrsf-thiel-20260916-01a0aa9f

- 원본 `local-15`, branch `refs/heads/codex/qmrsf-thiel-20260916-01a0aa9f`, HEAD `ebf92ef2b3806c290fb8a0b6306fd1fd9d697e90`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-worktrees/spin-cluster-20260821

- 원본 `local-15`, branch `refs/heads/claude/ci-spin-cluster-20260821`, HEAD `5e3d27d502bb27c08908e0a4464228bfc160653e`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/openqp-xray-signals-20260917

- 원본 `local-15`, branch `refs/heads/claude/mrsf-xray-signals-20260917`, HEAD `1f047d562c6b61cc0e0f4d864e9fb8f13125f195`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/worktrees/openqp-h4-original2-mac-20260812

- 원본 `local-15`, branch `refs/heads/codex/h4-original2-mac-20260812`, HEAD `79b470814a3fbadbd7cc26149ca0f773e527be71`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/worktrees/openqp-paldus-h4

- 원본 `local-15`, branch `refs/heads/codex/paldus-h4-endpoints-20260811`, HEAD `c4dd34f895a26ec1b4377960bf32cad40dd49da2`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260805_eqcheck_overleaf/openqp-uhf-grad-plan

- 원본 `local-15`, branch `detached`, HEAD `64597d4fdc1be046ac349d416f9615ec6816a79c`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260805_eqcheck_overleaf/openqp-umrsf-final

- 원본 `local-15`, branch `detached`, HEAD `50001b25457cbb8165824b287da3d52583511a61`
- 존재: True; status 항목: 0
## /Volumes/OpenQP.DB/claude/sessions/20260805_eqcheck_overleaf/openqp-umrsf-grad

- 원본 `local-15`, branch `detached`, HEAD `79e1e77bd48ab1f7ccc0e1a73b5745807c5da05b`
- 존재: True; status 항목: 0
## /Users/cheolhochoi/clone/pr-reviews/openqp-pr179

- 원본 `local-16`, branch `refs/heads/pr-179`, HEAD `0c3c29a29f911dbd494ba5fd5ffa0ebd44ebe255`
- 존재: True; status 항목: 2

```text
?? .hermes-pr179-review.md
?? build-review/
```
