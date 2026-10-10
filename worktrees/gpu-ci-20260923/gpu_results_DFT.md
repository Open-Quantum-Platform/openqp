# GPU DFT result (2026-07-04) — openqp-gpu vs gpu4pyscf, matched geometry

1x A100, cc-pVDZ / def2-universal-jkfit / BLYP, water_cluster(n), FP64.
Ours: spherical DF-SCF + Cartesian own-GPU XC via c2s bridge. gpu4pyscf DF-RKS.

| system   | ours E        | g4p E         | dE      | ours loop | g4p loop | ratio |
|----------|---------------|---------------|---------|-----------|----------|-------|
| (H2O)1   | -76.39793983  | (pyscf ref)   | 2.5e-11 | 12 ms     | -        | -     |
| (H2O)8   | -611.19916760 | -611.19916760 | 3.2e-10 | 959 ms    | 697 ms   | 0.73x |

DFT energies EXACT vs gpu4pyscf. Speed: ours currently 1.4x SLOWER (the opposite
of HF's 3-4.6x faster). Cause: the XC is wired with a HOST round-trip + naive
host triple-loop c2s (spherical<->cartesian) each iteration (~30 ms/iter on
200x200). The XC compute itself is on GPU and exact; only the bridge is on host.
FIX (clear): device-resident c2s (cublas) + device-pointer routec_vxc so the DFT
iteration never leaves the card, as HF does. Expected to recover the HF-class
advantage. XC validated exact standalone (SVWN/BLYP to 1e-14).

## Update: device-resident bridge + hybrid functional (2026-07-04, later)

| (H2O)8, 1xA100     | ours          | g4p           | wall  | per-iter |
|--------------------|---------------|---------------|-------|----------|
| BLYP (pure GGA)    | 565 ms (12it) | 697 ms (10it) | 1.2x  | 1.5x     |
| BHHLYP (50% exx)   | 473 ms (10it) | 854 ms (8it)  | 1.8x  | 2.3x     |

Energies identical (1e-10) in all cases. Pattern confirms the mechanism: our
advantage lives in exact exchange (occ-K). Adding 50% exx costs g4p +53% wall;
costs us ~nothing (occ-K nearly free) - we got FASTER (fewer iterations).
Remaining pure-GGA lever: grid/AO screening in the XC (ours evaluates the full
grid densely; g4p screens) - bounded known work, hybrids are the production case.
