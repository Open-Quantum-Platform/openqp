"""Per-system metadata for the SCF-convergence benchmark set.

`refs` maps an SCF reference type -> spin multiplicity exercised for that system.
`profile` selects the method x converger x guess matrix (see matrix.py / gen_inputs.py).
`geom` is relative to the benchmark root (parent of harness/).
"""

# profile -> determines functionals, convergers, guesses (resolved in matrix.py)
SYSTEMS = [
    # ---- Tier 0: sanity / no-regression ----
    dict(id="t0_h2o",     geom="geometries/t0_h2o.xyz",     charge=0, refs={"rhf": 1},            tier=0, profile="t0"),
    dict(id="t0_ch4",     geom="geometries/t0_ch4.xyz",     charge=0, refs={"rhf": 1},            tier=0, profile="t0"),
    dict(id="t0_nh3",     geom="geometries/t0_nh3.xyz",     charge=0, refs={"rhf": 1},            tier=0, profile="t0"),
    dict(id="t0_benzene", geom="geometries/t0_benzene.xyz", charge=0, refs={"rhf": 1},            tier=0, profile="t0"),
    dict(id="t0_o2",      geom="geometries/t0_o2.xyz",      charge=0, refs={"uhf": 3, "rohf": 3}, tier=0, profile="t0"),
    dict(id="t0_ch3",     geom="geometries/t0_ch3.xyz",     charge=0, refs={"uhf": 2, "rohf": 2}, tier=0, profile="t0"),

    # ---- Tier 1: routine + size ----
    dict(id="t1_glycine",  geom="geometries/t1_glycine.xyz",  charge=0,   refs={"rhf": 1}, tier=1, profile="t1small"),
    dict(id="t1_c8",       geom="geometries/t1_c8.xyz",       charge=0,   refs={"rhf": 1}, tier=1, profile="t1small"),
    dict(id="t1_c16",      geom="geometries/t1_c16.xyz",      charge=0,   refs={"rhf": 1}, tier=1, profile="t1small"),
    dict(id="t1_caffeine", geom="geometries/t1_caffeine.xyz", charge=0,   refs={"rhf": 1}, tier=1, profile="t1large"),
    dict(id="t1_porphine", geom="geometries/t1_porphine.xyz", charge=0,   refs={"rhf": 1}, tier=1, profile="t1large"),

    # ---- Tier 2: hard open-shell / small-gap ----
    dict(id="t2_cr_atom",    geom="geometries/t2_cr_atom.xyz",    charge=0, refs={"rohf": 7, "uhf": 7}, tier=2, profile="t2"),
    dict(id="t2_fe_atom",    geom="geometries/t2_fe_atom.xyz",    charge=0, refs={"rohf": 5, "uhf": 5}, tier=2, profile="t2"),
    dict(id="t2_cu_atom",    geom="geometries/t2_cu_atom.xyz",    charge=0, refs={"rohf": 2, "uhf": 2}, tier=2, profile="t2"),
    dict(id="t2_feh2o6_2",   geom="geometries/t2_feh2o6_2.xyz",   charge=2, refs={"uhf": 5},            tier=2, profile="t2"),
    dict(id="t2_feh2o6_3",   geom="geometries/t2_feh2o6_3.xyz",   charge=3, refs={"uhf": 6},            tier=2, profile="t2"),
    dict(id="t2_mnh2o6_2",   geom="geometries/t2_mnh2o6_2.xyz",   charge=2, refs={"uhf": 6},            tier=2, profile="t2"),
    dict(id="t2_n2_stretch", geom="geometries/t2_n2_stretch.xyz", charge=0, refs={"rhf": 1, "uhf": 1}, tier=2, profile="t2"),
    dict(id="t2_f2_stretch", geom="geometries/t2_f2_stretch.xyz", charge=0, refs={"rhf": 1, "uhf": 1}, tier=2, profile="t2"),

    # ---- Tier 3: pathological (DIIS diverges) ----
    dict(id="t3_cr2",     geom="geometries/t3_cr2.xyz",     charge=0,  refs={"uhf": 1},        tier=3, profile="t3"),
    dict(id="t3_mo2",     geom="geometries/t3_mo2.xyz",     charge=0,  refs={"uhf": 1},        tier=3, profile="t3"),
    dict(id="t3_cu2ac4",  geom="geometries/t3_cu2ac4.xyz",  charge=0,  refs={"uhf": 1},        tier=3, profile="t3"),
    dict(id="t3_fe2s2",   geom="geometries/t3_fe2s2.xyz",   charge=-2, refs={"uhf": 1},        tier=3, profile="t3"),
    dict(id="t3_feporph", geom="geometries/t3_feporph.xyz", charge=0,  refs={"uhf": 3},        tier=3, profile="t3"),
]

# matrix profiles: functionals, convergers, guesses per profile
# functional "" == pure Hartree-Fock (no XC); others are KS-DFT.
PROFILES = {
    "t0":      dict(functionals=["", "pbe", "b3lypv5"],
                    convergers=["cdiis", "ediis", "adiis", "vdiis", "soscf", "trah"],
                    guesses=["huckel", "sap"]),
    "t1small": dict(functionals=["", "b3lypv5"],
                    convergers=["cdiis", "ediis", "adiis", "vdiis", "soscf", "trah"],
                    guesses=["huckel", "sap"]),
    "t1large": dict(functionals=["", "b3lypv5"],
                    convergers=["cdiis", "vdiis", "soscf", "trah"],
                    guesses=["huckel"]),
    "t2":      dict(functionals=["", "pbe", "b3lypv5", "bhhlyp"],
                    convergers=["cdiis", "vdiis", "adiis", "soscf", "trah"],
                    guesses=["huckel", "sap"]),
    "t3":      dict(functionals=["", "b3lypv5", "bhhlyp"],
                    convergers=["cdiis", "vdiis", "soscf", "trah"],
                    guesses=["huckel", "sap"]),
}

# global SCF settings shared by all cells
COMMON = dict(conv="1.0e-8", maxit_diis=200, maxit_trah=100, maxdiis=10,
              vdiis_cdiis_switch=0.3, vdiis_vshift_switch=0.003, vshift_vdiis=0.2)
DFTGRID = dict(rad_npts=96, ang_npts=302)
